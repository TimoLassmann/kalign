#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <pybind11/numpy.h>
#include <kalign/kalign.h>
#include <vector>
#include <string>
#include <memory>
#include <stdexcept>
#include <fstream>
#include <cstdio>
#include <cstdlib>

// Include DSSim for sequence simulation and MSA structures
extern "C" {
    #include "msa_struct.h"
    #include "msa_alloc.h"
    #include "dssim.h"
}

namespace py = pybind11;

// Helper function to convert C strings to Python strings and free memory
std::vector<std::string> c_strings_to_python(char** c_strings, int count, int length) {
    std::vector<std::string> result;
    result.reserve(count);
    
    for (int i = 0; i < count; ++i) {
        if (c_strings[i]) {
            result.emplace_back(c_strings[i], length);
            free(c_strings[i]);
        }
    }
    free(c_strings);
    return result;
}

// Main alignment function
std::vector<std::string> align_sequences(
    const std::vector<std::string>& sequences,
    int seq_type = KALIGN_TYPE_UNDEFINED,
    float gap_open = -1.0f,
    float gap_extend = -1.0f, 
    float terminal_gap_extend = -1.0f,
    int n_threads = 1
) {
    if (sequences.empty()) {
        throw std::invalid_argument("Empty sequence list provided");
    }
    
    // Convert Python strings to C format
    std::vector<char*> seq_ptrs;
    std::vector<int> seq_lengths;
    seq_ptrs.reserve(sequences.size());
    seq_lengths.reserve(sequences.size());
    
    for (const auto& seq : sequences) {
        seq_ptrs.push_back(const_cast<char*>(seq.c_str()));
        seq_lengths.push_back(static_cast<int>(seq.length()));
    }
    
    // Call kalign function
    char** aligned_seqs = nullptr;
    int alignment_length = 0;
    
    int result = kalign(
        seq_ptrs.data(),
        seq_lengths.data(),
        static_cast<int>(sequences.size()),
        n_threads,
        seq_type,
        gap_open,
        gap_extend,
        terminal_gap_extend,
        &aligned_seqs,
        &alignment_length
    );
    
    if (result != 0) {
        throw std::runtime_error("Kalign alignment failed with error code: " + std::to_string(result));
    }
    
    if (!aligned_seqs) {
        throw std::runtime_error("Kalign returned null aligned sequences");
    }
    
    // Convert results back to Python
    return c_strings_to_python(aligned_seqs, static_cast<int>(sequences.size()), alignment_length);
}

// File-based alignment function
std::vector<std::string> align_from_file(
    const std::string& input_file,
    int seq_type = KALIGN_TYPE_UNDEFINED,
    float gap_open = -1.0f,
    float gap_extend = -1.0f,
    float terminal_gap_extend = -1.0f,
    int n_threads = 1
) {
    struct msa* msa_data = nullptr;
    
    // Read input file
    int result = kalign_read_input(const_cast<char*>(input_file.c_str()), &msa_data, 1);
    if (result != 0) {
        throw std::runtime_error("Failed to read input file: " + input_file);
    }
    
    // Check if msa_data is NULL - this happens when the file format cannot be detected
    if (!msa_data) {
        throw std::runtime_error("Could not detect valid sequence format in file: " + input_file);
    }
    
    // Perform alignment
    result = kalign_run(msa_data, n_threads, seq_type, gap_open, gap_extend, terminal_gap_extend);
    if (result != 0) {
        kalign_free_msa(msa_data);
        throw std::runtime_error("Kalign alignment failed with error code: " + std::to_string(result));
    }
    
    // For now, we'll need to write to a temporary file and read back
    // This is a limitation of the current C API that doesn't expose aligned sequences directly
    std::string temp_file = "/tmp/kalign_output.fa";
    result = kalign_write_msa(msa_data, const_cast<char*>(temp_file.c_str()), const_cast<char*>("fasta"));
    
    kalign_free_msa(msa_data);
    
    if (result != 0) {
        throw std::runtime_error("Failed to write alignment results");
    }
    
    // Read the temporary file back
    // This is a simplified approach - in practice, you'd want to parse the FASTA file properly
    std::ifstream file(temp_file);
    std::vector<std::string> aligned_sequences;
    std::string line, current_seq;
    
    while (std::getline(file, line)) {
        if (line.empty()) continue;
        
        if (line[0] == '>') {
            if (!current_seq.empty()) {
                aligned_sequences.push_back(current_seq);
                current_seq.clear();
            }
        } else {
            current_seq += line;
        }
    }
    
    if (!current_seq.empty()) {
        aligned_sequences.push_back(current_seq);
    }
    
    // Clean up temp file
    std::remove(temp_file.c_str());
    
    return aligned_sequences;
}

// Write alignment to file
void write_alignment(
    const std::vector<std::string>& sequences,
    const std::string& output_file,
    const std::string& format = "fasta"
) {
    // This would need to be implemented by converting sequences back to msa struct
    // For now, we'll provide a simple implementation
    throw std::runtime_error("write_alignment not yet implemented - use align_from_file instead");
}

// Generate test sequences using DSSim
std::vector<std::string> generate_test_sequences(
    int n_seq,
    int n_obs,
    bool dna,
    int length,
    int seed = 42
) {
    struct msa* msa_data = nullptr;
    
    // Call DSSim to generate sequences
    int result = dssim_get_fasta(&msa_data, n_seq, n_obs, dna ? 1 : 0, length, seed);
    
    if (result != 0) {
        throw std::runtime_error("DSSim sequence generation failed with error code: " + std::to_string(result));
    }
    
    if (!msa_data) {
        throw std::runtime_error("DSSim returned null MSA data");
    }
    
    // Convert MSA sequences to Python strings
    std::vector<std::string> sequences;
    sequences.reserve(msa_data->numseq);
    
    for (int i = 0; i < msa_data->numseq; ++i) {
        if (msa_data->sequences[i] && msa_data->sequences[i]->seq) {
            sequences.emplace_back(msa_data->sequences[i]->seq, msa_data->sequences[i]->len);
        }
    }
    
    // Clean up
    kalign_free_msa(msa_data);
    
    return sequences;
}

PYBIND11_MODULE(_core, m) {
    m.doc() = "Python bindings for Kalign multiple sequence alignment";
    
    // Main alignment function
    m.def("align", &align_sequences,
          py::arg("sequences"),
          py::arg("seq_type") = KALIGN_TYPE_UNDEFINED,
          py::arg("gap_open") = -1.0f,
          py::arg("gap_extend") = -1.0f,
          py::arg("terminal_gap_extend") = -1.0f,
          py::arg("n_threads") = 1,
          R"pbdoc(
          Align a list of sequences using Kalign.
          
          Parameters
          ----------
          sequences : list of str
              List of sequences to align
          seq_type : int, optional
              Sequence type (default: auto-detect)
          gap_open : float, optional
              Gap opening penalty (default: -1.0, uses Kalign defaults)
          gap_extend : float, optional
              Gap extension penalty (default: -1.0, uses Kalign defaults) 
          terminal_gap_extend : float, optional
              Terminal gap extension penalty (default: -1.0, uses Kalign defaults)
          n_threads : int, optional
              Number of threads to use (default: 1)
              
          Returns
          -------
          list of str
              Aligned sequences
          )pbdoc");
    
    // File-based alignment
    m.def("align_from_file", &align_from_file,
          py::arg("input_file"),
          py::arg("seq_type") = KALIGN_TYPE_UNDEFINED,
          py::arg("gap_open") = -1.0f,
          py::arg("gap_extend") = -1.0f,
          py::arg("terminal_gap_extend") = -1.0f,
          py::arg("n_threads") = 1,
          "Align sequences from a file");
    
    // Write alignment function (placeholder)
    m.def("write_alignment", &write_alignment,
          py::arg("sequences"),
          py::arg("output_file"),
          py::arg("format") = "fasta",
          "Write aligned sequences to file");
    
    // Generate test sequences using DSSim
    m.def("generate_test_sequences", &generate_test_sequences,
          py::arg("n_seq"),
          py::arg("n_obs"),
          py::arg("dna"),
          py::arg("length"),
          py::arg("seed") = 42,
          R"pbdoc(
          Generate test sequences using DSSim HMM-based simulator.
          
          Parameters
          ----------
          n_seq : int
              Number of sequences to generate
          n_obs : int
              Number of observed sequences for training the HMM
          dna : bool
              True for DNA sequences, False for protein sequences
          length : int
              Target sequence length
          seed : int, optional
              Random seed for reproducible results (default: 42)
              
          Returns
          -------
          list of str
              Generated sequences
          )pbdoc");
    
    // Constants for sequence types
    m.attr("DNA") = KALIGN_TYPE_DNA;
    m.attr("DNA_INTERNAL") = KALIGN_TYPE_DNA_INTERNAL;
    m.attr("RNA") = KALIGN_TYPE_RNA;
    m.attr("PROTEIN") = KALIGN_TYPE_PROTEIN;
    m.attr("PROTEIN_DIVERGENT") = KALIGN_TYPE_PROTEIN_DIVERGENT;
    m.attr("AUTO") = KALIGN_TYPE_UNDEFINED;
}