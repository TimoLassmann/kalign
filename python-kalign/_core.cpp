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

// File-based alignment function — returns (names, sequences)
std::pair<std::vector<std::string>, std::vector<std::string>> align_from_file(
    const std::string& input_file,
    int seq_type = KALIGN_TYPE_UNDEFINED,
    float gap_open = -1.0f,
    float gap_extend = -1.0f,
    float terminal_gap_extend = -1.0f,
    int n_threads = 1,
    int refine = KALIGN_REFINE_NONE,
    int adaptive_budget = 0
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
    result = kalign_run(msa_data, n_threads, seq_type, gap_open, gap_extend, terminal_gap_extend, refine, adaptive_budget);
    if (result != 0) {
        kalign_free_msa(msa_data);
        throw std::runtime_error("Kalign alignment failed with error code: " + std::to_string(result));
    }

    // Write to a temporary file so kalign_write_msa handles gap insertion
    const char* tmpdir = std::getenv("TMPDIR");
    if (!tmpdir) tmpdir = std::getenv("TMP");
    if (!tmpdir) tmpdir = std::getenv("TEMP");
    if (!tmpdir) tmpdir = "/tmp";
    std::string temp_file = std::string(tmpdir) + "/kalign_output.fa";
    result = kalign_write_msa(msa_data, const_cast<char*>(temp_file.c_str()), const_cast<char*>("fasta"));

    kalign_free_msa(msa_data);

    if (result != 0) {
        throw std::runtime_error("Failed to write alignment results");
    }

    // Parse the FASTA file, capturing both headers (names) and sequences
    std::ifstream file(temp_file);
    std::vector<std::string> names;
    std::vector<std::string> aligned_sequences;
    std::string line, current_name, current_seq;

    while (std::getline(file, line)) {
        if (line.empty()) continue;

        if (line[0] == '>') {
            if (!current_seq.empty()) {
                names.push_back(current_name);
                aligned_sequences.push_back(current_seq);
                current_seq.clear();
            }
            // Strip the '>' prefix; take everything up to the first whitespace as the name
            current_name = line.substr(1);
            auto ws = current_name.find_first_of(" \t");
            if (ws != std::string::npos) {
                current_name = current_name.substr(0, ws);
            }
        } else {
            current_seq += line;
        }
    }

    if (!current_seq.empty()) {
        names.push_back(current_name);
        aligned_sequences.push_back(current_seq);
    }

    // Clean up temp file
    std::remove(temp_file.c_str());

    return {names, aligned_sequences};
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

// Compare two MSA files (reference vs test), returning SP score
float compare_msa_files(const std::string& reference_file, const std::string& test_file) {
    struct msa* ref = nullptr;
    struct msa* test = nullptr;
    float score = 0.0f;

    int result = kalign_read_input(const_cast<char*>(reference_file.c_str()), &ref, 1);
    if (result != 0 || !ref) {
        throw std::runtime_error("Failed to read reference file: " + reference_file);
    }

    result = kalign_read_input(const_cast<char*>(test_file.c_str()), &test, 1);
    if (result != 0 || !test) {
        kalign_free_msa(ref);
        throw std::runtime_error("Failed to read test file: " + test_file);
    }

    result = kalign_msa_compare(ref, test, &score);
    kalign_free_msa(ref);
    kalign_free_msa(test);

    if (result != 0) {
        throw std::runtime_error("MSA comparison failed");
    }

    return score;
}

// Align sequences from input file and write result to output file, preserving all metadata
void align_file_to_file(
    const std::string& input_file,
    const std::string& output_file,
    const std::string& format = "fasta",
    int seq_type = KALIGN_TYPE_UNDEFINED,
    float gap_open = -1.0f,
    float gap_extend = -1.0f,
    float terminal_gap_extend = -1.0f,
    int n_threads = 1,
    int refine = KALIGN_REFINE_NONE,
    int adaptive_budget = 0
) {
    struct msa* msa_data = nullptr;

    int result = kalign_read_input(const_cast<char*>(input_file.c_str()), &msa_data, 1);
    if (result != 0 || !msa_data) {
        throw std::runtime_error("Failed to read input file: " + input_file);
    }

    result = kalign_run(msa_data, n_threads, seq_type, gap_open, gap_extend, terminal_gap_extend, refine, adaptive_budget);
    if (result != 0) {
        kalign_free_msa(msa_data);
        throw std::runtime_error("Alignment failed with error code: " + std::to_string(result));
    }

    result = kalign_write_msa(msa_data, const_cast<char*>(output_file.c_str()), const_cast<char*>(format.c_str()));
    kalign_free_msa(msa_data);

    if (result != 0) {
        throw std::runtime_error("Failed to write output file: " + output_file);
    }
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
    
    // File-based alignment — returns (names, sequences)
    m.def("align_from_file", &align_from_file,
          py::arg("input_file"),
          py::arg("seq_type") = KALIGN_TYPE_UNDEFINED,
          py::arg("gap_open") = -1.0f,
          py::arg("gap_extend") = -1.0f,
          py::arg("terminal_gap_extend") = -1.0f,
          py::arg("n_threads") = 1,
          py::arg("refine") = KALIGN_REFINE_NONE,
          py::arg("adaptive_budget") = 0,
          "Align sequences from a file. Returns (names, sequences) tuple.");
    
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
    
    // Compare two MSA files
    m.def("compare", &compare_msa_files,
          py::arg("reference_file"),
          py::arg("test_file"),
          R"pbdoc(
          Compare two multiple sequence alignments and return SP score.

          Parameters
          ----------
          reference_file : str
              Path to reference alignment file
          test_file : str
              Path to test alignment file

          Returns
          -------
          float
              SP score (0-100)
          )pbdoc");

    // Align file to file (preserves sequence names/metadata)
    m.def("align_file_to_file", &align_file_to_file,
          py::arg("input_file"),
          py::arg("output_file"),
          py::arg("format") = "fasta",
          py::arg("seq_type") = KALIGN_TYPE_UNDEFINED,
          py::arg("gap_open") = -1.0f,
          py::arg("gap_extend") = -1.0f,
          py::arg("terminal_gap_extend") = -1.0f,
          py::arg("n_threads") = 1,
          py::arg("refine") = KALIGN_REFINE_NONE,
          py::arg("adaptive_budget") = 0,
          R"pbdoc(
          Align sequences from input file and write to output file.

          Unlike align_from_file, this preserves all sequence metadata
          (names, descriptions) which is required for MSA comparison.

          Parameters
          ----------
          input_file : str
              Path to input sequence file
          output_file : str
              Path to output alignment file
          format : str, optional
              Output format: "fasta", "msf", "clu" (default: "fasta")
          seq_type : int, optional
              Sequence type (default: auto-detect)
          gap_open : float, optional
              Gap opening penalty
          gap_extend : float, optional
              Gap extension penalty
          terminal_gap_extend : float, optional
              Terminal gap extension penalty
          n_threads : int, optional
              Number of threads (default: 1)
          )pbdoc");

    // Constants for sequence types
    m.attr("DNA") = KALIGN_TYPE_DNA;
    m.attr("DNA_INTERNAL") = KALIGN_TYPE_DNA_INTERNAL;
    m.attr("RNA") = KALIGN_TYPE_RNA;
    m.attr("PROTEIN") = KALIGN_TYPE_PROTEIN;
    m.attr("PROTEIN_DIVERGENT") = KALIGN_TYPE_PROTEIN_DIVERGENT;
    m.attr("AUTO") = KALIGN_TYPE_UNDEFINED;

    // Constants for refinement modes
    m.attr("REFINE_NONE") = KALIGN_REFINE_NONE;
    m.attr("REFINE_ALL") = KALIGN_REFINE_ALL;
    m.attr("REFINE_CONFIDENT") = KALIGN_REFINE_CONFIDENT;
}