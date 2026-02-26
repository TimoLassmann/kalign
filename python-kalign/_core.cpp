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
    #include "msa_op.h"
    #include "msa_cmp.h"
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

// Helper to extract confidence data from MSA before freeing
static py::object extract_confidence(struct msa* msa_data, int numseq) {
    if (!msa_data->col_confidence) {
        return py::none();
    }

    int alnlen = msa_data->alnlen;

    // Per-column confidence
    py::list col_conf;
    for (int c = 0; c < alnlen; c++) {
        col_conf.append(msa_data->col_confidence[c]);
    }

    // Per-residue confidence
    py::list res_conf;
    for (int i = 0; i < numseq; i++) {
        py::list row;
        if (msa_data->sequences[i]->confidence) {
            for (int c = 0; c < alnlen; c++) {
                row.append(msa_data->sequences[i]->confidence[c]);
            }
        }
        res_conf.append(row);
    }

    py::dict result;
    result["column_confidence"] = col_conf;
    result["residue_confidence"] = res_conf;
    return result;
}

// Shared alignment routing helper — selects the appropriate kalign function
// based on the combination of parameters provided.
static int run_alignment(struct msa* msa_data, int n_threads, int seq_type,
                         float gap_open, float gap_extend, float terminal_gap_extend,
                         int refine, int adaptive_budget,
                         int ensemble, uint64_t ensemble_seed,
                         float dist_scale, float vsm_amax,
                         int min_support, int realign,
                         const std::string& save_poar,
                         const std::string& load_poar,
                         float use_seq_weights = -1.0f,
                         int consistency_anchors = 0, float consistency_weight = 2.0f)
{
    if (!load_poar.empty()) {
        return kalign_consensus_from_poar(msa_data, load_poar.c_str(),
                                          min_support > 0 ? min_support : 2);
    } else if (ensemble > 0) {
        const char* save_path = save_poar.empty() ? nullptr : save_poar.c_str();
        return kalign_ensemble(msa_data, n_threads, seq_type, ensemble,
                               gap_open, gap_extend, terminal_gap_extend,
                               ensemble_seed, min_support, save_path,
                               refine, dist_scale, vsm_amax, realign,
                               use_seq_weights,
                               consistency_anchors, consistency_weight);
    } else if (realign > 0) {
        return kalign_run_realign(msa_data, n_threads, seq_type, gap_open, gap_extend, terminal_gap_extend, refine, adaptive_budget, dist_scale, vsm_amax, realign, use_seq_weights, consistency_anchors, consistency_weight);
    } else if (consistency_anchors > 0) {
        return kalign_run_seeded(msa_data, n_threads, seq_type, gap_open, gap_extend, terminal_gap_extend, refine, adaptive_budget, 0, 0.0f, dist_scale, vsm_amax, use_seq_weights, consistency_anchors, consistency_weight);
    } else if (dist_scale > 0.0f || vsm_amax >= 0.0f || use_seq_weights >= 0.0f) {
        return kalign_run_dist_scale(msa_data, n_threads, seq_type, gap_open, gap_extend, terminal_gap_extend, refine, adaptive_budget, dist_scale, vsm_amax, use_seq_weights);
    } else {
        return kalign_run(msa_data, n_threads, seq_type, gap_open, gap_extend, terminal_gap_extend, refine, adaptive_budget);
    }
}

// Main alignment function
py::object align_sequences(
    const std::vector<std::string>& sequences,
    int seq_type = KALIGN_TYPE_UNDEFINED,
    float gap_open = -1.0f,
    float gap_extend = -1.0f,
    float terminal_gap_extend = -1.0f,
    int n_threads = 1,
    int refine = KALIGN_REFINE_NONE,
    int ensemble = 0,
    int min_support = 0,
    float seq_weights = -1.0f,
    int consistency_anchors = 0, float consistency_weight = 2.0f,
    float vsm_amax = -1.0f,
    int realign = 0,
    uint64_t ensemble_seed = 42
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

    if (n_threads < 1) {
        n_threads = 1;
    }

    // Build msa struct from input arrays
    struct msa* msa_data = nullptr;
    int result = kalign_arr_to_msa(seq_ptrs.data(), seq_lengths.data(),
                                    static_cast<int>(sequences.size()), &msa_data);
    if (result != 0 || !msa_data) {
        throw std::runtime_error("Failed to create MSA from input sequences");
    }
    msa_data->quiet = 1;

    // Route to appropriate alignment function
    result = run_alignment(msa_data, n_threads, seq_type,
                           gap_open, gap_extend, terminal_gap_extend,
                           refine, 0,
                           ensemble, ensemble_seed,
                           0.0f, vsm_amax,
                           min_support, realign,
                           "", "",
                           seq_weights,
                           consistency_anchors, consistency_weight);

    if (result != 0) {
        kalign_free_msa(msa_data);
        throw std::runtime_error("Kalign alignment failed with error code: " + std::to_string(result));
    }

    // Extract confidence data before converting to arrays (which frees the MSA)
    py::object confidence = extract_confidence(msa_data, static_cast<int>(sequences.size()));

    // Extract aligned sequences
    char** aligned_seqs = nullptr;
    int alignment_length = 0;
    result = kalign_msa_to_arr(msa_data, &aligned_seqs, &alignment_length);
    kalign_free_msa(msa_data);

    if (result != 0 || !aligned_seqs) {
        throw std::runtime_error("Failed to extract aligned sequences");
    }

    // Convert results back to Python
    auto seqs = c_strings_to_python(aligned_seqs, static_cast<int>(sequences.size()), alignment_length);

    // If ensemble was used and confidence data exists, return tuple
    if (ensemble > 0 && !confidence.is_none()) {
        return py::make_tuple(seqs, confidence);
    }

    // Otherwise return just sequences (backward compat)
    return py::cast(seqs);
}

// File-based alignment function — returns (names, sequences) or (names, sequences, confidence)
py::object align_from_file(
    const std::string& input_file,
    int seq_type = KALIGN_TYPE_UNDEFINED,
    float gap_open = -1.0f,
    float gap_extend = -1.0f,
    float terminal_gap_extend = -1.0f,
    int n_threads = 1,
    int refine = KALIGN_REFINE_NONE,
    int adaptive_budget = 0,
    int ensemble = 0,
    uint64_t ensemble_seed = 42,
    float dist_scale = 0.0f,
    float vsm_amax = -1.0f,
    int min_support = 0,
    int realign = 0,
    const std::string& save_poar = "",
    const std::string& load_poar = "",
    float seq_weights = -1.0f,
    int consistency_anchors = 0, float consistency_weight = 2.0f
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
    result = run_alignment(msa_data, n_threads, seq_type,
                           gap_open, gap_extend, terminal_gap_extend,
                           refine, adaptive_budget,
                           ensemble, ensemble_seed,
                           dist_scale, vsm_amax,
                           min_support, realign,
                           save_poar, load_poar,
                           seq_weights,
                           consistency_anchors, consistency_weight);
    if (result != 0) {
        kalign_free_msa(msa_data);
        throw std::runtime_error("Kalign alignment failed with error code: " + std::to_string(result));
    }

    // Extract confidence data before writing (which doesn't preserve it)
    py::object confidence = extract_confidence(msa_data, msa_data->numseq);

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

    // If confidence data exists, return 3-tuple
    if (!confidence.is_none()) {
        return py::make_tuple(names, aligned_sequences, confidence);
    }

    return py::make_tuple(names, aligned_sequences);
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

// Compare two MSA files returning detailed POAR scores
py::dict compare_detailed_files(const std::string& reference_file,
                                const std::string& test_file,
                                float max_gap_frac = 0.2f) {
    struct msa* ref = nullptr;
    struct msa* test = nullptr;
    struct poar_score score;

    int result = kalign_read_input(const_cast<char*>(reference_file.c_str()), &ref, 1);
    if (result != 0 || !ref) {
        throw std::runtime_error("Failed to read reference file: " + reference_file);
    }

    result = kalign_read_input(const_cast<char*>(test_file.c_str()), &test, 1);
    if (result != 0 || !test) {
        kalign_free_msa(ref);
        throw std::runtime_error("Failed to read test file: " + test_file);
    }

    result = kalign_msa_compare_detailed(ref, test, max_gap_frac, &score);
    kalign_free_msa(ref);
    kalign_free_msa(test);

    if (result != 0) {
        throw std::runtime_error("Detailed MSA comparison failed");
    }

    py::dict d;
    d["recall"] = score.recall;
    d["precision"] = score.precision;
    d["f1"] = score.f1;
    d["tc"] = score.tc;
    d["ref_pairs"] = score.ref_pairs;
    d["test_pairs"] = score.test_pairs;
    d["common_pairs"] = score.common;
    return d;
}

// Compare two MSA files with an explicit column mask
py::dict compare_detailed_with_mask_files(const std::string& reference_file,
                                          const std::string& test_file,
                                          py::list column_mask) {
    struct msa* ref = nullptr;
    struct msa* test = nullptr;
    struct poar_score score;

    int result = kalign_read_input(const_cast<char*>(reference_file.c_str()), &ref, 1);
    if (result != 0 || !ref) {
        throw std::runtime_error("Failed to read reference file: " + reference_file);
    }

    result = kalign_read_input(const_cast<char*>(test_file.c_str()), &test, 1);
    if (result != 0 || !test) {
        kalign_free_msa(ref);
        throw std::runtime_error("Failed to read test file: " + test_file);
    }

    // Convert py::list to int array
    int n_cols = static_cast<int>(column_mask.size());
    std::vector<int> mask(n_cols);
    for (int i = 0; i < n_cols; i++) {
        mask[i] = column_mask[i].cast<int>();
    }

    result = kalign_msa_compare_with_mask(ref, test, mask.data(), n_cols, &score);
    kalign_free_msa(ref);
    kalign_free_msa(test);

    if (result != 0) {
        throw std::runtime_error("Detailed MSA comparison with mask failed");
    }

    py::dict d;
    d["recall"] = score.recall;
    d["precision"] = score.precision;
    d["f1"] = score.f1;
    d["tc"] = score.tc;
    d["ref_pairs"] = score.ref_pairs;
    d["test_pairs"] = score.test_pairs;
    d["common_pairs"] = score.common;
    return d;
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
    int adaptive_budget = 0,
    int ensemble = 0,
    uint64_t ensemble_seed = 42,
    float dist_scale = 0.0f,
    float vsm_amax = -1.0f,
    int min_support = 0,
    int realign = 0,
    const std::string& save_poar = "",
    const std::string& load_poar = "",
    float seq_weights = -1.0f,
    int consistency_anchors = 0, float consistency_weight = 2.0f
) {
    struct msa* msa_data = nullptr;

    int result = kalign_read_input(const_cast<char*>(input_file.c_str()), &msa_data, 1);
    if (result != 0 || !msa_data) {
        throw std::runtime_error("Failed to read input file: " + input_file);
    }

    result = run_alignment(msa_data, n_threads, seq_type,
                           gap_open, gap_extend, terminal_gap_extend,
                           refine, adaptive_budget,
                           ensemble, ensemble_seed,
                           dist_scale, vsm_amax,
                           min_support, realign,
                           save_poar, load_poar,
                           seq_weights,
                           consistency_anchors, consistency_weight);
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
          py::arg("refine") = KALIGN_REFINE_NONE,
          py::arg("ensemble") = 0,
          py::arg("min_support") = 0,
          py::arg("seq_weights") = -1.0f,
          py::arg("consistency_anchors") = 0,
          py::arg("consistency_weight") = 2.0f,
          py::arg("vsm_amax") = -1.0f,
          py::arg("realign") = 0,
          py::arg("ensemble_seed") = (uint64_t)42,
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
          refine : int, optional
              Refinement mode (default: REFINE_NONE)
          ensemble : int, optional
              Number of ensemble runs (default: 0 = off)
          min_support : int, optional
              Explicit consensus threshold (default: 0 = auto)
          vsm_amax : float, optional
              Variable scoring matrix amplitude (default: -1.0, uses Kalign defaults)
          realign : int, optional
              Number of realignment iterations (default: 0 = off)
          ensemble_seed : int, optional
              RNG seed for ensemble (default: 42)

          Returns
          -------
          list of str or tuple
              When ensemble > 0: (aligned_seqs, confidence_dict)
              Otherwise: aligned sequences
          )pbdoc");
    
    // File-based alignment — returns (names, sequences) or (names, sequences, confidence)
    m.def("align_from_file", &align_from_file,
          py::arg("input_file"),
          py::arg("seq_type") = KALIGN_TYPE_UNDEFINED,
          py::arg("gap_open") = -1.0f,
          py::arg("gap_extend") = -1.0f,
          py::arg("terminal_gap_extend") = -1.0f,
          py::arg("n_threads") = 1,
          py::arg("refine") = KALIGN_REFINE_NONE,
          py::arg("adaptive_budget") = 0,
          py::arg("ensemble") = 0,
          py::arg("ensemble_seed") = (uint64_t)42,
          py::arg("dist_scale") = 0.0f,
          py::arg("vsm_amax") = -1.0f,
          py::arg("min_support") = 0,
          py::arg("realign") = 0,
          py::arg("save_poar") = "",
          py::arg("load_poar") = "",
          py::arg("seq_weights") = -1.0f,
          py::arg("consistency_anchors") = 0,
          py::arg("consistency_weight") = 2.0f,
          "Align sequences from a file. Returns (names, sequences) or (names, sequences, confidence) tuple.");

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

    // Detailed MSA comparison (POAR recall/precision/F1/TC)
    m.def("compare_detailed", &compare_detailed_files,
          py::arg("reference_file"),
          py::arg("test_file"),
          py::arg("max_gap_frac") = 0.2f,
          R"pbdoc(
          Compare two MSAs returning detailed POAR scores.

          Parameters
          ----------
          reference_file : str
              Path to reference alignment file
          test_file : str
              Path to test alignment file
          max_gap_frac : float, optional
              Max gap fraction for scored columns (default: 0.2 for bali_score compat).
              Use -1.0 to score all columns.

          Returns
          -------
          dict
              Keys: recall, precision, f1, tc, ref_pairs, test_pairs, common_pairs
          )pbdoc");

    // Detailed MSA comparison with explicit column mask
    m.def("compare_detailed_with_mask", &compare_detailed_with_mask_files,
          py::arg("reference_file"),
          py::arg("test_file"),
          py::arg("column_mask"),
          R"pbdoc(
          Compare two MSAs with an explicit column mask.

          Parameters
          ----------
          reference_file : str
              Path to reference alignment file
          test_file : str
              Path to test alignment file
          column_mask : list of int
              Binary mask (0/1) for each column in the reference alignment.
              Only columns with mask=1 are scored for recall/TC.

          Returns
          -------
          dict
              Keys: recall, precision, f1, tc, ref_pairs, test_pairs, common_pairs
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
          py::arg("ensemble") = 0,
          py::arg("ensemble_seed") = (uint64_t)42,
          py::arg("dist_scale") = 0.0f,
          py::arg("vsm_amax") = -1.0f,
          py::arg("min_support") = 0,
          py::arg("realign") = 0,
          py::arg("save_poar") = "",
          py::arg("load_poar") = "",
          py::arg("seq_weights") = -1.0f,
          py::arg("consistency_anchors") = 0,
          py::arg("consistency_weight") = 2.0f,
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
    m.attr("PROTEIN_PFASUM43") = KALIGN_TYPE_PROTEIN_PFASUM43;
    m.attr("PROTEIN_PFASUM60") = KALIGN_TYPE_PROTEIN_PFASUM60;
    m.attr("PROTEIN_PFASUM_AUTO") = KALIGN_TYPE_PROTEIN_PFASUM_AUTO;
    m.attr("PROTEIN_DIVERGENT") = KALIGN_TYPE_PROTEIN_DIVERGENT;
    m.attr("AUTO") = KALIGN_TYPE_UNDEFINED;

    // Constants for refinement modes
    m.attr("REFINE_NONE") = KALIGN_REFINE_NONE;
    m.attr("REFINE_ALL") = KALIGN_REFINE_ALL;
    m.attr("REFINE_CONFIDENT") = KALIGN_REFINE_CONFIDENT;
    m.attr("REFINE_INLINE") = KALIGN_REFINE_INLINE;
}