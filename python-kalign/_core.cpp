#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <pybind11/numpy.h>
#include <kalign/kalign.h>
#include <kalign/kalign_config.h>
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

// Ensemble with per-run parameters — playground for optimization.
// Each run gets its own gap penalties, matrix type, and tree noise.
// Now uses kalign_align_full with per-run configs.
void ensemble_custom_file_to_file(
    const std::string& input_file,
    const std::string& output_file,
    const std::vector<float>& run_gpo,
    const std::vector<float>& run_gpe,
    const std::vector<float>& run_tgpe,
    const std::vector<float>& run_noise,
    const std::vector<int>& run_types = {},
    const std::string& format = "fasta",
    int seq_type = KALIGN_TYPE_PROTEIN,
    uint64_t seed = 42,
    int min_support = 0,
    int refine = KALIGN_REFINE_NONE,
    float vsm_amax = -1.0f,
    int realign = 0,
    float seq_weights = -1.0f,
    int n_threads = 1,
    int consistency_anchors = 0,
    float consistency_weight = 2.0f,
    // Per-run overrides: when non-empty, override the shared value per-run.
    // Same pattern as run_types: empty = use shared value for all runs.
    const std::vector<float>& run_vsm_amax = {},
    const std::vector<float>& run_seq_weights = {},
    const std::vector<int>& run_refine = {},
    const std::vector<int>& run_realign = {},
    const std::vector<int>& run_consistency_anchors = {},
    const std::vector<float>& run_consistency_weight = {},
    int consistency_merge = 0,
    float consistency_merge_weight = 2.0f
) {
    int n_runs = static_cast<int>(run_gpo.size());
    if (n_runs < 1) {
        throw std::invalid_argument("Must provide at least 1 run");
    }
    if (static_cast<int>(run_gpe.size()) != n_runs ||
        static_cast<int>(run_tgpe.size()) != n_runs ||
        static_cast<int>(run_noise.size()) != n_runs) {
        throw std::invalid_argument("All per-run arrays must have the same length");
    }
    // Validate optional per-run arrays: must be empty or same length as run_gpo
    if (!run_types.empty() && static_cast<int>(run_types.size()) != n_runs) {
        throw std::invalid_argument("run_types must be empty or same length as run_gpo");
    }
    if (!run_vsm_amax.empty() && static_cast<int>(run_vsm_amax.size()) != n_runs) {
        throw std::invalid_argument("run_vsm_amax must be empty or same length as run_gpo");
    }
    if (!run_seq_weights.empty() && static_cast<int>(run_seq_weights.size()) != n_runs) {
        throw std::invalid_argument("run_seq_weights must be empty or same length as run_gpo");
    }
    if (!run_refine.empty() && static_cast<int>(run_refine.size()) != n_runs) {
        throw std::invalid_argument("run_refine must be empty or same length as run_gpo");
    }
    if (!run_realign.empty() && static_cast<int>(run_realign.size()) != n_runs) {
        throw std::invalid_argument("run_realign must be empty or same length as run_gpo");
    }
    if (!run_consistency_anchors.empty() && static_cast<int>(run_consistency_anchors.size()) != n_runs) {
        throw std::invalid_argument("run_consistency_anchors must be empty or same length as run_gpo");
    }
    if (!run_consistency_weight.empty() && static_cast<int>(run_consistency_weight.size()) != n_runs) {
        throw std::invalid_argument("run_consistency_weight must be empty or same length as run_gpo");
    }

    struct msa* msa_data = nullptr;
    int result = kalign_read_input(const_cast<char*>(input_file.c_str()), &msa_data, 1);
    if (result != 0 || !msa_data) {
        throw std::runtime_error("Failed to read input file: " + input_file);
    }

    /* Build per-run configs.  Each optional per-run array overrides the
       shared scalar when non-empty, following the run_types pattern. */
    std::vector<kalign_run_config> runs(n_runs);
    for (int k = 0; k < n_runs; k++) {
        runs[k] = kalign_run_config_defaults();
        runs[k].matrix = (!run_types.empty()) ? run_types[k] : seq_type;
        runs[k].gpo = run_gpo[k];
        runs[k].gpe = run_gpe[k];
        runs[k].tgpe = run_tgpe[k];
        runs[k].tree_seed = seed + static_cast<uint64_t>(k);
        runs[k].tree_noise = run_noise[k];
        runs[k].vsm_amax = (!run_vsm_amax.empty()) ? run_vsm_amax[k] : vsm_amax;
        runs[k].dist_scale = 0.0f;
        runs[k].seq_weights = (!run_seq_weights.empty()) ? run_seq_weights[k] : seq_weights;
        runs[k].refine = (!run_refine.empty()) ? run_refine[k] : refine;
        runs[k].realign = (!run_realign.empty()) ? run_realign[k] : realign;
        runs[k].consistency_anchors = (!run_consistency_anchors.empty()) ? run_consistency_anchors[k] : consistency_anchors;
        runs[k].consistency_weight = (!run_consistency_weight.empty()) ? run_consistency_weight[k] : consistency_weight;
    }

    struct kalign_ensemble_config ens = kalign_ensemble_config_defaults();
    ens.min_support = min_support;
    ens.consistency_merge = consistency_merge;
    ens.consistency_merge_weight = consistency_merge_weight;

    result = kalign_align_full(msa_data, runs.data(), n_runs, &ens, n_threads);
    if (result != 0) {
        kalign_free_msa(msa_data);
        throw std::runtime_error("Ensemble alignment failed with error code: " + std::to_string(result));
    }

    result = kalign_write_msa(msa_data, const_cast<char*>(output_file.c_str()),
                               const_cast<char*>(format.c_str()));
    kalign_free_msa(msa_data);

    if (result != 0) {
        throw std::runtime_error("Failed to write output file: " + output_file);
    }
}

// Align in-memory sequences using a named mode preset.
// Detects biotype from sequences and delegates to C preset system.
py::object align_mode(
    const std::vector<std::string>& sequences,
    const std::string& mode,
    int seq_type = KALIGN_TYPE_UNDEFINED,
    float gap_open = -1.0f,
    float gap_extend = -1.0f,
    float terminal_gap_extend = -1.0f,
    int n_threads = 1
) {
    if (sequences.empty()) {
        throw std::invalid_argument("Empty sequence list provided");
    }

    std::vector<char*> seq_ptrs;
    std::vector<int> seq_lengths;
    seq_ptrs.reserve(sequences.size());
    seq_lengths.reserve(sequences.size());
    for (const auto& seq : sequences) {
        seq_ptrs.push_back(const_cast<char*>(seq.c_str()));
        seq_lengths.push_back(static_cast<int>(seq.length()));
    }
    if (n_threads < 1) n_threads = 1;

    struct msa* msa_data = nullptr;
    int result = kalign_arr_to_msa(seq_ptrs.data(), seq_lengths.data(),
                                    static_cast<int>(sequences.size()), &msa_data);
    if (result != 0 || !msa_data) {
        throw std::runtime_error("Failed to create MSA from input sequences");
    }
    msa_data->quiet = 1;

    /* Force biotype if caller specified one */
    if (seq_type == KALIGN_MATRIX_DNA || seq_type == KALIGN_MATRIX_DNA_INTERNAL ||
        seq_type == KALIGN_MATRIX_NUC_1PAM || seq_type == KALIGN_MATRIX_NUC_20PAM ||
        seq_type == KALIGN_MATRIX_NUC_200PAM) {
        msa_data->biotype = ALN_BIOTYPE_DNA;
    } else if (seq_type == KALIGN_MATRIX_RNA) {
        msa_data->biotype = ALN_BIOTYPE_DNA;  /* RNA uses DNA biotype internally */
    } else if (seq_type != KALIGN_MATRIX_AUTO && seq_type != KALIGN_TYPE_UNDEFINED) {
        msa_data->biotype = ALN_BIOTYPE_PROTEIN;
    }

    /* Detect biotype if not already set */
    if (msa_data->biotype == ALN_BIOTYPE_UNDEF) {
        result = detect_alphabet(msa_data);
        if (result != 0) {
            kalign_free_msa(msa_data);
            throw std::runtime_error("Failed to detect sequence type");
        }
    }

    /* Get preset configs */
    struct kalign_run_config runs[KALIGN_MAX_PRESET_RUNS];
    struct kalign_ensemble_config ens;
    int n_runs = 0;

    result = kalign_get_mode_preset(mode.c_str(), msa_data->biotype,
                                     runs, &n_runs, &ens);
    if (result != 0) {
        kalign_free_msa(msa_data);
        throw std::invalid_argument("Unknown mode: " + mode);
    }

    /* Apply user gap penalty overrides to all runs */
    for (int k = 0; k < n_runs; k++) {
        if (gap_open >= 0.0f) runs[k].gpo = gap_open;
        if (gap_extend >= 0.0f) runs[k].gpe = gap_extend;
        if (terminal_gap_extend >= 0.0f) runs[k].tgpe = terminal_gap_extend;
    }

    result = kalign_align_full(msa_data, runs, n_runs,
                                n_runs > 1 ? &ens : nullptr, n_threads);
    if (result != 0) {
        kalign_free_msa(msa_data);
        throw std::runtime_error("Alignment failed with error code: " + std::to_string(result));
    }

    py::object confidence = extract_confidence(msa_data, static_cast<int>(sequences.size()));

    char** aligned_seqs = nullptr;
    int alignment_length = 0;
    result = kalign_msa_to_arr(msa_data, &aligned_seqs, &alignment_length);
    kalign_free_msa(msa_data);

    if (result != 0 || !aligned_seqs) {
        throw std::runtime_error("Failed to extract aligned sequences");
    }

    auto seqs = c_strings_to_python(aligned_seqs, static_cast<int>(sequences.size()), alignment_length);

    if (n_runs > 1 && !confidence.is_none()) {
        return py::make_tuple(seqs, confidence);
    }
    return py::cast(seqs);
}

// Align from file using a named mode preset, returning (names, sequences).
py::object align_from_file_mode(
    const std::string& input_file,
    const std::string& mode,
    int seq_type = KALIGN_TYPE_UNDEFINED,
    float gap_open = -1.0f,
    float gap_extend = -1.0f,
    float terminal_gap_extend = -1.0f,
    int n_threads = 1
) {
    struct msa* msa_data = nullptr;
    int result = kalign_read_input(const_cast<char*>(input_file.c_str()), &msa_data, 1);
    if (result != 0 || !msa_data) {
        throw std::runtime_error("Failed to read input file: " + input_file);
    }

    /* Force biotype if caller specified one */
    if (seq_type == KALIGN_MATRIX_DNA || seq_type == KALIGN_MATRIX_DNA_INTERNAL) {
        msa_data->biotype = ALN_BIOTYPE_DNA;
    } else if (seq_type == KALIGN_MATRIX_RNA) {
        msa_data->biotype = ALN_BIOTYPE_DNA;
    } else if (seq_type != KALIGN_MATRIX_AUTO && seq_type != KALIGN_TYPE_UNDEFINED) {
        msa_data->biotype = ALN_BIOTYPE_PROTEIN;
    }

    if (msa_data->biotype == ALN_BIOTYPE_UNDEF) {
        result = detect_alphabet(msa_data);
        if (result != 0) {
            kalign_free_msa(msa_data);
            throw std::runtime_error("Failed to detect sequence type");
        }
    }

    struct kalign_run_config runs[KALIGN_MAX_PRESET_RUNS];
    struct kalign_ensemble_config ens;
    int n_runs = 0;

    result = kalign_get_mode_preset(mode.c_str(), msa_data->biotype,
                                     runs, &n_runs, &ens);
    if (result != 0) {
        kalign_free_msa(msa_data);
        throw std::invalid_argument("Unknown mode: " + mode);
    }

    for (int k = 0; k < n_runs; k++) {
        if (gap_open >= 0.0f) runs[k].gpo = gap_open;
        if (gap_extend >= 0.0f) runs[k].gpe = gap_extend;
        if (terminal_gap_extend >= 0.0f) runs[k].tgpe = terminal_gap_extend;
    }

    result = kalign_align_full(msa_data, runs, n_runs,
                                n_runs > 1 ? &ens : nullptr, n_threads);
    if (result != 0) {
        kalign_free_msa(msa_data);
        throw std::runtime_error("Alignment failed with error code: " + std::to_string(result));
    }

    py::object confidence = extract_confidence(msa_data, msa_data->numseq);

    /* Write to temp file to get gap-inserted FASTA output */
    const char* tmpdir = std::getenv("TMPDIR");
    if (!tmpdir) tmpdir = std::getenv("TMP");
    if (!tmpdir) tmpdir = std::getenv("TEMP");
    if (!tmpdir) tmpdir = "/tmp";
    std::string temp_file = std::string(tmpdir) + "/kalign_output.fa";
    result = kalign_write_msa(msa_data, const_cast<char*>(temp_file.c_str()),
                               const_cast<char*>("fasta"));
    kalign_free_msa(msa_data);

    if (result != 0) {
        throw std::runtime_error("Failed to write alignment results");
    }

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
            current_name = line.substr(1);
            auto ws = current_name.find_first_of(" \t");
            if (ws != std::string::npos) current_name = current_name.substr(0, ws);
        } else {
            current_seq += line;
        }
    }
    if (!current_seq.empty()) {
        names.push_back(current_name);
        aligned_sequences.push_back(current_seq);
    }
    std::remove(temp_file.c_str());

    if (!confidence.is_none()) {
        return py::make_tuple(names, aligned_sequences, confidence);
    }
    return py::make_tuple(names, aligned_sequences);
}

// Align file-to-file using a named mode preset (fast/default/recall/accurate).
// The C library provides NSGA-III optimized presets per biotype.
void align_file_to_file_mode(
    const std::string& input_file,
    const std::string& output_file,
    const std::string& mode,
    const std::string& format = "fasta",
    int n_threads = 1,
    int seq_type = KALIGN_TYPE_UNDEFINED,
    float gap_open = -1.0f,
    float gap_extend = -1.0f,
    float terminal_gap_extend = -1.0f
) {
    struct msa* msa_data = nullptr;
    int result = kalign_read_input(const_cast<char*>(input_file.c_str()), &msa_data, 1);
    if (result != 0 || !msa_data) {
        throw std::runtime_error("Failed to read input file: " + input_file);
    }

    /* Force biotype if caller specified one */
    if (seq_type == KALIGN_MATRIX_DNA || seq_type == KALIGN_MATRIX_DNA_INTERNAL) {
        msa_data->biotype = ALN_BIOTYPE_DNA;
    } else if (seq_type == KALIGN_MATRIX_RNA) {
        msa_data->biotype = ALN_BIOTYPE_DNA;
    } else if (seq_type != KALIGN_MATRIX_AUTO && seq_type != KALIGN_TYPE_UNDEFINED) {
        msa_data->biotype = ALN_BIOTYPE_PROTEIN;
    }

    /* Detect biotype from sequences so we can select the right preset grid slot */
    if (msa_data->biotype == ALN_BIOTYPE_UNDEF) {
        result = detect_alphabet(msa_data);
        if (result != 0) {
            kalign_free_msa(msa_data);
            throw std::runtime_error("Failed to detect sequence type");
        }
    }

    struct kalign_run_config runs[KALIGN_MAX_PRESET_RUNS];
    struct kalign_ensemble_config ens;
    int n_runs = 0;

    result = kalign_get_mode_preset(mode.c_str(), msa_data->biotype,
                                     runs, &n_runs, &ens);
    if (result != 0) {
        kalign_free_msa(msa_data);
        throw std::invalid_argument("Unknown mode: " + mode);
    }

    /* Apply user gap penalty overrides to all runs */
    for (int k = 0; k < n_runs; k++) {
        if (gap_open >= 0.0f) runs[k].gpo = gap_open;
        if (gap_extend >= 0.0f) runs[k].gpe = gap_extend;
        if (terminal_gap_extend >= 0.0f) runs[k].tgpe = terminal_gap_extend;
    }

    result = kalign_align_full(msa_data, runs, n_runs,
                                n_runs > 1 ? &ens : nullptr, n_threads);
    if (result != 0) {
        kalign_free_msa(msa_data);
        throw std::runtime_error("Alignment failed with error code: " + std::to_string(result));
    }

    result = kalign_write_msa(msa_data, const_cast<char*>(output_file.c_str()),
                               const_cast<char*>(format.c_str()));
    kalign_free_msa(msa_data);

    if (result != 0) {
        throw std::runtime_error("Failed to write output file: " + output_file);
    }
}

PYBIND11_MODULE(_core, m) {
    m.doc() = "Python bindings for Kalign multiple sequence alignment";
    
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

    // Ensemble with per-run parameters (optimization playground)
    m.def("ensemble_custom_file_to_file", &ensemble_custom_file_to_file,
          py::arg("input_file"),
          py::arg("output_file"),
          py::arg("run_gpo"),
          py::arg("run_gpe"),
          py::arg("run_tgpe"),
          py::arg("run_noise"),
          py::arg("run_types") = std::vector<int>{},
          py::arg("format") = "fasta",
          py::arg("seq_type") = KALIGN_TYPE_PROTEIN,
          py::arg("seed") = (uint64_t)42,
          py::arg("min_support") = 0,
          py::arg("refine") = KALIGN_REFINE_NONE,
          py::arg("vsm_amax") = -1.0f,
          py::arg("realign") = 0,
          py::arg("seq_weights") = -1.0f,
          py::arg("n_threads") = 1,
          py::arg("consistency_anchors") = 0,
          py::arg("consistency_weight") = 2.0f,
          py::arg("run_vsm_amax") = std::vector<float>{},
          py::arg("run_seq_weights") = std::vector<float>{},
          py::arg("run_refine") = std::vector<int>{},
          py::arg("run_realign") = std::vector<int>{},
          py::arg("run_consistency_anchors") = std::vector<int>{},
          py::arg("run_consistency_weight") = std::vector<float>{},
          py::arg("consistency_merge") = 0,
          py::arg("consistency_merge_weight") = 2.0f,
          R"pbdoc(
          Ensemble alignment with per-run parameters.

          Each run gets its own gap penalties, matrix type, and tree noise.
          Additional parameters can optionally be varied per-run by passing
          arrays of the same length as run_gpo. When empty (default), the
          shared scalar value is used for all runs.

          Parameters
          ----------
          input_file : str
              Path to input sequence file
          output_file : str
              Path to output alignment file
          run_gpo : list of float
              Per-run gap open penalties (length = n_runs)
          run_gpe : list of float
              Per-run gap extend penalties
          run_tgpe : list of float
              Per-run terminal gap extend penalties
          run_noise : list of float
              Per-run tree noise sigma values
          run_types : list of int, optional
              Per-run matrix types. Empty = use seq_type for all runs.
          run_vsm_amax : list of float, optional
              Per-run VSM amplitude. Empty = use vsm_amax for all runs.
          run_seq_weights : list of float, optional
              Per-run profile rebalancing weight. Empty = use seq_weights for all.
          run_refine : list of int, optional
              Per-run refinement mode (REFINE_* constants). Empty = use refine for all.
          run_realign : list of int, optional
              Per-run realign iterations. Empty = use realign for all.
          run_consistency_anchors : list of int, optional
              Per-run consistency rounds. Empty = use consistency_anchors for all.
          run_consistency_weight : list of float, optional
              Per-run consistency weight. Empty = use consistency_weight for all.
          )pbdoc");

    // In-memory alignment using a named mode preset
    m.def("align", &align_mode,
          py::arg("sequences"),
          py::arg("mode"),
          py::arg("seq_type") = KALIGN_TYPE_UNDEFINED,
          py::arg("gap_open") = -1.0f,
          py::arg("gap_extend") = -1.0f,
          py::arg("terminal_gap_extend") = -1.0f,
          py::arg("n_threads") = 1,
          "Align sequences using a named mode preset (fast/default/recall/accurate).");
    m.def("align_mode", &align_mode,
          py::arg("sequences"),
          py::arg("mode"),
          py::arg("seq_type") = KALIGN_TYPE_UNDEFINED,
          py::arg("gap_open") = -1.0f,
          py::arg("gap_extend") = -1.0f,
          py::arg("terminal_gap_extend") = -1.0f,
          py::arg("n_threads") = 1,
          "Alias for align(). Align sequences using a named mode preset.");

    // File alignment returning (names, sequences) using a named mode preset
    m.def("align_from_file", &align_from_file_mode,
          py::arg("input_file"),
          py::arg("mode"),
          py::arg("seq_type") = KALIGN_TYPE_UNDEFINED,
          py::arg("gap_open") = -1.0f,
          py::arg("gap_extend") = -1.0f,
          py::arg("terminal_gap_extend") = -1.0f,
          py::arg("n_threads") = 1,
          "Align from file using a named mode preset. Returns (names, sequences) tuple.");
    m.def("align_from_file_mode", &align_from_file_mode,
          py::arg("input_file"),
          py::arg("mode"),
          py::arg("seq_type") = KALIGN_TYPE_UNDEFINED,
          py::arg("gap_open") = -1.0f,
          py::arg("gap_extend") = -1.0f,
          py::arg("terminal_gap_extend") = -1.0f,
          py::arg("n_threads") = 1,
          "Alias for align_from_file(). Align from file using a named mode preset.");

    // File-to-file alignment using a named mode preset
    m.def("align_file_to_file", &align_file_to_file_mode,
          py::arg("input_file"),
          py::arg("output_file"),
          py::arg("mode"),
          py::arg("format") = "fasta",
          py::arg("n_threads") = 1,
          py::arg("seq_type") = KALIGN_TYPE_UNDEFINED,
          py::arg("gap_open") = -1.0f,
          py::arg("gap_extend") = -1.0f,
          py::arg("terminal_gap_extend") = -1.0f,
          "Align file to file using a named mode preset (fast/default/recall/accurate).");
    m.def("align_file_to_file_mode", &align_file_to_file_mode,
          py::arg("input_file"),
          py::arg("output_file"),
          py::arg("mode"),
          py::arg("format") = "fasta",
          py::arg("n_threads") = 1,
          py::arg("seq_type") = KALIGN_TYPE_UNDEFINED,
          py::arg("gap_open") = -1.0f,
          py::arg("gap_extend") = -1.0f,
          py::arg("terminal_gap_extend") = -1.0f,
          "Alias for align_file_to_file(). Align file to file using a named mode preset.");

    // Matrix constants (canonical names)
    m.attr("MATRIX_AUTO") = KALIGN_MATRIX_AUTO;
    m.attr("MATRIX_PFASUM43") = KALIGN_MATRIX_PFASUM43;
    m.attr("MATRIX_PFASUM60") = KALIGN_MATRIX_PFASUM60;
    m.attr("MATRIX_CORBLOSUM66") = KALIGN_MATRIX_CORBLOSUM66;
    m.attr("MATRIX_DNA") = KALIGN_MATRIX_DNA;
    m.attr("MATRIX_DNA_INTERNAL") = KALIGN_MATRIX_DNA_INTERNAL;
    m.attr("MATRIX_RNA") = KALIGN_MATRIX_RNA;
    m.attr("MATRIX_NUC_1PAM") = KALIGN_MATRIX_NUC_1PAM;
    m.attr("MATRIX_NUC_20PAM") = KALIGN_MATRIX_NUC_20PAM;
    m.attr("MATRIX_NUC_200PAM") = KALIGN_MATRIX_NUC_200PAM;

    // Backward compat: old names point to new values
    m.attr("DNA") = KALIGN_TYPE_DNA;
    m.attr("DNA_INTERNAL") = KALIGN_TYPE_DNA_INTERNAL;
    m.attr("RNA") = KALIGN_TYPE_RNA;
    m.attr("PROTEIN") = KALIGN_TYPE_PROTEIN;
    m.attr("PROTEIN_PFASUM43") = KALIGN_TYPE_PROTEIN_PFASUM43;
    m.attr("PROTEIN_PFASUM60") = KALIGN_TYPE_PROTEIN_PFASUM60;
    m.attr("PROTEIN_PFASUM_AUTO") = KALIGN_TYPE_PROTEIN_PFASUM_AUTO;
    m.attr("PROTEIN_DIVERGENT") = KALIGN_TYPE_PROTEIN_DIVERGENT;
    m.attr("PROTEIN_CORBLOSUM66") = KALIGN_TYPE_PROTEIN_CORBLOSUM66;
    m.attr("AUTO") = KALIGN_TYPE_UNDEFINED;

    // Constants for refinement modes
    m.attr("REFINE_NONE") = KALIGN_REFINE_NONE;
    m.attr("REFINE_ALL") = KALIGN_REFINE_ALL;
    m.attr("REFINE_CONFIDENT") = KALIGN_REFINE_CONFIDENT;
    m.attr("REFINE_INLINE") = KALIGN_REFINE_INLINE;

}