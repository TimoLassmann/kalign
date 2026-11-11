#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <pybind11/numpy.h>
#include <kalign/kalign.h>
#include <kalign/kalign_config.h>
#include <vector>
#include <string>
#include <memory>
#include <stdexcept>
#include <cstdlib>

// Include DSSim for sequence simulation and MSA structures
extern "C" {
    #include "msa_struct.h"
    #include "msa_alloc.h"
    #include "msa_op.h"
    #include "msa_cmp.h"
    #include "aln_add.h"
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

// Append the C library's first error message on this thread (the root cause) to a
// failure message, and reset it for the next call.
static std::string with_reason(const std::string& what, const char* fallback = "") {
    const char* r = kalign_last_error();
    std::string msg = (r && r[0]) ? what + ": " + r : (fallback[0] ? what + ": " + fallback : what);
    kalign_clear_error();
    return msg;
}

// Compare two MSA files (reference vs test), returning SP score
float compare_msa_files(const std::string& reference_file, const std::string& test_file) {
    struct msa* ref = nullptr;
    struct msa* test = nullptr;
    float score = 0.0f;

    kalign_clear_error();
    int result = kalign_read_input(const_cast<char*>(reference_file.c_str()), &ref, 1);
    if (result != 0 || !ref) {
        throw std::runtime_error(with_reason("Failed to read reference file: " + reference_file, result == 0 ? "no sequences found" : ""));
    }

    kalign_clear_error();
    result = kalign_read_input(const_cast<char*>(test_file.c_str()), &test, 1);
    if (result != 0 || !test) {
        kalign_free_msa(ref);
        throw std::runtime_error(with_reason("Failed to read test file: " + test_file, result == 0 ? "no sequences found" : ""));
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

    kalign_clear_error();
    int result = kalign_read_input(const_cast<char*>(reference_file.c_str()), &ref, 1);
    if (result != 0 || !ref) {
        throw std::runtime_error(with_reason("Failed to read reference file: " + reference_file, result == 0 ? "no sequences found" : ""));
    }

    kalign_clear_error();
    result = kalign_read_input(const_cast<char*>(test_file.c_str()), &test, 1);
    if (result != 0 || !test) {
        kalign_free_msa(ref);
        throw std::runtime_error(with_reason("Failed to read test file: " + test_file, result == 0 ? "no sequences found" : ""));
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

    kalign_clear_error();
    int result = kalign_read_input(const_cast<char*>(reference_file.c_str()), &ref, 1);
    if (result != 0 || !ref) {
        throw std::runtime_error(with_reason("Failed to read reference file: " + reference_file, result == 0 ? "no sequences found" : ""));
    }

    kalign_clear_error();
    result = kalign_read_input(const_cast<char*>(test_file.c_str()), &test, 1);
    if (result != 0 || !test) {
        kalign_free_msa(ref);
        throw std::runtime_error(with_reason("Failed to read test file: " + test_file, result == 0 ? "no sequences found" : ""));
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

// Ensemble with per-run parameters.
// DEVELOPMENT INTERFACE: this function was used to fit the Kalign 4
// settings. It is not part of the supported API and may change. Its
// defaults are those of the low-level library, not of the fast or accurate
// mode. Use align, align_from_file or align_file_to_file with mode instead.
// Each run gets its own gap penalties, matrix type, and tree noise.
// Uses kalign_align_full with per-run configs.
py::object ensemble_custom_file_to_file(
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
    float consistency_merge_weight = 2.0f,
    int vsm_mode = 0,
    float vsm_kappa = 1.0f,
    float vsm_dzero = 1.0f,
    const std::vector<float>& run_vsm_kappa = {},
    const std::vector<float>& run_vsm_dzero = {},
    int fill_in = 0,
    float fill_anchor_frac = 0.0f,
    int fill_support = 1,
    int guide_dist = KALIGN_GUIDE_DIST_BPM,
    const std::vector<int>& run_guide_dist = {},
    const std::vector<int>& run_tree_method = {},
    const std::vector<int>& run_tree_dist = {},
    const std::vector<int>& run_tree_alpha = {},
    const std::vector<int>& run_tree_maxlen = {},
    const std::vector<int>& run_tree_leaf = {},
    const std::vector<int>& run_tree_kmeans_dist = {},
    const std::vector<int>& run_realign_tree = {},
    const std::vector<int>& run_gap_model = {},
    const std::vector<int>& run_vsm_mode = {},
    const std::vector<float>& run_tree_sl_noise = {}
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
    if (!run_vsm_kappa.empty() && static_cast<int>(run_vsm_kappa.size()) != n_runs) {
        throw std::invalid_argument("run_vsm_kappa must be empty or same length as run_gpo");
    }
    if (!run_vsm_dzero.empty() && static_cast<int>(run_vsm_dzero.size()) != n_runs) {
        throw std::invalid_argument("run_vsm_dzero must be empty or same length as run_gpo");
    }
    if (vsm_mode != 0 && vsm_mode != 1) {
        throw std::invalid_argument("vsm_mode must be 0 (legacy) or 1 (per-merge)");
    }
    if (!run_vsm_mode.empty() && static_cast<int>(run_vsm_mode.size()) != n_runs) {
        throw std::invalid_argument("run_vsm_mode must be empty or same length as run_gpo");
    }
    for (int v : run_vsm_mode) {
        if (v != 0 && v != 1) {
            throw std::invalid_argument("run_vsm_mode must be 0 (legacy) or 1 (per-merge)");
        }
    }
    if (fill_in < KALIGN_FILL_OFF || fill_in > KALIGN_FILL_ATTACH) {
        throw std::invalid_argument("fill_in must be 0 (off), 1 (realign), 2 (project) or 3 (attach)");
    }
    if (!(fill_anchor_frac >= 0.0f && fill_anchor_frac <= 1.0f)) {
        throw std::invalid_argument("fill_anchor_frac must be in [0, 1]");
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
    if (!run_guide_dist.empty() && static_cast<int>(run_guide_dist.size()) != n_runs) {
        throw std::invalid_argument("run_guide_dist must be empty or same length as run_gpo");
    }
    if (!run_tree_method.empty() && static_cast<int>(run_tree_method.size()) != n_runs) {
        throw std::invalid_argument("run_tree_method must be empty or same length as run_gpo");
    }
    if (!run_tree_dist.empty() && static_cast<int>(run_tree_dist.size()) != n_runs) {
        throw std::invalid_argument("run_tree_dist must be empty or same length as run_gpo");
    }
    if (!run_tree_alpha.empty() && static_cast<int>(run_tree_alpha.size()) != n_runs) {
        throw std::invalid_argument("run_tree_alpha must be empty or same length as run_gpo");
    }
    if (!run_tree_maxlen.empty() && static_cast<int>(run_tree_maxlen.size()) != n_runs) {
        throw std::invalid_argument("run_tree_maxlen must be empty or same length as run_gpo");
    }
    if (!run_tree_leaf.empty() && static_cast<int>(run_tree_leaf.size()) != n_runs) {
        throw std::invalid_argument("run_tree_leaf must be empty or same length as run_gpo");
    }
    if (!run_tree_kmeans_dist.empty() && static_cast<int>(run_tree_kmeans_dist.size()) != n_runs) {
        throw std::invalid_argument("run_tree_kmeans_dist must be empty or same length as run_gpo");
    }
    if (!run_realign_tree.empty() && static_cast<int>(run_realign_tree.size()) != n_runs) {
        throw std::invalid_argument("run_realign_tree must be empty or same length as run_gpo");
    }
    if (!run_gap_model.empty() && static_cast<int>(run_gap_model.size()) != n_runs) {
        throw std::invalid_argument("run_gap_model must be empty or same length as run_gpo");
    }
    for (int v : run_gap_model) {
        if (v != 0 && v != 1) {
            throw std::invalid_argument("run_gap_model must be 0 (Kalign gap rule) or 1 (member gap model)");
        }
    }
    if (!run_tree_sl_noise.empty() && static_cast<int>(run_tree_sl_noise.size()) != n_runs) {
        throw std::invalid_argument("run_tree_sl_noise must be empty or same length as run_gpo");
    }
    for (float v : run_tree_sl_noise) {
        if (!(v >= 0.0f)) {
            throw std::invalid_argument("run_tree_sl_noise must be >= 0");
        }
    }
    {
        std::vector<int> gd = run_guide_dist;
        gd.push_back(guide_dist);
        for (int v : gd) {
            if (v != KALIGN_GUIDE_DIST_BPM && v != KALIGN_GUIDE_DIST_PAIRALN) {
                throw std::invalid_argument("guide_dist must be 0 (BPM) or 1 (pairwise alignment)");
            }
        }
    }

    struct msa* msa_data = nullptr;
    kalign_clear_error();
    int result = kalign_read_input(const_cast<char*>(input_file.c_str()), &msa_data, 1);
    if (result != 0 || !msa_data) {
        throw std::runtime_error(with_reason("Failed to read input file: " + input_file, result == 0 ? "no sequences found" : ""));
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
        runs[k].vsm_mode = (!run_vsm_mode.empty()) ? run_vsm_mode[k] : vsm_mode;
        runs[k].vsm_kappa = (!run_vsm_kappa.empty()) ? run_vsm_kappa[k] : vsm_kappa;
        runs[k].vsm_dzero = (!run_vsm_dzero.empty()) ? run_vsm_dzero[k] : vsm_dzero;
        runs[k].guide_dist = (!run_guide_dist.empty()) ? run_guide_dist[k] : guide_dist;
        if (!run_tree_method.empty()) runs[k].tree_method = run_tree_method[k];
        if (!run_tree_dist.empty()) runs[k].tree_dist = run_tree_dist[k];
        if (!run_tree_alpha.empty()) runs[k].tree_alpha = run_tree_alpha[k];
        if (!run_tree_maxlen.empty()) runs[k].tree_maxlen = run_tree_maxlen[k];
        if (!run_tree_leaf.empty()) runs[k].tree_leaf = run_tree_leaf[k];
        if (!run_tree_kmeans_dist.empty()) runs[k].tree_kmeans_dist = run_tree_kmeans_dist[k];
        if (!run_realign_tree.empty()) runs[k].realign_tree = run_realign_tree[k];
        if (!run_gap_model.empty()) runs[k].gap_model = run_gap_model[k];
        if (!run_tree_sl_noise.empty()) runs[k].tree_sl_noise = run_tree_sl_noise[k];
    }

    struct kalign_ensemble_config ens = kalign_ensemble_config_defaults();
    ens.min_support = min_support;
    ens.consistency_merge = consistency_merge;
    ens.consistency_merge_weight = consistency_merge_weight;
    ens.fill_in = fill_in;
    ens.fill_anchor_frac = fill_anchor_frac;
    ens.fill_support = fill_support;

    kalign_clear_error();
    result = kalign_align_full(msa_data, runs.data(), n_runs, &ens, n_threads);
    if (result != 0) {
        kalign_free_msa(msa_data);
        throw std::runtime_error(with_reason("Ensemble alignment failed with error code: " + std::to_string(result)));
    }

    /* Confidence (None for a 1-member run) and the row names in the same order, taken before the write. */
    py::object confidence = extract_confidence(msa_data, msa_data->numseq);
    py::list names;
    for (int i = 0; i < msa_data->numseq; i++) {
        names.append(std::string(msa_data->sequences[i]->name));
    }

    kalign_clear_error();
    result = kalign_write_msa(msa_data, const_cast<char*>(output_file.c_str()),
                               const_cast<char*>(format.c_str()));
    kalign_free_msa(msa_data);

    if (result != 0) {
        throw std::runtime_error(with_reason("Failed to write output file: " + output_file));
    }
    if (confidence.is_none()) {
        return py::none();
    }
    py::dict out = confidence.cast<py::dict>();
    out["names"] = names;
    return out;
}

// Align in-memory sequences using a named mode preset.
// Sequence type (seq_type), shared by every entry point: the same names and
// the same behaviour as the C command line's --type (kalign_parse_seq_type /
// kalign_set_seq_type). "auto" keeps the detected type; "protein", "dna" and
// "rna" override it and select the protein or nucleotide presets. A
// UserWarning is issued when the override differs from the detected type.
// Does not free msa; callers free it if this throws.
static void apply_seq_type(struct msa* msa, const std::string& seq_type, const std::string& what)
{
    int st = KALIGN_SEQ_AUTO;
    if (kalign_parse_seq_type(seq_type.c_str(), &st) != 0) {
        throw std::invalid_argument("Invalid seq_type: '" + seq_type +
                                    "'. Must be one of: auto, protein, dna, rna");
    }
    std::string detected = kalign_seq_type_name(msa);
    int overridden = 0;
    if (kalign_set_seq_type(msa, st, &overridden) != 0) {
        throw std::invalid_argument("Invalid seq_type: '" + seq_type + "'");
    }
    if (overridden) {
        std::string msg = what + ": sequences look like " + detected + ", aligning as " +
                          kalign_seq_type_name(msa) + " as requested (seq_type='" + seq_type + "').";
        if (PyErr_WarnEx(PyExc_UserWarning, msg.c_str(), 2) != 0) {
            throw py::error_already_set();
        }
    }
}

static void apply_seq_type_or_free(struct msa* msa, const std::string& seq_type, const std::string& what)
{
    try {
        apply_seq_type(msa, seq_type, what);
    } catch (...) {
        kalign_free_msa(msa);
        throw;
    }
}

// Detects biotype from sequences and delegates to C preset system.
py::object align_mode(
    const std::vector<std::string>& sequences,
    const std::string& mode,
    const std::string& seq_type = "auto",
    float gap_open = -1.0f,
    float gap_extend = -1.0f,
    float terminal_gap_extend = -1.0f,
    int n_threads = 1,
    int fill_in = -1,
    float fill_anchor_frac = -1.0f,
    int fill_support = 0
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
    kalign_clear_error();
    int result = kalign_arr_to_msa(seq_ptrs.data(), seq_lengths.data(),
                                    static_cast<int>(sequences.size()), &msa_data);
    if (result != 0 || !msa_data) {
        throw std::runtime_error(with_reason("Failed to create MSA from input sequences"));
    }
    msa_data->quiet = 1;

    apply_seq_type_or_free(msa_data, seq_type, "sequences");

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
    if (!kalign_mode_size_ok(mode.c_str(), msa_data->biotype, msa_data->numseq)) {
        int n = msa_data->numseq;
        kalign_free_msa(msa_data);
        throw std::invalid_argument("Mode '" + mode + "' is not feasible above " +
                                    std::to_string(KALIGN_ACCURATE_MAX_SEQ) + " sequences (input: " +
                                    std::to_string(n) + "). Use mode='fast'.");
    }

    /* Apply user gap penalty overrides to all runs */
    for (int k = 0; k < n_runs; k++) {
        if (gap_open >= 0.0f) runs[k].gpo = gap_open;
        if (gap_extend >= 0.0f) runs[k].gpe = gap_extend;
        if (terminal_gap_extend >= 0.0f) runs[k].tgpe = terminal_gap_extend;
    }
    /* Fill-in overrides (negative = keep preset) */
    if (fill_in >= 0) ens.fill_in = fill_in;
    if (fill_anchor_frac >= 0.0f) ens.fill_anchor_frac = fill_anchor_frac;
    if (fill_support >= 1) ens.fill_support = fill_support;

    kalign_clear_error();
    result = kalign_align_full(msa_data, runs, n_runs,
                                n_runs > 1 ? &ens : nullptr, n_threads);
    if (result != 0) {
        kalign_free_msa(msa_data);
        throw std::runtime_error(with_reason("Alignment failed with error code: " + std::to_string(result)));
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
    const std::string& seq_type = "auto",
    float gap_open = -1.0f,
    float gap_extend = -1.0f,
    float terminal_gap_extend = -1.0f,
    int n_threads = 1,
    int fill_in = -1,
    float fill_anchor_frac = -1.0f,
    int fill_support = 0
) {
    struct msa* msa_data = nullptr;
    kalign_clear_error();
    int result = kalign_read_input(const_cast<char*>(input_file.c_str()), &msa_data, 1);
    if (result != 0 || !msa_data) {
        throw std::runtime_error(with_reason("Failed to read input file: " + input_file, result == 0 ? "no sequences found" : ""));
    }

    apply_seq_type_or_free(msa_data, seq_type, input_file);

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
    if (!kalign_mode_size_ok(mode.c_str(), msa_data->biotype, msa_data->numseq)) {
        int n = msa_data->numseq;
        kalign_free_msa(msa_data);
        throw std::invalid_argument("Mode '" + mode + "' is not feasible above " +
                                    std::to_string(KALIGN_ACCURATE_MAX_SEQ) + " sequences (input: " +
                                    std::to_string(n) + "). Use mode='fast'.");
    }

    for (int k = 0; k < n_runs; k++) {
        if (gap_open >= 0.0f) runs[k].gpo = gap_open;
        if (gap_extend >= 0.0f) runs[k].gpe = gap_extend;
        if (terminal_gap_extend >= 0.0f) runs[k].tgpe = terminal_gap_extend;
    }
    /* Fill-in overrides (negative = keep preset) */
    if (fill_in >= 0) ens.fill_in = fill_in;
    if (fill_anchor_frac >= 0.0f) ens.fill_anchor_frac = fill_anchor_frac;
    if (fill_support >= 1) ens.fill_support = fill_support;

    kalign_clear_error();
    result = kalign_align_full(msa_data, runs, n_runs,
                                n_runs > 1 ? &ens : nullptr, n_threads);
    if (result != 0) {
        kalign_free_msa(msa_data);
        throw std::runtime_error(with_reason("Alignment failed with error code: " + std::to_string(result)));
    }

    py::object confidence = extract_confidence(msa_data, msa_data->numseq);

    int numseq = msa_data->numseq;
    int alnlen = msa_data->alnlen;

    std::vector<std::string> names;
    std::vector<std::string> aligned_sequences;
    names.reserve(numseq);
    aligned_sequences.reserve(numseq);

    for (int i = 0; i < numseq; i++) {
        names.emplace_back(msa_data->sequences[i]->name);
        aligned_sequences.emplace_back(msa_data->sequences[i]->seq, alnlen);
    }

    kalign_free_msa(msa_data);

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
    const std::string& seq_type = "auto",
    float gap_open = -1.0f,
    float gap_extend = -1.0f,
    float terminal_gap_extend = -1.0f,
    int fill_in = -1,
    float fill_anchor_frac = -1.0f,
    int fill_support = 0
) {
    struct msa* msa_data = nullptr;
    kalign_clear_error();
    int result = kalign_read_input(const_cast<char*>(input_file.c_str()), &msa_data, 1);
    if (result != 0 || !msa_data) {
        throw std::runtime_error(with_reason("Failed to read input file: " + input_file, result == 0 ? "no sequences found" : ""));
    }

    apply_seq_type_or_free(msa_data, seq_type, input_file);

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
    if (!kalign_mode_size_ok(mode.c_str(), msa_data->biotype, msa_data->numseq)) {
        int n = msa_data->numseq;
        kalign_free_msa(msa_data);
        throw std::invalid_argument("Mode '" + mode + "' is not feasible above " +
                                    std::to_string(KALIGN_ACCURATE_MAX_SEQ) + " sequences (input: " +
                                    std::to_string(n) + "). Use mode='fast'.");
    }

    /* Apply user gap penalty overrides to all runs */
    for (int k = 0; k < n_runs; k++) {
        if (gap_open >= 0.0f) runs[k].gpo = gap_open;
        if (gap_extend >= 0.0f) runs[k].gpe = gap_extend;
        if (terminal_gap_extend >= 0.0f) runs[k].tgpe = terminal_gap_extend;
    }
    /* Fill-in overrides (negative = keep preset) */
    if (fill_in >= 0) ens.fill_in = fill_in;
    if (fill_anchor_frac >= 0.0f) ens.fill_anchor_frac = fill_anchor_frac;
    if (fill_support >= 1) ens.fill_support = fill_support;

    kalign_clear_error();
    result = kalign_align_full(msa_data, runs, n_runs,
                                n_runs > 1 ? &ens : nullptr, n_threads);
    if (result != 0) {
        kalign_free_msa(msa_data);
        throw std::runtime_error(with_reason("Alignment failed with error code: " + std::to_string(result)));
    }

    kalign_clear_error();
    result = kalign_write_msa(msa_data, const_cast<char*>(output_file.c_str()),
                               const_cast<char*>(format.c_str()));
    kalign_free_msa(msa_data);

    if (result != 0) {
        throw std::runtime_error(with_reason("Failed to write output file: " + output_file));
    }
}

PYBIND11_MODULE(_core, m) {
    m.doc() =
        "Python bindings for Kalign multiple sequence alignment.\n\n"
        "Supported API: align, align_from_file and align_file_to_file (each "
        "takes mode 'fast' or 'accurate'), add_to_alignment_file, and the "
        "AUTO, PROTEIN, DNA and RNA sequence type constants. Most users "
        "should call these through the kalign package. The helpers "
        "generate_test_sequences, compare, compare_detailed and "
        "compare_detailed_with_mask are utilities for testing and "
        "benchmarking.\n\n"
        "Development interface: ensemble_custom_file_to_file (per-run arrays "
        "for gap penalties, noise, refinement, realignment, sequence weights, "
        "consistency, VSM and guide-tree options), the REFINE_*, FILL_* and "
        "MATRIX_* constants and the sequence type constants other than AUTO, "
        "PROTEIN, DNA and RNA (such as PROTEIN_DIVERGENT and DNA_INTERNAL) "
        "were used to fit the Kalign 4 settings. They are "
        "not part of the supported API and may change without notice. "
        "Settings built this way have not been validated.";
    
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

    // Ensemble with per-run parameters (development interface, see above)
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
          py::arg("vsm_mode") = 0,
          py::arg("vsm_kappa") = 1.0f,
          py::arg("vsm_dzero") = 1.0f,
          py::arg("run_vsm_kappa") = std::vector<float>{},
          py::arg("run_vsm_dzero") = std::vector<float>{},
          py::arg("fill_in") = 0,
          py::arg("fill_anchor_frac") = 0.0f,
          py::arg("fill_support") = 1,
          py::arg("guide_dist") = static_cast<int>(KALIGN_GUIDE_DIST_BPM),
          py::arg("run_guide_dist") = std::vector<int>{},
          py::arg("run_tree_method") = std::vector<int>{},
          py::arg("run_tree_dist") = std::vector<int>{},
          py::arg("run_tree_alpha") = std::vector<int>{},
          py::arg("run_tree_maxlen") = std::vector<int>{},
          py::arg("run_tree_leaf") = std::vector<int>{},
          py::arg("run_tree_kmeans_dist") = std::vector<int>{},
          py::arg("run_realign_tree") = std::vector<int>{},
          py::arg("run_gap_model") = std::vector<int>{},
          py::arg("run_vsm_mode") = std::vector<int>{},
          py::arg("run_tree_sl_noise") = std::vector<float>{},
          R"pbdoc(
          Ensemble alignment with per-run parameters.

          Development interface. This function was used to fit the Kalign 4
          settings. It is not part of the supported API and may change without
          notice. Its refine, realign, seq_weights, consistency, VSM and guide-tree
          options can build configurations that were never fitted or validated.
          To align sequences, use align, align_from_file or align_file_to_file with
          mode="fast" or mode="accurate".

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
          vsm_mode : int, optional
              0 = legacy VSM (offset from vsm_amax and mean distance to anchors,
              the legacy behaviour); 1 = per-merge VSM: offset
              min(kappa,1) * max(0, 2*gpe + mu_bg) * max(0, 1 - d/dzero), with d the
              Kimura-corrected distance between the two groups being merged.
              vsm_amax is ignored when vsm_mode = 1.
          vsm_kappa, vsm_dzero : float, optional
              Per-merge VSM parameters (fraction of the offset bound, 0-1; distance in
              substitutions/site at which the offset reaches 0).
          run_vsm_mode : list of int, optional
              Per-run vsm_mode. Empty = use vsm_mode for all runs.
          run_vsm_kappa, run_vsm_dzero : list of float, optional
              Per-run overrides. Empty = use the shared value for all runs.
          guide_dist : int, optional
              Guide-tree distance. 0 = alignment-free BPM edit distance (default
              behaviour); 1 = Kimura (protein) / Jukes-Cantor (DNA/RNA) corrected
              distance from independent pairwise global alignments with the run's
              matrix and gap penalties. Clustering is unchanged.
          run_guide_dist : list of int, optional
              Per-run guide_dist. Empty = use guide_dist for all runs.
          run_tree_method, run_tree_dist, run_tree_alpha, run_tree_maxlen, run_tree_leaf,
          run_tree_kmeans_dist, run_realign_tree : list of int, optional
              Per-run guide-tree options, see struct tree_opts in lib/src/msa_struct.h.
              Empty or -1 = the default guide tree.
          run_gap_model : list of int, optional
              Per-run gap model: 0 = Kalign gap rule, 1 = member gap model.
              Empty = 0 for all runs.
          run_tree_sl_noise : list of float, optional
              Per-run sigma of the noise on the hybrid tree's single-linkage subtrees
              (tree_method 3), keyed by the run's tree seed. Empty = 0 (no noise).

          Returns
          -------
          dict or None
              For an ensemble with confidence: column_confidence, residue_confidence and
              names (row order of residue_confidence); None otherwise.
          )pbdoc");

    // In-memory alignment using a named mode preset
    m.def("align", &align_mode,
          py::arg("sequences"),
          py::arg("mode"),
          py::arg("seq_type") = "auto",
          py::arg("gap_open") = -1.0f,
          py::arg("gap_extend") = -1.0f,
          py::arg("terminal_gap_extend") = -1.0f,
          py::arg("n_threads") = 1,
          py::arg("fill_in") = -1,
          py::arg("fill_anchor_frac") = -1.0f,
          py::arg("fill_support") = 0,
          "Align sequences using a named mode preset (fast/default/recall/accurate).");
    m.def("align_mode", &align_mode,
          py::arg("sequences"),
          py::arg("mode"),
          py::arg("seq_type") = "auto",
          py::arg("gap_open") = -1.0f,
          py::arg("gap_extend") = -1.0f,
          py::arg("terminal_gap_extend") = -1.0f,
          py::arg("n_threads") = 1,
          py::arg("fill_in") = -1,
          py::arg("fill_anchor_frac") = -1.0f,
          py::arg("fill_support") = 0,
          "Alias for align(). Align sequences using a named mode preset.");

    // File alignment returning (names, sequences) using a named mode preset
    m.def("align_from_file", &align_from_file_mode,
          py::arg("input_file"),
          py::arg("mode"),
          py::arg("seq_type") = "auto",
          py::arg("gap_open") = -1.0f,
          py::arg("gap_extend") = -1.0f,
          py::arg("terminal_gap_extend") = -1.0f,
          py::arg("n_threads") = 1,
          py::arg("fill_in") = -1,
          py::arg("fill_anchor_frac") = -1.0f,
          py::arg("fill_support") = 0,
          "Align from file using a named mode preset. Returns (names, sequences) tuple.");
    m.def("align_from_file_mode", &align_from_file_mode,
          py::arg("input_file"),
          py::arg("mode"),
          py::arg("seq_type") = "auto",
          py::arg("gap_open") = -1.0f,
          py::arg("gap_extend") = -1.0f,
          py::arg("terminal_gap_extend") = -1.0f,
          py::arg("n_threads") = 1,
          py::arg("fill_in") = -1,
          py::arg("fill_anchor_frac") = -1.0f,
          py::arg("fill_support") = 0,
          "Alias for align_from_file(). Align from file using a named mode preset.");

    // File-to-file alignment using a named mode preset
    m.def("align_file_to_file", &align_file_to_file_mode,
          py::arg("input_file"),
          py::arg("output_file"),
          py::arg("mode"),
          py::arg("format") = "fasta",
          py::arg("n_threads") = 1,
          py::arg("seq_type") = "auto",
          py::arg("gap_open") = -1.0f,
          py::arg("gap_extend") = -1.0f,
          py::arg("terminal_gap_extend") = -1.0f,
          py::arg("fill_in") = -1,
          py::arg("fill_anchor_frac") = -1.0f,
          py::arg("fill_support") = 0,
          "Align file to file using a named mode preset (fast/default/recall/accurate).");
    m.def("align_file_to_file_mode", &align_file_to_file_mode,
          py::arg("input_file"),
          py::arg("output_file"),
          py::arg("mode"),
          py::arg("format") = "fasta",
          py::arg("n_threads") = 1,
          py::arg("seq_type") = "auto",
          py::arg("gap_open") = -1.0f,
          py::arg("gap_extend") = -1.0f,
          py::arg("terminal_gap_extend") = -1.0f,
          py::arg("fill_in") = -1,
          py::arg("fill_anchor_frac") = -1.0f,
          py::arg("fill_support") = 0,
          "Alias for align_file_to_file(). Align file to file using a named mode preset.");

    m.def("add_to_alignment_file", [](
            const std::string& existing_file,
            const std::string& new_seqs_file,
            const std::string& output_file,
            const std::string& format,
            int n_threads,
            const std::string& seq_type
    ) {
        struct msa* existing = nullptr;
        struct msa* new_seqs = nullptr;

        kalign_clear_error();
        int result = kalign_read_input(
            const_cast<char*>(existing_file.c_str()), &existing, 1);
        if (result != 0 || !existing) {
            throw std::runtime_error(with_reason("Failed to read existing alignment: " + existing_file, result == 0 ? "no sequences found" : ""));
        }

        kalign_clear_error();
        result = kalign_read_sequences(
            const_cast<char*>(new_seqs_file.c_str()), &new_seqs, 1);
        if (result != 0 || !new_seqs) {
            kalign_free_msa(existing);
            throw std::runtime_error(with_reason("Failed to read new sequences: " + new_seqs_file, result == 0 ? "no sequences found" : ""));
        }

        try {
            apply_seq_type(existing, seq_type, existing_file);
            apply_seq_type(new_seqs, seq_type, new_seqs_file);
        } catch (...) {
            kalign_free_msa(existing);
            kalign_free_msa(new_seqs);
            throw;
        }

        kalign_clear_error();
        result = kalign_add_sequences(existing, new_seqs, n_threads);
        if (result != 0) {
            kalign_free_msa(existing);
            kalign_free_msa(new_seqs);
            throw std::runtime_error(with_reason("Failed to add sequences to alignment"));
        }

        kalign_free_msa(new_seqs);

        kalign_clear_error();
        result = kalign_write_msa(existing,
            const_cast<char*>(output_file.c_str()),
            const_cast<char*>(format.c_str()));
        kalign_free_msa(existing);
        if (result != 0) {
            throw std::runtime_error(with_reason("Failed to write output alignment"));
        }
    },
    py::arg("existing_file"),
    py::arg("new_seqs_file"),
    py::arg("output_file"),
    py::arg("format") = "fasta",
    py::arg("n_threads") = 1,
    py::arg("seq_type") = "auto",
    "Add new sequences to an existing alignment without re-aligning existing sequences.");

    // The matrix, REFINE_* and FILL_* constants, and the sequence type constants
    // other than AUTO, PROTEIN, DNA and RNA, belong to the development interface
    // used to fit the Kalign 4 settings. They are not part of the supported API
    // and may change. AUTO, PROTEIN, DNA and RNA are supported values for
    // seq_type. Use align with mode and seq_type.

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

    // Sequence type constants. AUTO, PROTEIN, DNA and RNA are supported values
    // for seq_type, the same as the strings "auto"/"protein"/"dna"/"rna".
    // DNA_INTERNAL and the PROTEIN_* variants are development interface. Note
    // that ensemble_custom_file_to_file defaults to PROTEIN, not AUTO.
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

    // Constants for refinement modes (development interface, used only by
    // ensemble_custom_file_to_file)
    m.attr("REFINE_NONE") = KALIGN_REFINE_NONE;
    m.attr("REFINE_ALL") = KALIGN_REFINE_ALL;
    m.attr("REFINE_CONFIDENT") = KALIGN_REFINE_CONFIDENT;
    m.attr("REFINE_INLINE") = KALIGN_REFINE_INLINE;

    // Consensus fill-in modes (accurate ensemble only; development interface)
    m.attr("FILL_OFF") = KALIGN_FILL_OFF;
    m.attr("FILL_REALIGN") = KALIGN_FILL_REALIGN;
    m.attr("FILL_PROJECT") = KALIGN_FILL_PROJECT;
    m.attr("FILL_ATTACH") = KALIGN_FILL_ATTACH;

}