#include <math.h>
#include <string.h>

#include "tldevel.h"

#define PROBMSA_DIRICHLET_IMPORT
#include "probmsa_dirichlet.h"

/*
 * Sjölander 20-component Dirichlet mixture prior (recode3.20comp).
 * Alphabet: ACDEFGHIKLMNPQRSTVWY (indices 0-19).
 * Source: MEME suite glam2_recode3_20comp.c (UCSC SAM toolkit).
 */

static const double MIX_WEIGHTS[DIRI_NUM_COMPONENTS] = {
        0.176513, 0.207622, 0.0669246, 0.0868259, 0.0593123,
        0.0358616, 0.03427, 0.0428319, 0.047875, 0.0466614,
        0.0283695, 0.0301127, 0.0233828, 0.034662, 0.0270202,
        0.0226822, 0.00898452, 0.00716226, 0.00710292, 0.00582299
};

/* Alpha values: [component][amino_acid] in ACDEFGHIKLMNPQRSTVWY order */
static const double MIX_ALPHA[DIRI_NUM_COMPONENTS][DIRI_ALPHA_SIZE] = {
        /* Component 0: general polar/charged (w=0.177) */
        { 0.668781, 0.072626, 1.08101, 0.930936, 0.1235,
          0.982387, 0.366701, 0.0586812, 0.849019, 0.185073,
          0.0648062, 0.979541, 0.519704, 0.555784, 0.636552,
          1.04987, 0.60931, 0.129148, 0.0676805, 0.201772 },
        /* Component 1: background/general (w=0.208, largest) */
        { 0.880009, 0.197393, 0.303744, 0.497357, 0.657492,
          0.391988, 0.348159, 0.936094, 0.527206, 1.42844,
          0.449302, 0.368825, 0.384576, 0.439775, 0.581516,
          0.622624, 0.747064, 1.08007, 0.235012, 0.606928 },
        /* Component 2: hydrophobic core L,I,F,M,V (w=0.067) */
        { 0.153384, 0.0520756, 0.0073824, 0.0158439, 0.428964,
          0.025533, 0.0185789, 0.845361, 0.0282996, 2.42256,
          0.424296, 0.0190716, 0.0313429, 0.0274578, 0.0252186,
          0.028514, 0.0519217, 0.522946, 0.0279653, 0.0664755 },
        /* Component 3: charged/polar exposed E,K,R,Q (w=0.087) */
        { 0.794007, 0.0130659, 0.624236, 1.85769, 0.0290214,
          0.115707, 0.123504, 0.22099, 1.52605, 0.341371,
          0.111114, 0.308302, 0.263545, 0.953727, 0.933444,
          0.554741, 0.604551, 0.396451, 0.00823516, 0.0420054 },
        /* Component 4: aliphatic V,I,A,L (w=0.059) */
        { 0.740015, 0.187165, 0.0213261, 0.0456854, 0.118944,
          0.0633687, 0.0170331, 1.06684, 0.0380614, 0.733524,
          0.138456, 0.0300644, 0.0718692, 0.0240143, 0.0301022,
          0.0862989, 0.367283, 1.70735, 0.0113856, 0.045079 },
        /* Component 5: basic R,K (w=0.036) */
        { 0.15978, 0.0261585, 0.0505181, 0.125524, 0.0350331,
          0.102549, 0.157461, 0.0795041, 1.26261, 0.189383,
          0.0550608, 0.171028, 0.0844169, 0.290476, 1.44604,
          0.129158, 0.138972, 0.0851144, 0.0159134, 0.0637679 },
        /* Component 6: acidic D,E (w=0.034) */
        { 0.308434, 0.0137217, 1.69731, 1.92422, 0.0361113,
          0.162357, 0.07232, 0.0487895, 0.236135, 0.0809074,
          0.0286236, 0.213663, 0.181631, 0.320245, 0.104878,
          0.218398, 0.141668, 0.0747719, 0.0141705, 0.0453433 },
        /* Component 7: near-zero (highly conserved positions) (w=0.043) */
        { 0.00260287, 9.99856e-06, 0.00631292, 0.00445502, 0.00274753,
          1.03886e-05, 1.02839e-05, 0.000913052, 0.0029241, 0.00353485,
          0.00105128, 0.00338172, 1.04172e-05, 0.00173574, 0.00459583,
          0.00274255, 0.00247625, 0.00175366, 1.02411e-05, 0.00288489 },
        /* Component 8: A,G,S rich (w=0.048) */
        { 1.61043, 0.15522, 0.0378292, 0.0498243, 0.0406484,
          0.529136, 0.0217524, 0.040597, 0.0413396, 0.100193,
          0.0509779, 0.0357917, 0.0931204, 0.0367156, 0.0330646,
          0.529587, 0.196607, 0.230878, 0.00909518, 0.0329275 },
        /* Component 9: glycine-dominated G=3.11 (w=0.047) */
        { 0.15525, 0.0136827, 0.0857138, 0.0508316, 0.0151451,
          3.10555, 0.027169, 0.0140491, 0.0654038, 0.0257501,
          0.00901049, 0.127437, 0.0423873, 0.0345064, 0.0477247,
          0.12452, 0.0341196, 0.0230637, 0.00930115, 0.0187464 },
        /* Component 10: S,T (w=0.028) */
        { 0.225739, 0.0684326, 0.101072, 0.0813791, 0.0298832,
          0.0915218, 0.0336807, 0.0833114, 0.0931673, 0.0731542,
          0.0419314, 0.230216, 0.087446, 0.0694702, 0.0751969,
          1.13857, 1.63158, 0.179083, 0.00912576, 0.0311963 },
        /* Component 11: I,V,L hydrophobic (w=0.030) */
        { 1.3431e-06, 0.0166656, 0.00743068, 1.34592e-06, 0.169076,
          0.00407061, 0.00714122, 1.98221, 0.017522, 0.816669,
          0.114773, 0.00678027, 0.0106392, 0.0100244, 0.0158968,
          0.00879658, 0.0399043, 1.81043, 0.0150671, 0.0517801 },
        /* Component 12: D,N (w=0.023) */
        { 0.063525, 0.0288391, 1.09265, 0.0959581, 0.00965196,
          0.216914, 0.0730986, 0.0207325, 0.08719, 0.0315107,
          0.00960293, 0.752755, 0.059914, 0.0445321, 0.0312317,
          0.287327, 0.116896, 0.0249446, 0.00701663, 0.0305859 },
        /* Component 13: proline-dominated P=3.70 (w=0.035) */
        { 0.294281, 0.019271, 0.12293, 0.162747, 0.0373667,
          0.145029, 0.0412349, 0.0815261, 0.157594, 0.151631,
          0.021412, 0.0601581, 3.6966, 0.0809085, 0.101856,
          0.23533, 0.135424, 0.140532, 0.00900473, 0.0321389 },
        /* Component 14: aromatic Y,F,H,W (w=0.027) */
        { 0.0844832, 0.0584945, 0.0411628, 0.045719, 0.847822,
          0.0590839, 0.250253, 0.0675757, 0.0562614, 0.168617,
          0.0439737, 0.0794234, 0.028301, 0.0305672, 0.0598024,
          0.0798202, 0.0585385, 0.0858243, 0.227395, 1.30336 },
        /* Component 15: F,L,Y,W aromatic hydrophobic (w=0.023) */
        { 0.0634034, 0.0246167, 1.3443e-06, 0.00389272, 1.12953,
          0.00796028, 1.35032e-06, 0.233395, 1.34466e-06, 0.541933,
          0.101309, 1.36412e-06, 0.027467, 0.00704479, 0.00802297,
          0.0248977, 0.0276933, 0.185467, 0.183309, 0.516892 },
        /* Component 16: Q,E,H,R (w=0.009) */
        { 0.123696, 0.0454619, 0.0386434, 0.351847, 0.0560181,
          0.0439442, 0.223229, 0.01302, 0.148699, 0.19001,
          0.120964, 0.098734, 5.90055e-06, 0.554971, 0.219233,
          0.0453885, 0.0564686, 0.0614792, 0.0410248, 0.0800036 },
        /* Component 17: cysteine-dominated C=3.19 (w=0.007) */
        { 0.0212037, 3.18769, 0.00745627, 0.00382411, 0.00691924,
          0.0126233, 1.34375e-06, 0.00724293, 0.00522979, 0.00785563,
          0.00489521, 0.0105326, 0.0136265, 0.00505819, 0.00677712,
          0.0251744, 0.0235516, 0.0371462, 0.00187667, 0.0038893 },
        /* Component 18: tryptophan-dominated W=1.84 (w=0.007) */
        { 0.0229376, 0.00427768, 0.00959934, 0.013608, 0.182277,
          0.0227654, 0.0157344, 0.0226783, 0.011561, 0.0803491,
          0.0154283, 0.00899225, 0.00980608, 0.00600945, 0.0342359,
          0.0216842, 0.0189306, 0.0223176, 1.83914, 0.154565 },
        /* Component 19: histidine-dominated H=1.03 (w=0.006) */
        { 2.16602e-06, 2.16245e-06, 0.0198496, 2.17942e-06, 0.0246741,
          2.47051e-06, 1.02563, 0.0131152, 2.16539e-06, 0.00637704,
          2.1414e-06, 0.0839371, 0.0168135, 0.0438887, 0.0252951,
          0.0235533, 0.0130626, 0.00797507, 2.16433e-06, 0.0545531 }
};

int probmsa_diri_init(struct probmsa_diri_mix** out, double prior_scale)
{
        struct probmsa_diri_mix* m = NULL;
        int k, i;

        ASSERT(prior_scale > 0.0, "prior_scale must be > 0");

        MMALLOC(m, sizeof(struct probmsa_diri_mix));

        m->n_comp = DIRI_NUM_COMPONENTS;
        m->alpha_size = DIRI_ALPHA_SIZE;

        for(k = 0; k < DIRI_NUM_COMPONENTS; k++){
                m->q[k] = MIX_WEIGHTS[k];
                m->alpha_sum[k] = 0.0;
                for(i = 0; i < DIRI_ALPHA_SIZE; i++){
                        m->alpha[k][i] = MIX_ALPHA[k][i] * prior_scale;
                        m->log_gamma_alpha[k][i] = lgamma(m->alpha[k][i]);
                        m->alpha_sum[k] += m->alpha[k][i];
                }
                m->log_gamma_alpha_sum[k] = lgamma(m->alpha_sum[k]);
        }

        *out = m;
        return OK;
ERROR:
        return FAIL;
}

void probmsa_diri_free(struct probmsa_diri_mix* m)
{
        if(m){
                MFREE(m);
        }
}

/*
 * Compute log Prob(n | alpha_k) using equation 51 from Sjölander 1996:
 *
 * log P(n|alpha) = lgamma(|n|+1) + lgamma(|alpha|) - lgamma(|n|+|alpha|)
 *                  + sum_i [lgamma(n_i + alpha_i) - lgamma(n_i + 1) - lgamma(alpha_i)]
 *
 * We drop the lgamma(|n|+1) and lgamma(n_i+1) terms since they cancel
 * in the posterior ratio (equation 16). This gives:
 *
 * log P(n|alpha) ∝ lgamma(|alpha|) - lgamma(|n|+|alpha|)
 *                  + sum_i [lgamma(n_i + alpha_i) - lgamma(alpha_i)]
 */
static double log_prob_counts_given_component(const struct probmsa_diri_mix* m,
                                               int k,
                                               const double* counts,
                                               double n_total)
{
        double logp;
        int i;

        logp = m->log_gamma_alpha_sum[k] - lgamma(n_total + m->alpha_sum[k]);

        for(i = 0; i < DIRI_ALPHA_SIZE; i++){
                if(counts[i] > 0.0){
                        logp += lgamma(counts[i] + m->alpha[k][i]) - m->log_gamma_alpha[k][i];
                }
                /* When counts[i] == 0, lgamma(0 + alpha) - lgamma(alpha) = 0 */
        }

        return logp;
}

int probmsa_diri_joint_matrix(const struct probmsa_diri_mix* m,
                               double joint[][DIRI_ALPHA_SIZE],
                               double* background)
{
        int k, s, t;

        /* Zero output */
        for(s = 0; s < DIRI_ALPHA_SIZE; s++){
                background[s] = 0.0;
                for(t = 0; t < DIRI_ALPHA_SIZE; t++){
                        joint[s][t] = 0.0;
                }
        }

        /* joint[s][t] = Σ_k w_k * f_k_s * f_k_t
         * bg[s] = Σ_k w_k * f_k_s
         * where f_k_s = α_k_s / |α_k| */
        for(k = 0; k < m->n_comp; k++){
                double w = m->q[k];
                double asum = m->alpha_sum[k];
                for(s = 0; s < DIRI_ALPHA_SIZE; s++){
                        double fs = m->alpha[k][s] / asum;
                        background[s] += w * fs;
                        for(t = 0; t < DIRI_ALPHA_SIZE; t++){
                                double ft = m->alpha[k][t] / asum;
                                joint[s][t] += w * fs * ft;
                        }
                }
        }

        return OK;
}

int probmsa_diri_posterior_mean(const struct probmsa_diri_mix* m,
                                const double* counts, double n_total,
                                double* prob_out,
                                double* concentration_out)
{
        double log_comp_prob[DIRI_NUM_COMPONENTS];
        double comp_weight[DIRI_NUM_COMPONENTS];
        double max_log, sum_weight;
        int k, i;

        /* Compute log P(n|alpha_k) + log(q_k) for each component */
        max_log = -INFINITY;
        for(k = 0; k < DIRI_NUM_COMPONENTS; k++){
                log_comp_prob[k] = log(m->q[k])
                        + log_prob_counts_given_component(m, k, counts, n_total);
                if(log_comp_prob[k] > max_log){
                        max_log = log_comp_prob[k];
                }
        }

        /* Convert to normalized weights using log-sum-exp for numerical stability */
        sum_weight = 0.0;
        for(k = 0; k < DIRI_NUM_COMPONENTS; k++){
                comp_weight[k] = exp(log_comp_prob[k] - max_log);
                sum_weight += comp_weight[k];
        }
        for(k = 0; k < DIRI_NUM_COMPONENTS; k++){
                comp_weight[k] /= sum_weight;
        }

        /* Posterior mean: p_i = sum_k w_k * (n_i + alpha_ki) / (|n| + |alpha_k|) */
        for(i = 0; i < DIRI_ALPHA_SIZE; i++){
                prob_out[i] = 0.0;
        }

        for(k = 0; k < DIRI_NUM_COMPONENTS; k++){
                if(comp_weight[k] < 1e-15) continue;  /* skip negligible components */
                double denom = n_total + m->alpha_sum[k];
                for(i = 0; i < DIRI_ALPHA_SIZE; i++){
                        prob_out[i] += comp_weight[k] * (counts[i] + m->alpha[k][i]) / denom;
                }
        }

        /* Effective concentration: weighted sum of (n_total + |alpha_k|) */
        if(concentration_out){
                double conc = 0.0;
                for(k = 0; k < DIRI_NUM_COMPONENTS; k++){
                        conc += comp_weight[k] * (n_total + m->alpha_sum[k]);
                }
                *concentration_out = conc;
        }

        return OK;
}
