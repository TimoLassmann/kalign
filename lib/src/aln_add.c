/* aln_add.c — Add new sequences to an existing alignment.
 *
 * Builds a consensus profile from the existing aligned sequences, then
 * aligns each new sequence against it via seq-to-profile Hirschberg DP.
 * The existing sequences are NOT modified — only gap characters are
 * inserted into the new sequences to fit the existing column structure.
 */

#include "tldevel.h"
#include <ctype.h>
#include <string.h>
#include <float.h>

#include "msa_struct.h"
#include "msa_alloc.h"
#include "msa_op.h"
#include "msa_check.h"
#include "alphabet.h"

#include "aln_param.h"
#include "aln_struct.h"
#include "aln_mem.h"
#include "aln_setup.h"
#include "aln_controller.h"
#include "kalign/kalign.h"

#ifdef USE_THREADPOOL
#include "threadpool/threadpool.h"
#endif

#define ALN_ADD_IMPORT
#include "aln_add.h"

/* Map a character (from finalized alignment) to internal index 0-22.
   Returns -1 for gaps and unknown characters.
   Uses the standard kalign protein alphabet: ARNDCQEGHILKMFPSTWYVBZX */
static int char_to_internal(char c, int biotype)
{
        static const int protein_map[128] = {
                /*   0-15 */ -1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,
                /*  16-31 */ -1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,
                /*  32-47 */ -1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,
                /*  48-63 */ -1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,
                /* @ABCDE */  -1, 0,20, 4, 3, 6,13, 7, 8, 9,-1,11,10,12, 2,-1,
                /* PQRSTU */  14, 5, 1,15,16,22,19,17,22,18,21,-1,-1,-1,-1,-1,
                /* `abcde */  -1, 0,20, 4, 3, 6,13, 7, 8, 9,-1,11,10,12, 2,-1,
                /* pqrstu */  14, 5, 1,15,16,22,19,17,22,18,21,-1,-1,-1,-1,-1,
        };
        static const int dna_map[128] = {
                -1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,
                -1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,
                -1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,
                -1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,
                -1, 0,-1, 1,-1,-1,-1, 2,-1,-1,-1,-1,-1,-1,-1,-1,
                -1,-1,-1,-1, 3, 3,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,
                -1, 0,-1, 1,-1,-1,-1, 2,-1,-1,-1,-1,-1,-1,-1,-1,
                -1,-1,-1,-1, 3, 3,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,
        };
        if(c < 0 || c == '-') return -1;
        if(biotype == ALN_BIOTYPE_DNA){
                return dna_map[(unsigned char)c];
        }
        return protein_map[(unsigned char)c];
}

/* Build a consensus profile from an existing finalized alignment.
   Profile has alnlen columns in the standard 64-float-per-column format.
   Residue frequencies are summed; substitution scores are weighted averages. */
static int build_consensus_profile(struct msa* existing,
                                    struct aln_param* ap,
                                    float** prof_out)
{
        float* prof = NULL;
        float** subm = ap->subm;
        float gpo = ap->gpo;
        float gpe = ap->gpe;
        float tgpe = ap->tgpe;
        int alnlen = existing->alnlen;
        int numseq = existing->numseq;
        int biotype = existing->biotype;
        int i, j, col;

        MMALLOC(prof, sizeof(float) * (alnlen + 2) * 64);

        /* Point to trailing boundary row */
        prof += (64 * (alnlen + 1));

        /* Trailing boundary */
        for(i = 0; i < 64; i++) prof[i] = 0;
        prof[23 + 32] = -gpo;
        prof[24 + 32] = -gpe;
        prof[25 + 32] = -tgpe;
        prof[55] = -gpo;
        prof[56] = -gpe;
        prof[57] = -tgpe;

        /* Fill each column from alnlen-1 down to 0 */
        for(col = alnlen - 1; col >= 0; col--){
                prof -= 64;
                for(j = 0; j < 64; j++) prof[j] = 0;

                /* Count residue frequencies at this column */
                float freq[23];
                for(j = 0; j < 23; j++) freq[j] = 0.0f;
                float total = 0.0f;

                for(i = 0; i < numseq; i++){
                        char c = existing->sequences[i]->seq[col];
                        int idx = char_to_internal(c, biotype);
                        if(idx >= 0 && idx < 23){
                                freq[idx] += 1.0f;
                                total += 1.0f;
                        }
                }

                /* Store frequency counts in positions 0-22 */
                for(j = 0; j < 23; j++){
                        prof[j] = freq[j];
                }

                /* Compute substitution scores as weighted average over
                   observed residues: score[j] = sum(freq[c] * subm[c][j]) / total */
                prof += 32;
                if(total > 0.0f){
                        for(j = 0; j < 23; j++){
                                float score = 0.0f;
                                for(int c = 0; c < 23; c++){
                                        if(freq[c] > 0.0f){
                                                score += freq[c] * subm[c][j];
                                        }
                                }
                                prof[j] = score / total;
                        }
                }
                /* Gap penalties */
                prof[23] = -gpo;
                prof[24] = -gpe;
                prof[25] = -tgpe;
                prof -= 32;

                /* Store base penalties for set_gap_penalties_n */
                prof[55] = -gpo;
                prof[56] = -gpe;
                prof[57] = -tgpe;
        }

        /* Leading boundary */
        prof -= 64;
        for(i = 0; i < 64; i++) prof[i] = 0;
        prof[23 + 32] = -gpo;
        prof[24 + 32] = -gpe;
        prof[25 + 32] = -tgpe;
        prof[55] = -gpo;
        prof[56] = -gpe;
        prof[57] = -tgpe;

        *prof_out = prof;
        return OK;
ERROR:
        return FAIL;
}

/* Align one new sequence against the consensus profile and produce
   a gapped sequence string matching the existing alignment columns.
   Returns a newly allocated string (caller must free). */
static int align_one_to_profile(struct aln_param* ap,
                                 float* cons_profile,
                                 int profile_len,
                                 int n_existing,
                                 const uint8_t* new_seq,
                                 int new_len,
                                 char* new_seq_chars,
                                 char** gapped_out)
{
        struct aln_mem* m = NULL;
        char* gapped = NULL;
        int i, c;
        int pos_seq;

        RUN(alloc_aln_mem(&m, 256));
        m->ap = ap;
        m->mode = ALN_MODE_FULL;
        m->run_parallel = 0;
        m->flip_threshold = 0.0F;
        m->flip_trial = 0;
        m->flip_stride = 1;
        m->flip_counter = 0;
        m->flip_mask = 0;
        m->margin_sum = 0.0F;
        m->margin_count = 0;

        /* Profile is always "seq1" (the longer axis in Hirschberg).
           New sequence is "seq2". */
        m->len_a = profile_len;
        m->len_b = new_len;
        m->enda = profile_len;
        m->endb = new_len;

        m->seq1 = NULL;  /* not a sequence — it's a profile */
        m->seq2 = new_seq;
        m->prof1 = cons_profile;
        m->prof2 = NULL;
        m->sip = n_existing;
        m->consistency = NULL;

        m->f[0].a = 0.0F;
        m->f[0].ga = -FLT_MAX;
        m->f[0].gb = -FLT_MAX;
        m->b[0].a = 0.0F;
        m->b[0].ga = -FLT_MAX;
        m->b[0].gb = -FLT_MAX;

        /* Scale gap penalties in profile by n_existing */
        RUN(set_gap_penalties_n(cons_profile, profile_len, n_existing));

        RUN(init_alnmem(m));
        aln_runner(m);
        RUN(add_gap_info_to_path_n(m));

        /* Build gapped sequence from alignment path.
           path[0] = alignment length
           path[c]: 0=match, &1=gap in profile (insertion in new seq — SKIP in strict mode),
                    &2=gap in new seq (insert '-'), 3=end */
        MMALLOC(gapped, sizeof(char) * (m->path[0] + 2));

        pos_seq = 0;
        i = 0;
        c = 1;
        while(m->path[c] != 3){
                if(m->path[c] == 0){
                        /* Match: new seq residue aligns to profile column */
                        if(pos_seq < new_len){
                                gapped[i] = new_seq_chars[pos_seq];
                        }else{
                                gapped[i] = '-';
                        }
                        pos_seq++;
                        i++;
                }else if(m->path[c] & 1){
                        /* Gap in profile = insertion in new seq.
                           Strict mode: skip this residue (don't add new columns). */
                        pos_seq++;
                        /* Don't increment i — residue is dropped */
                }else if(m->path[c] & 2){
                        /* Gap in new seq: insert gap at this profile column */
                        gapped[i] = '-';
                        i++;
                }
                c++;
        }
        gapped[i] = '\0';

        /* Verify: gapped length should equal profile_len (strict mode) */
        if(i != profile_len){
                /* May differ slightly — pad or truncate */
                while(i < profile_len){
                        gapped[i] = '-';
                        i++;
                }
                gapped[profile_len] = '\0';
        }

        free_aln_mem(m);
        *gapped_out = gapped;
        return OK;
ERROR:
        if(m) free_aln_mem(m);
        if(gapped) MFREE(gapped);
        return FAIL;
}

int kalign_add_sequences(struct msa* existing,
                          struct msa* new_seqs,
                          int n_threads)
{
        struct aln_param* ap = NULL;
        float* cons_profile = NULL;
        int i;
        int numseq_existing;
        int numseq_new;
        int alnlen;

        ASSERT(existing != NULL, "No existing alignment");
        ASSERT(new_seqs != NULL, "No new sequences");
        /* Finalize if not already — ensures seq->seq has gap characters and alnlen is set */
        if(existing->aligned != ALN_STATUS_FINAL){
                if(existing->aligned == ALN_STATUS_ALIGNED){
                        RUN(finalise_alignment(existing));
                }else{
                        ERROR_MSG("Existing MSA must be aligned (status=%d)", existing->aligned);
                }
        }

        numseq_existing = existing->numseq;
        numseq_new = new_seqs->numseq;
        alnlen = existing->alnlen;

        if(numseq_new == 0){
                return OK;  /* nothing to add */
        }

        /* Detect biotype if needed */
        if(existing->biotype == ALN_BIOTYPE_UNDEF){
                RUN(detect_alphabet(existing));
        }
        if(new_seqs->biotype == ALN_BIOTYPE_UNDEF){
                new_seqs->biotype = existing->biotype;
        }

        /* Encode new sequences to internal representation */
        if(existing->biotype == ALN_BIOTYPE_DNA){
                new_seqs->L = ALPHA_defDNA;
                RUN(convert_msa_to_internal(new_seqs, ALPHA_defDNA));
        }else{
                new_seqs->L = ALPHA_ambigiousPROTEIN;
                RUN(convert_msa_to_internal(new_seqs, ALPHA_ambigiousPROTEIN));
        }

        /* Init alignment parameters */
        RUN(aln_param_init(&ap, existing->biotype, n_threads,
                           KALIGN_MATRIX_AUTO, -1.0f, -1.0f, -1.0f));

        /* Build consensus profile from existing alignment */
        RUN(build_consensus_profile(existing, ap, &cons_profile));

        /* Align each new sequence to the consensus profile */
        for(i = 0; i < numseq_new; i++){
                char* gapped = NULL;

                RUN(align_one_to_profile(
                        ap, cons_profile, alnlen, numseq_existing,
                        new_seqs->sequences[i]->s,
                        new_seqs->sequences[i]->len,
                        new_seqs->sequences[i]->seq,
                        &gapped));

                /* Replace the new sequence's seq with the gapped version */
                MFREE(new_seqs->sequences[i]->seq);
                new_seqs->sequences[i]->seq = gapped;
        }

        /* Append new sequences to existing MSA (manual, not merge_msa which
           may change alignment properties). */
        {
                int total = numseq_existing + numseq_new;
                if(total > existing->alloc_numseq){
                        MREALLOC(existing->sequences,
                                 sizeof(struct msa_seq*) * total);
                        existing->alloc_numseq = total;
                }
                for(i = 0; i < numseq_new; i++){
                        /* Transfer ownership of the new sequence from new_seqs to existing */
                        existing->sequences[numseq_existing + i] = new_seqs->sequences[i];
                        new_seqs->sequences[i] = NULL;
                        /* Update len to alnlen (gapped length) */
                        existing->sequences[numseq_existing + i]->len =
                                (int)strlen(existing->sequences[numseq_existing + i]->seq);
                }
                existing->numseq = total;
        }
        existing->alnlen = alnlen;
        existing->aligned = ALN_STATUS_FINAL;

        /* Cleanup */
        /* Free the profile — need to adjust pointer back to allocation start */
        {
                float* prof_base = cons_profile; /* already points to start (leading boundary) */
                MFREE(prof_base);
        }
        aln_param_free(ap);

        return OK;
ERROR:
        if(cons_profile) MFREE(cons_profile);
        if(ap) aln_param_free(ap);
        return FAIL;
}
