#ifndef PICK_ANCHOR_H
#define PICK_ANCHOR_H


#ifdef PICK_ANCHOR_IMPORT
#define EXTERN
#else
#ifdef __cplusplus
#define EXTERN extern "C"
#else
#define EXTERN extern
#endif
#endif

struct msa;
EXTERN int* pick_anchor(struct msa* msa, int* n);
EXTERN int* pick_anchor_n(struct msa* msa, int requested, int* n);
/* pick_anchor plus the N x K matrix of BPM distances it computed (-1 where
   not computed); pass it to d_estimation_cached, free with
   pick_anchor_dist_free. */
EXTERN int* pick_anchor_dist(struct msa* msa, int* n, float*** dist_out);
EXTERN void pick_anchor_dist_free(float** dist, int numseq);


#undef PICK_ANCHOR_IMPORT
#undef EXTERN

#endif
