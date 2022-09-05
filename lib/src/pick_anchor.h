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


#undef PICK_ANCHOR_IMPORT
#undef EXTERN

#endif
