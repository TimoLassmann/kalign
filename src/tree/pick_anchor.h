#ifndef PICK_ANCHOR_H
#define PICK_ANCHOR_H

/* #include "global.h" */
#include "msa.h"
/* #include "alignment_parameters.h" */
/* kalign_extern */

#ifndef kalign_extern
#ifdef __cplusplus
#define kalign_extern extern "C"
#else
#define kalign_extern extern
#endif
#endif

kalign_extern int* pick_anchor(struct msa* msa, int* n);

#endif
