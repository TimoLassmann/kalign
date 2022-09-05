#ifndef LIBKALIGN_H
#define LIBKALIGN_H

#ifdef __cplusplus
#define kalign_extern extern "C"
#else
#define kalign_extern extern
#endif



/* #include "mod_tldevel.h" */
/* #include "mod_interface.h" */
/* #include "task.h" */
/* #include "msa.h" */
/* #include "mod_msaio.h" */
/* #include "mod_tree.h" */
/* #include "mod_aln.h" */
/* #include "alphabet.h" */



kalign_extern int kalign(char **seq, int *len,int numseq, char ***aligned, int *out_aln_len);

#undef LIBKALIGN_IMPORT
#undef EXTERN


#endif
