#ifndef MSA_IO_H
#define MSA_IO_H

#ifdef MSA_IO_IMPORT
#define EXTERN
#else
#define EXTERN extern
#endif

#define FORMAT_DETECT_FAIL -1
#define FORMAT_FA 1
#define FORMAT_MSF 2
#define FORMAT_CLU 3

struct msa;


EXTERN int kalign_read_input(char* infile, struct msa** msa,int quiet);
EXTERN int kalign_write_msa(struct msa* msa, char* outfile, char* format);



#undef MSA_IO_IMPORT
#undef EXTERN


#endif
