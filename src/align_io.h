#ifndef ALIGN_IO_H
#define ALIGN_IO_H

#include <unistd.h>

#define SEEK_START 0
#define SEEK_END 2



extern struct alignment* detect_and_read_sequences(struct parameters* param);
extern void free_aln(struct alignment* aln);

#endif
