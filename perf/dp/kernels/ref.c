/* ref -- the current Kalign DP code (lib/src/aln_seqseq.c,
 * aln_seqprofile.c, aln_profileprofile.c), linked from the regular,
 * uninstrumented kalign_static library.  This is the baseline and the oracle.
 */
#include "kdp_kernel.h"
#include "aln_seqseq.h"
#include "aln_seqprofile.h"
#include "aln_profileprofile.h"

static const struct kdp_kernel ref_kernel = {
        .name = "ref",
        .description = "Kalign 3.7.0 reference DP (lib/src/aln_*.c)",
        .fwd = { aln_seqseq_foward, aln_seqprofile_foward, aln_profileprofile_foward },
        .bwd = { aln_seqseq_backward, aln_seqprofile_backward, aln_profileprofile_backward },
        .meet = { aln_seqseq_meetup, aln_seqprofile_meetup, aln_profileprofile_meetup },
};
KDP_REGISTER_KERNEL(ref_kernel)
