/* #include "tldevel.h" */
#include "stdio.h"
#include "string.h"
/* #include "tlrng.h" */

#include <stddef.h>

/* #include "strnlen_compat.h" */

int main(int argc, char *argv[])
{
        /* struct rng_state* rng = NULL; */
        /* rng = init_rng(0); */
        /* fprintf(stdout,"Kalign (%s)\n", KALIGN_PACKAGE_VERSION); */
        /* fprintf(stdout,"HEllo random: %f\n", tl_random_double(rng)); */

        int len = strnlen("FAFAF", 40);
        fprintf(stdout,"Len: %d\n",len);
        return 0;
}
