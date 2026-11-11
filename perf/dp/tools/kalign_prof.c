/* kalign_prof -- profiling driver: read -> mode preset -> align -> write.
 *
 * Built twice:
 *   kalign_prof   linked to kalign_perf_static (KALIGN_PERF_TIMERS): writes a
 *                 JSON breakdown of where time goes (-j file).
 *   kalign_plain  linked to the regular kalign_static: wall times only
 *                 (overhead / output-parity reference).
 *
 * Usage: kalign_prof -i in.fa -o out.fa [-m fast|default|accurate] [-n threads]
 *                    [-f fasta] [-j timers.json]
 * The alignment is the same as `kalign -i in.fa -o out.fa --mode <m> -n <t>`
 * (same preset lookup, no parameter overrides).
 */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include <unistd.h>

#include "kalign/kalign.h"
#include "kalign/kalign_config.h"

#ifdef KALIGN_PERF_TIMERS
#include "kperf.h"
#endif

static double now_s(void)
{
        struct timespec ts;
        clock_gettime(CLOCK_MONOTONIC, &ts);
        return (double)ts.tv_sec + 1e-9 * (double)ts.tv_nsec;
}

static void usage(void)
{
        fprintf(stderr, "usage: kalign_prof -i in.fa -o out.fa [-m mode] [-n threads] [-f format] [-j timers.json]\n");
}

static int run(char* in, char* out, const char* mode, int nthreads, char* fmt, const char* json)
{
        struct msa* msa = NULL;
        struct kalign_run_config runs[KALIGN_MAX_PRESET_RUNS];
        struct kalign_ensemble_config ens;
        int n_runs = 0;
        double t0, t1, t2, t3;
        int rc = 1;

        t0 = now_s();
        {
#ifdef KALIGN_PERF_TIMERS
                KPERF_REGION(KP_DRIVER);
                {
                        KPERF_REGION(KP_IO_READ);
#endif
                        if(kalign_read_input(in, &msa, 1) != 0){
                                fprintf(stderr, "read failed: %s\n", in);
                                return 1;
                        }
#ifdef KALIGN_PERF_TIMERS
                }
#endif
                t1 = now_s();
                if(kalign_get_mode_preset(mode, kalign_msa_get_biotype(msa), runs, &n_runs, &ens) != 0){
                        fprintf(stderr, "unknown mode %s\n", mode);
                        kalign_free_msa(msa);
                        return 1;
                }
                if(kalign_align_full(msa, runs, n_runs, &ens, nthreads) != 0){
                        fprintf(stderr, "alignment failed\n");
                        kalign_free_msa(msa);
                        return 1;
                }
                t2 = now_s();
#ifdef KALIGN_PERF_TIMERS
                {
                        KPERF_REGION(KP_IO_WRITE);
#endif
                        if(kalign_write_msa(msa, out, fmt) != 0){
                                fprintf(stderr, "write failed\n");
                                kalign_free_msa(msa);
                                return 1;
                        }
#ifdef KALIGN_PERF_TIMERS
                }
#endif
                t3 = now_s();
        }
        rc = 0;
        fprintf(stderr, "kalign_prof: mode=%s n_runs=%d read=%.3fs align=%.3fs write=%.3fs total=%.3fs\n",
                mode, n_runs, t1 - t0, t2 - t1, t3 - t2, t3 - t0);
        if(json){
                FILE* f = fopen(json, "w");
                if(!f){
                        fprintf(stderr, "cannot write %s\n", json);
                        rc = 1;
                }else{
                        fprintf(f, "{\"mode\": \"%s\", \"n_runs\": %d, \"nthreads\": %d, \"input\": \"%s\",\n"
                                " \"wall\": {\"read\": %.6f, \"align\": %.6f, \"write\": %.6f, \"total\": %.6f},\n"
                                " \"instrumented\": %d,\n \"timers\": ",
                                mode, n_runs, nthreads, in, t1 - t0, t2 - t1, t3 - t2, t3 - t0,
#ifdef KALIGN_PERF_TIMERS
                                1
#else
                                0
#endif
                                );
#ifdef KALIGN_PERF_TIMERS
                        kperf_write_json(f);
#else
                        fprintf(f, "null\n");
#endif
                        fprintf(f, "}\n");
                        fclose(f);
                }
        }
        kalign_free_msa(msa);
        return rc;
}

int main(int argc, char** argv)
{
        char* in = NULL;
        char* out = NULL;
        const char* mode = "default";
        char* fmt = NULL; /* same default as the kalign CLI */
        const char* json = NULL;
        int nthreads = 1;
        int c;
        while((c = getopt(argc, argv, "i:o:m:n:f:j:h")) != -1){
                switch(c){
                case 'i': in = optarg; break;
                case 'o': out = optarg; break;
                case 'm': mode = optarg; break;
                case 'n': nthreads = atoi(optarg); break;
                case 'f': fmt = optarg; break;
                case 'j': json = optarg; break;
                default: usage(); return 1;
                }
        }
        if(!in || !out){
                usage();
                return 1;
        }
        return run(in, out, mode, nthreads, fmt, json);
}
