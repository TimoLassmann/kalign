#include "tldevel.h"
#include "kalign/kalign.h"

#include <string.h>
#ifdef _WIN32
#include <windows.h>
#else
#include <unistd.h>
#endif
#ifdef HAVE_OPENMP
#include <omp.h>
#endif
#define PARAMETERS_IMPORT
#include "parameters.h"

static int get_default_thread_count(void)
{
        int cores = 1;

#ifdef HAVE_OPENMP
        cores = omp_get_num_procs();
#elif defined(_WIN32)
        SYSTEM_INFO sysinfo;
        GetSystemInfo(&sysinfo);
        cores = sysinfo.dwNumberOfProcessors;
#elif defined(_SC_NPROCESSORS_ONLN)
        cores = sysconf(_SC_NPROCESSORS_ONLN);
        if (cores <= 0) cores = 1;
#endif

        if (cores > 1) cores = cores - 1;
        if (cores > 16) cores = 16;
        if (cores < 1) cores = 1;

        return cores;
}

struct parameters*init_param(void)
{
        struct parameters* param = NULL;
        MMALLOC(param, sizeof(struct parameters));
        param->infile = NULL;
        param->num_infiles = 0;
        param->input = NULL;
        param->outfile = NULL;
        param->format = NULL;
        param->help_flag = 0;
        param->dump_internal = 0;
        param->type = -1;
        param->gpo = -1.0;
        param->gpe = -1.0;
        param->tgpe = -1.0;
        param->nthreads = get_default_thread_count();
        param->min_support = 0;
        param->load_poar = NULL;
        param->mode = NULL;
        param->quiet = 0;
        param->out_format = 0;
        param->reformat = 0;
        param->rename = 0;
        param->clean = 0;
        param->unalign = 0;
        return param;
ERROR:
        free_parameters(param);
        return NULL;
}

void free_parameters(struct parameters* param)
{
        if(param){
                if(param->num_infiles){
                        MFREE(param->infile);
                }
                MFREE(param);
        }
}

int check_msa_format_string(char* format)
{
        int ok = 0;
        if(format){
                if(strstr(format,"msf")){
                        ok = 1;
                }else if(strstr(format,"clu")){
                        ok = 1;
                }else if(strstr(format,"fasta")){
                        ok = 1;
                }else if(strstr(format,"fa")){
                        ok = 1;
                }else{
                        ok = 0;
                        ERROR_MSG("Format %s not recognized.",format);
                }
                if(!ok){
                        ERROR_MSG("Format %s not recognized.",format);
                }
        }
        return OK;
ERROR:
        return FAIL;
}
