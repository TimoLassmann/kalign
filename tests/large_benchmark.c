#include "tldevel.h"
#include "tlmisc.h"

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include <sys/stat.h>
#include <dirent.h>
#include <errno.h>

#include <unistd.h>
#include <getopt.h>

#include "kalign/kalign.h"

#define MAX_PATH_LEN 4096

struct parameters {
        char* target_dir;
        int verbose;
        int dryrun;
        int help_flag;
};

struct aln_case {
        char* dir;
        char* name;
        float score;
        int alloc_len;
        int valid;
};

struct bench_collection {
        struct aln_case** l;
        int n;
        int n_alloc;
};


static int run_test_aln(struct aln_case *tcase);
static int bench_collection_add_case(struct bench_collection* b, char* dir, char* name);
static int bench_collection_alloc(struct bench_collection **bench_collection, int n);
static int bench_collection_resize(struct bench_collection *b);
static void bench_collection_free(struct bench_collection *b);


static int aln_case_alloc(struct aln_case **aln_case);
static void aln_case_free(struct aln_case *a);

static int run_bench(struct parameters *p);
static int process_dir(struct parameters *p,char* dir,struct bench_collection* b);


static int detect_likely_alignment(char *name, int *is_aln);

static int param_sanity_check(struct parameters *p);
static int print_help(char * argv[]);
static int param_print(struct parameters *p);
static int param_init(struct parameters **param);
static void param_free(struct parameters *p);

static int local_dir_exists(const char* name);

#define OPT_VERBOSE 1
#define OPT_DRY 2

int main(int argc, char *argv[])
{

        struct parameters* param = NULL;
        int c;

        RUN(param_init(&param));

        while (1){
                static struct option long_options[] ={
                        {"dir",  required_argument, 0, 'd'},
                        {"verbose",   no_argument,0, OPT_VERBOSE},
                        {"dry-run",   no_argument,0,OPT_DRY},
                        {"help",   no_argument,0,'h'},
                        {0, 0, 0, 0}
                };

                int option_index = 0;

                c = getopt_long_only (argc, argv,"d:h",long_options, &option_index);

                /* Detect the end of the options. */
                if (c == -1){
                        break;
                }
                switch(c) {
                case OPT_DRY:
                        param->dryrun = 1;
                        break;
                case OPT_VERBOSE:
                        param->verbose = 1;
                        break;
                case 'd':
                        param->target_dir = optarg;
                        break;
                        case 'h':
                        param->help_flag = 1;
                        break;
                default:
                        abort ();
                }
        }
        if(param->help_flag || argc < 2){
                RUN(print_help(argv));
                param_free(param);
                return EXIT_SUCCESS;
        }

        RUN(param_sanity_check(param));
        RUN(param_print(param));

        RUN(run_bench(param));
        param_free(param);

        return EXIT_SUCCESS;
ERROR:
        param_free(param);
        return EXIT_FAILURE;
}

int run_bench(struct parameters *p)
{
        struct bench_collection* b = NULL;
        ASSERT(p != NULL,"No parameters");
        RUN(bench_collection_alloc(&b, 1024));
        RUN(process_dir(p, p->target_dir,b));

        for(int i = 0; i < b->n;i++){
                struct aln_case *c = NULL;
                c = b->l[i];
                if(p->dryrun){
                        LOG_MSG("DryRun: would run kalign on: %s",c->name );
                }else{
                        RUN(run_test_aln(c));

                        if(p->verbose){
                                LOG_MSG("Case %d:\t%-*s\t%d\t%f ",i, 20,c->name,c->valid, c->score);
                        }
                }
        }
        fprintf(stdout,"Alignment,Score\n");
        for(int i = 0; i < b->n;i++){
                struct aln_case *c = NULL;
                c = b->l[i];
                if(c->valid){
                        fprintf(stdout,"%s,%f\n",c->name,c->score);
                }
        }
        bench_collection_free(b);
        return OK;
ERROR:
        bench_collection_free(b);
        return FAIL;
}

int run_test_aln(struct aln_case *tcase)
{
        char path[MAX_PATH_LEN];
        struct msa* r = NULL;
        struct msa* t = NULL;

        snprintf(path, MAX_PATH_LEN, "%s/%s",tcase->dir,tcase->name);

        RUN(kalign_read_input(path, &t,1));

        if(!t){
                kalign_free_msa(t);
                tcase->valid = 0;
                return OK;
        }
        /* reference  */
        RUN(kalign_read_input(path, &r,1));
        kalign_run(t,16 , -1, -1, -1 , -1);



        kalign_msa_compare(r, t, &tcase->score);
        /* if(tcase->score >= 0.5){ */
        /*         LOG_MSG("Ref:"); */
        /*          kalign_write_msa(r, NULL , "clu"); */
        /*          LOG_MSG("test:"); */
        /*          kalign_write_msa(t, NULL , "clu"); */
        /* } */
        kalign_free_msa(r);
        kalign_free_msa(t);

        return OK;
ERROR:
        kalign_free_msa(r);
        kalign_free_msa(t);
        return FAIL;
}

int process_dir(struct parameters *p,char* dir,struct bench_collection* b)
{
        struct dirent *dp;
        DIR *d_ptr;
        if ((d_ptr = opendir(dir)) == NULL){
                ERROR_MSG("Can't open %s", dir);
        }

        while ((dp = readdir(d_ptr)) != NULL){
                /* struct stat stbuf ; */
                int is_aln = 0;
                if(dp->d_type == DT_DIR){
                        if(dp->d_name[0] != '.'){
                                char* tmp = NULL;
                                MMALLOC(tmp, sizeof(char) * MAX_PATH_LEN);
                                snprintf(tmp,MAX_PATH_LEN,"%s/%s", dir, dp->d_name);
                                RUN(process_dir(p, tmp,b));
                                MFREE(tmp);
                        }
                }else if(dp->d_type == DT_REG){
                        /* LOG_MSG("File: %s   %s",dir, dp->d_name); */
                        char* tmp = NULL;
                        MMALLOC(tmp, sizeof(char) * MAX_PATH_LEN);
                        snprintf(tmp,MAX_PATH_LEN,"%s/%s", dir, dp->d_name);
                        detect_likely_alignment(tmp, &is_aln);
                        if(is_aln){
                                RUN(bench_collection_add_case(b, dir, dp->d_name));
                        }
                        MFREE(tmp);
                }else if(dp->d_type == DT_LNK){
                        /* LOG_MSG("LINK: %s   %s",dir, dp->d_name); */
                        char* tmp = NULL;
                        char* rp = NULL;
                        MMALLOC(tmp, sizeof(char) * MAX_PATH_LEN);
                        snprintf(tmp,MAX_PATH_LEN,"%s/%s", dir, dp->d_name);
                        rp = realpath(tmp, NULL);
                        if(rp!= NULL){
                                /* LOG_MSG("resolved: %s", rp); */

                                detect_likely_alignment(rp, &is_aln);
                                if(is_aln){
                                        RUN(bench_collection_add_case(b, dir, dp->d_name));
                                }
                                free(rp);
                        }
                        MFREE(tmp);
                }
        }
        closedir(d_ptr);
        return OK;
ERROR:
        if(d_ptr){
                closedir(d_ptr);
        }
        return FAIL;
}

int detect_likely_alignment(char *name, int *is_aln)
{
        int len = strnlen(name,MAX_PATH_LEN);
        int ret = 0;
        if(strncmp(name + len - 4, ".msf", 4) == 0){
                /* LOG_MSG("%s is aln", name); */
                ret = 1;
        }else if(strncmp(name + len - 4, ".vie", 4) == 0){
                /* LOG_MSG("%s is aln", name); */
                ret = 1;
        }

        *is_aln = ret;
        return OK;
}


int print_help(char * argv[])
{
        const char usage[] = " -d <target dir> ";
        char* basename = NULL;


        RUN(tlfilename(argv[0], &basename));

        fprintf(stdout,"\nUsage: %s %s\n\n",basename ,usage);
        fprintf(stdout,"Options:\n\n");

        fprintf(stdout,"%*s%-*s: %s %s\n",3,"",MESSAGE_MARGIN-3,"-d","Target directory" ,"[]"  );
        fprintf(stdout,"%*s%-*s: %s %s\n",3,"",MESSAGE_MARGIN-3,"-verbose","Target directory" ,"[off]"  );
        fprintf(stdout,"%*s%-*s: %s %s\n",3,"",MESSAGE_MARGIN-3,"-dry-run","Target directory" ,"[off]"  );
        fprintf(stdout,"%*s%-*s: %s %s\n",3,"",MESSAGE_MARGIN-3,"-h/--help","Print help." ,"[off]"  );
        fprintf(stdout,"\n\n");
        if(basename){
                MFREE(basename);
        }
        return OK;
ERROR:
        if(basename){
                MFREE(basename);
        }
        return FAIL;
}


int param_print(struct parameters *p)
{

        fprintf(stderr,"%*s%-*s: %s\n",3,"",MESSAGE_MARGIN-3,"-d", p->target_dir);
        fprintf(stderr,"%*s%-*s: %d\n",3,"",MESSAGE_MARGIN-3,"-verbose", p->verbose);
        fprintf(stderr,"%*s%-*s: %d\n",3,"",MESSAGE_MARGIN-3,"-dry-run", p->dryrun);
        return OK;

}

int param_sanity_check(struct parameters *p)
{
        ASSERT(p->target_dir != NULL, " No target dir use -d .... ");
        RUN(local_dir_exists(p->target_dir));
        return OK;
ERROR:
        return FAIL;
}

int local_dir_exists(const char* name)
{
        ASSERT(name != NULL," No directory name!");
        DIR* dir = opendir(name);
        if (dir) {
                /* Directory exists. */
                closedir(dir);
        } else if (ENOENT == errno) {
                ERROR_MSG("Directory %s does not exist!", name);
                /* Directory does not exist. */
        } else {
                ERROR_MSG("Opening directory %s failed!", name);
                /* opendir() failed for some other reason. */
        }
        return OK;
ERROR:
        return FAIL;
}


int param_init(struct parameters **param)
{
        struct parameters* p = NULL;
        MMALLOC(p, sizeof(struct parameters));
        p->target_dir = NULL;
        p->verbose = 0;
        p->dryrun = 0;
        p->help_flag = 0;
        *param = p;

        return OK;
ERROR:
        param_free(p);
        return FAIL;
}

void param_free(struct parameters *p)
{
        if(p){
                MFREE(p);
        }
}


int bench_collection_add_case(struct bench_collection* b, char* dir, char* name)
{


        struct aln_case* a = NULL;
        int idx = 0;

        idx = b->n;
        a = b->l[idx];

        strncpy(a->name, name, a->alloc_len);
        strncpy(a->dir, dir, a->alloc_len);

        b->n++;

        if(b->n == b->n_alloc){
                RUN(bench_collection_resize(b));
        }
        return OK;
ERROR:
        return FAIL;
}


int bench_collection_alloc(struct bench_collection **bench_collection, int n)
{
        struct bench_collection* b = NULL;

        MMALLOC(b, sizeof(struct bench_collection));
        b->n = 0;
        b->n_alloc = n;
        b->l = NULL;

        MMALLOC(b->l, sizeof(struct aln_cases*) * b->n_alloc);

        for(int i = 0; i < b->n_alloc;i++){
                b->l[i] = NULL;
                RUN(aln_case_alloc(&b->l[i]));
        }
        *bench_collection = b;
        return OK;
ERROR:
        return FAIL;
}

int bench_collection_resize(struct bench_collection *b)
{

        int old_size = b->n_alloc;
        int new_size = b->n_alloc + b->n_alloc / 2;

        MREALLOC(b->l, sizeof(struct aln_cases*) * new_size);

        for(int i = old_size; i < new_size;i++){
                b->l[i] = NULL;
                RUN(aln_case_alloc(&b->l[i]));
        }

        b->n_alloc = new_size;
        return OK;
ERROR:
        return FAIL;
}


void bench_collection_free(struct bench_collection *b)
{
        if(b){
                for(int i = 0; i < b->n_alloc;i++){
                        aln_case_free(b->l[i]);
                }
                MFREE(b->l);
                MFREE(b);
        }
}

int aln_case_alloc(struct aln_case **aln_case)
{
        struct aln_case* a = NULL;
        MMALLOC(a, sizeof(struct aln_case));
        a->name = NULL;
        a->dir = NULL;
        a->alloc_len = MAX_PATH_LEN;
        a->valid  = 1;
        MMALLOC(a->dir, sizeof(char) * a->alloc_len);
        MMALLOC(a->name, sizeof(char) * a->alloc_len);
        a->score = -1.0;
        *aln_case = a;
        return OK;
ERROR:
        aln_case_free(a);
        return FAIL;
}

void aln_case_free(struct aln_case *a)
{
        if(a){
                if(a->dir){
                        MFREE(a->dir);
                }
                if(a->name){
                        MFREE(a->name);
                }
                MFREE(a);
        }

}
