
#include "tldevel.h"
#include "rng.h"
#include <getopt.h>

#include "thr_pool.h"
#include "alphabet.h"
#include "hash_table.h"
#include "align_io.h"
#include "alignment_parameters.h"
#include "estimate_aln_param.h"
#include "tree_building.h"
#include "bisectingKmeans.h"
#include "alignment.h"
#include "weave_alignment.h"
#include "misc.h"

#define ALPHABET_LEN 21


struct jobs{
        char* in;
        char* param;
        int param_index;
        double SP;
        double TC;
        int job_number;
};

struct thread_data{
        struct jobs** jobs;
        struct batch_train* bt;
        int num_threads;
        int num_jobs;
        int id;
};

struct individual{
        unsigned int* code;
        double* param;
        double score;
};

struct pbil_data{
        struct individual** population;
        struct rng_state* rng;
        struct individual* best;
        double* bit_prob;
        double* min;
        double* max;
        double bit_mutation_prob;
        double bit_mutation_shift;
        int bits_per_param;
        double* step;
        int n_bits;
        int mu;                 /* number of samples */
        int lambda;             /* selected */
        int num_gen;
        int num_param;
        int num_threads;
        double gamma;           /* learning rate */
};
HT_GLOBAL_INIT(TESTINT, int*)

struct batch_train{
        struct alignment** aln; /* each element can be  */
        struct pbil_data* pbil;
        double** scores;        /* threads work on separate row */
        double* num_pairs;
        HT_TYPE(TESTINT )* ht;
        int num_alignments;
};

int fill_pair_hash(struct batch_train* bt);
/* Memory functions  */
struct pbil_data* init_pbil_data(int num_param,int num_bits,int num_gen,int sampled, int selected);
void free_pbil_data(struct pbil_data* d);

/* core pbil functions */
int sample_pop(struct pbil_data* d);
int write_kalign_parameter_files(struct pbil_data* d);
int update_pbil(struct pbil_data* d);
int mutate_prob_vector(struct pbil_data* d);

int print_best(struct pbil_data* d, char* out);


/* Misc  */
unsigned int BinaryToGray(unsigned int num);
unsigned int GrayToBinary32(unsigned int num);
int print_help(char **argv);
/* objective function */
int eval_batch(struct batch_train* bt,struct thr_pool* pool,int num_threads);
int eval(struct pbil_data* d,char** infile, int num_infiles,struct thr_pool* pool);
void* run_kalign_thread(void *threadarg);
void* run_kalign_batch_thread(void *threadarg);
int fill_sampled_parameters(struct aln_param* ap, double* values,int L);
int score_aln(struct alignment* aln, int aln_id,  HT_TYPE(TESTINT )* ht, double pairs,double* score);

int run_kalign_and_score(char* infile,char* p_file, double* SP, double* TC,char* name);

/* seeding */
int set_pbil_based_on_pop(struct pbil_data* d);
int init_pop_from_seed(struct pbil_data*d, char* infile);
int read_aln_param_from_file(struct pbil_data*d, char* infile);

/* for testing */
int random_score(struct pbil_data* d);

int sort_pop_by_score(const void *a, const void *b);


#define OPT_NGEN 1
#define OPT_MU 2
#define OPT_LAMBDA 3
#define OPT_NTHREADS 4

int main(int argc, char *argv[])
{
        int i;
        //double SP,TC;
        int num_param;

        struct batch_train* bt = NULL;
        struct pbil_data* d = NULL;

        struct thr_pool* pool;



        int num_infiles = 0;
        char** infile = NULL;
        char* outfile = NULL;
        char* seedfile = NULL;
        int c = 0;
        int n_gen,mu,lambda,num_threads,help;

        n_gen = 1000;
        mu = 100;
        lambda = 1;
        num_threads = 8;
        help = 0;
        print_program_header(argv, "Optimises alignment parameters.");
        //int help = 0;
        while (1){
                static struct option long_options[] ={
                        {"ngen",  required_argument, 0, OPT_NGEN},
                        {"popsize",  required_argument, 0, OPT_MU},
                        {"keep",  required_argument, 0, OPT_LAMBDA},
                        {"nthreads",  required_argument, 0, OPT_NTHREADS},
                        {"seed",  required_argument, 0, 's'},
                        {"out",  required_argument, 0, 'o'},
                        {"help",0,0,'h'},
                        {0, 0, 0, 0}
                };
                int option_index = 0;
                c = getopt_long_only (argc, argv,"hao:s:",long_options, &option_index);

                if (c == -1){
                        break;
                }
                switch(c) {
                case OPT_NGEN:
                        n_gen= atoi(optarg);
                        break;
                case OPT_MU:
                        mu = atoi(optarg);
                        break;
                case OPT_LAMBDA:
                        lambda = atoi(optarg);
                        break;
                case OPT_NTHREADS:
                        num_threads = atoi(optarg);
                        break;
                case 's':
                        seedfile = optarg;
                        break;
                case 'o':
                        outfile = optarg;
                        break;
                case 'h':

                        help = 1;
                        break;
                default:
                        ERROR_MSG("not recognized");
                        break;
                }
        }

        if(help){
                print_help(argv);
                return EXIT_SUCCESS;
        }

        if(!outfile){
                print_help(argv);
                ERROR_MSG("no outfile suffix - use -o <>");

        }
        num_infiles = 0;
        if (optind < argc){

                //fprintf(stderr,"EXTRA :%d\n",argc - optind);
                num_infiles = argc-optind;
                MMALLOC(infile, sizeof(char*) * num_infiles);
                c = 0;
                while (optind < argc){
                        infile[c] =  argv[optind++];
                        c++;
                }
        }


        if(!num_infiles){
                print_help(argv);
                ERROR_MSG("No input files found");

        }



        RUNP(pool = thr_pool_create(num_threads, num_threads, 0, 0));

        num_param = (ALPHABET_LEN * (ALPHABET_LEN-1)) / 2 + ALPHABET_LEN + 3;

        LOG_MSG("Starting run with parameters:");

        LOG_MSG("%d generations.", n_gen);
        LOG_MSG("%d popsize.", mu);
        LOG_MSG("%d keep.", lambda);
        LOG_MSG("%d threads.", num_threads);



        MMALLOC(bt, sizeof(struct batch_train));
        bt->num_alignments = num_infiles;
        bt->num_pairs = NULL;
        bt->aln = NULL;
        bt->scores = NULL;
        bt->scores = galloc(bt->scores,num_infiles, mu,0.0);
        MMALLOC(bt->aln, sizeof(struct alignment*)* bt->num_alignments);
        MMALLOC(bt->num_pairs, sizeof(double) * bt->num_alignments);
        for(i = 0; i < bt->num_alignments;i++){
                RUNP(bt->aln[i] = read_alignment(infile[i]));
                bt->num_pairs[i] = 0.0;
        }

        LOG_MSG("Read all alignments");
        RUN(fill_pair_hash(bt));


        for(i = 0; i < bt->num_alignments;i++){
                RUN(dealign(bt->aln[i]));
                RUN(convert_alignment_to_internal(bt->aln[i], defPROTEIN));

        }



        RUNP(d = init_pbil_data(num_param, 16, n_gen, mu, lambda));
        d->num_threads = num_threads;

        if(seedfile){
                RUN(init_pop_from_seed(d, seedfile));
        }
        bt->pbil = d;
        //d->lambda = 5;

        for(i = 0; i < d->num_gen;i++){
                RUN(sample_pop(d));

                RUN(write_kalign_parameter_files(d));
                eval_batch(bt, pool, d->num_threads);
//                exit(0);
                //RUN(eval(d, infile, num_infiles, pool));
                //random_score(d);
                RUN(update_pbil(d));
                RUN(mutate_prob_vector(d));
                RUN(print_best(d, outfile));
        }



        //RUN(run_kalign_and_score("~/data/bb3_release/RV30/BB30001.xml","test_kalign_param.txt", &SP, &TC));
        //fprintf(stdout,"%f %f\n",SP,TC);
        for(i = 0; i < bt->num_alignments;i++){
                free_aln(bt->aln[i]);


        }

        HT_FREE(TESTINT, bt->ht);

        MFREE(bt->aln);
        MFREE(bt->num_pairs);
        gfree(bt->scores);
        MFREE(bt);
        free_pbil_data(d);

        thr_pool_destroy( pool);

        MFREE(infile);
        return EXIT_SUCCESS;
ERROR:
        if(infile){
                MFREE(infile);
        }
        return EXIT_FAILURE;
}

int print_help(char **argv)
{
        const char usage[] = "  < *.xml | *.msf > -o <outfile suffix>";
        fprintf(stdout,"\nUsage: %s [-options] %s\n\n",basename(argv[0]) ,usage);
        fprintf(stdout,"Options:\n\n");
        fprintf(stdout,"%*s%-*s: %s %s\n",3,"",MESSAGE_MARGIN-3,"--ngen","Number of generations." ,"[1000]"  );
        fprintf(stdout,"%*s%-*s: %s %s\n",3,"",MESSAGE_MARGIN-3,"--popsize","Size of sampled population." ,"[100]"  );
        fprintf(stdout,"%*s%-*s: %s %s\n",3,"",MESSAGE_MARGIN-3,"--keep","Number of best solutions to keep." ,"[1]"  );
        fprintf(stdout,"%*s%-*s: %s %s\n",3,"",MESSAGE_MARGIN-3,"--nthreads","Number of threads." ,"[8]"  );
        fprintf(stdout,"%*s%-*s: %s %s\n",3,"",MESSAGE_MARGIN-3,"--seed","File containing a starting solution." ,"[]"  );

        return OK;
}

int eval_batch(struct batch_train* bt,struct thr_pool* pool,int num_threads)
{
        struct thread_data** td = NULL;
        int i,j;
        int status;

        /* clear scores */
        for ( i = 0; i < bt->num_alignments;i++){
                for(j = 0; j < bt->pbil->mu;j++){
                        bt->scores[i][j]= 0.0;
                }
        }

        MMALLOC(td, sizeof(struct thread_data*)* num_threads);

        for(i = 0; i < num_threads;i++){
                td[i] = NULL;
                MMALLOC(td[i],sizeof(struct thread_data));
                td[i]->id = i;
                td[i]->bt = bt;
                //td[i]->jobs = jobs;
                //td[i]->num_jobs = num_jobs;
                td[i]->num_threads = num_threads;
        }
        for(i = 0; i < num_threads ;i++){
                if((status = thr_pool_queue(pool, run_kalign_batch_thread, td[i])) == -1) ERROR_MSG("Adding job to queue failed.");
        }
        thr_pool_wait(pool);


        for(i = 0; i < bt->pbil->mu;i++){
                bt->pbil->population[i]->score = 0.0;
                for(j = 0; j < bt->num_alignments;j++){

                        bt->pbil->population[i]->score += bt->scores[j][i];
                }
        }

        for(i = 0; i < bt->pbil->mu;i++){
                bt->pbil->population[i]->score /= (double)  bt->num_alignments;
                //fprintf(stdout,"SCORE: %d %f\n",i,d->population[i]->score);
        }



        for(i = 0; i < num_threads;i++){

                MFREE(td[i]);//,sizeof(struct thread_data));
        }

        MFREE(td);
        return OK;
ERROR:
        return FAIL;
}


void* run_kalign_batch_thread(void *threadarg)
{
        struct thread_data *data;
        struct alignment* aln = NULL;
        struct aln_param* ap = NULL;
        int** map = NULL;       /* holds all alignment paths  */
        int num_jobs;
        int id;
        int i;
        int j;
        int c;

        data = (struct thread_data *) threadarg;
        ASSERT(data != NULL, "No data");

        id = data->id;


        for(i = 0; i < data->bt->num_alignments;i++){

                if(i % data->num_threads == id){
                        aln = data->bt->aln[i];

                        for(j = 0; j < data->bt->pbil->mu;j++){
                                /* dealign */
                                dealign(aln);

                                /* run kalign */
                                RUNP(ap = init_ap(aln->numseq,aln->L));
                                RUN(fill_sampled_parameters(ap, data->bt->pbil->population[j]->param, aln->L));
                                RUN(build_tree_kmeans(aln,ap));
                                RUNP(map = hirschberg_alignment(aln, ap));
                                RUN(weave(aln , map, ap->tree));

                                RUN(score_aln(aln, i, data->bt->ht, data->bt->num_pairs[i], &data->bt->scores[i][j]));
                                //fprintf(stdout,"%f ",  data->bt->scores[i][j]);

                                for(c = 0; c < aln->num_profiles ;c++){
                                        if(map[c]){
                                                MFREE(map[c]);
                                        }
                                }

                                free_ap(ap);
                                MFREE(map);

                                map = NULL;

                                int c;
                                for( c = aln->numseq;c < aln->num_profiles;c++){
                                        MFREE(aln->sip[c]);
                                        aln->sip[c] = NULL;
                                }

                        }
                        //fprintf(stdout,"\n");
                }
        }

        return NULL;
ERROR:
        return NULL;
}

int fill_sampled_parameters(struct aln_param* ap, double* values,int L)
{

        int i,j;

        int  m_pos = 0;
        for (i = 0;i < L;i++){
                for (j = 0;j <= i;j++){
                        ap->subm[i][j] = values[m_pos];
                        ap->subm[j][i] = ap->subm[i][j];
                        m_pos++;
                }
        }

        ap->gpo = values[m_pos];
        m_pos++;
        ap->gpe = values[m_pos];
        m_pos++;

        ap->tgpe = values[m_pos];
        return OK;
}

int score_aln(struct alignment* aln, int aln_id,  HT_TYPE(TESTINT )* ht, double pairs,double* score)
{

        uint8_t* aligned_a = NULL;
        uint8_t* aligned_b = NULL;

        int numseq;
        int i,j,c,f,g,res;

        int aln_len = 0;

        double common = 0.0;

        int* tmp = NULL;
        hash_table_node_TESTINT_t* hashnode = NULL;

        tmp = NULL;
        tmp= galloc(tmp,5);

        numseq = aln->numseq;

        for(i = 0; i <= aln->sl[0];i++){
                aln_len += aln->gaps[0][i];
        }
        aln_len += aln->sl[0];
        //LOG_MSG("Aln len: %d.",aln_len);

        MMALLOC(aligned_a, sizeof(uint8_t) * aln_len);
        MMALLOC(aligned_b, sizeof(uint8_t) * aln_len);



        for(i = 0; i < numseq;i++){
                RUN(make_aliged_seq(aligned_a, aln->s[i], aln->gaps[i], aln->sl[i]));
                for(j = i+1;j < numseq;j++){
                        RUN(make_aliged_seq(aligned_b, aln->s[j], aln->gaps[j], aln->sl[j]));

                        f = 0;
                        g = 0;
                        for(c = 0;c < aln_len;c++){
                                res= 0;
                                if(aligned_a[c] != 255){
                                        res |= 1;
                                }
                                if(aligned_b[c] != 255){
                                        res |=2;
                                }
                                switch(res){
                                case 3:
                                        tmp[0] = aln_id;
                                        tmp[1] = i;
                                        tmp[2] = j;
                                        tmp[3] = f;
                                        tmp[4] = g;
                                        hashnode = HT_SEARCH(TESTINT, ht, tmp);
                                        if(hashnode){
                                                common += 1.0;
                                        }

                                        f++;
                                        g++;
                                        break;
                                case 2:
                                        g++;
                                        break;
                                case 1:
                                        f++;
                                        break;
                                default:
                                        break;
                                }


                        }
                }
        }
        gfree(tmp);
        //MFREE(tmp);
        MFREE(aligned_a);
        MFREE(aligned_b);
//        LOG_MSG("Alignment: %d , %f %f\n", aln_id,common,pairs);
        *score = common /pairs;
        return OK;
ERROR:
        return FAIL;
}


int fill_pair_hash(struct batch_train* bt)
{
        struct alignment* aln = NULL;
        int i,j,c,f,g,a,res;
        int* tmp;

        bt->ht = NULL;
        bt->ht = HT_INIT(TESTINT,63556);
        for(i = 0; i < bt->num_alignments;i++){
                aln = bt->aln[i];
                bt->num_pairs[i] = 0.0;
                for(j = 0; j < aln->numseq;j++){
                        for(c = j+1;c < aln->numseq;c++){
                                f = 0;
                                g = 0;
                                for(a = 0; a < aln->sl[0];a++){
                                        res= 0;
                                        if(isalpha(aln->seq[j][a])){
                                                res |= 1;
                                        }
                                        if(isalpha(aln->seq[c][a])){
                                                res |=2;
                                        }
                                        switch(res){
                                        case 3:
                                                tmp = NULL;
                                                tmp= galloc(tmp,5);
                                                tmp[0] = i;
                                                tmp[1] = j;
                                                tmp[2] = c;
                                                tmp[3] = f;
                                                tmp[4] = g;
                                                RUN(HT_INSERT(TESTINT,bt->ht,tmp,NULL));
                                                bt->num_pairs[i] += 1.0;
                                                f++;
                                                g++;
                                                break;
                                        case 2:
                                                g++;
                                                break;
                                        case 1:
                                                f++;
                                                break;
                                        default:
                                                break;
                                        }

                                }

                        }
                }

        }
/*        HT_PRINT(TESTINT,bt->ht);
        HT_FREE(TESTINT,bt->ht);
        exit(0);*/
        return OK;
ERROR:
        return FAIL;
}

int eval(struct pbil_data* d,char** infile, int num_infiles,struct thr_pool* pool)
{

        struct jobs** jobs = NULL;
        struct thread_data** td = NULL;
        char** param_file_name_buffer = NULL;
        int i,j,c;
        int num_jobs;
        int status;


        /* fill parameter namesl */
        param_file_name_buffer = galloc(param_file_name_buffer,d->mu,BUFFER_LEN,0);
        for(i = 0; i < d->mu;i++){
                snprintf(param_file_name_buffer[i] , BUFFER_LEN, "kalign_param_%d.txt", i+1);
                d->population[i]->score =0.0;
        }

        num_jobs = num_infiles * d->mu;

        MMALLOC(jobs, sizeof(struct jobs*) * num_jobs);
        for(i = 0; i < num_jobs;i++){
                jobs[i] = NULL;
                MMALLOC(jobs[i], sizeof(struct jobs));
                jobs[i]->SP = 0.0;
                jobs[i]->TC = 0.0;
                jobs[i]->in = NULL;
                jobs[i]->param = NULL;
                jobs[i]->job_number = i;
        }
        c = 0;
        for(i = 0; i < num_infiles;i++){

                for(j = 0; j < d->mu;j++){
                        jobs[c]->in = infile[i];
                        jobs[c]->param = param_file_name_buffer[j];
                        jobs[c]->param_index = j;
                        c++;
                }
        }

        MMALLOC(td, sizeof(struct thread_data*)* d->num_threads);

        for(i = 0; i < d->num_threads;i++){
                td[i] = NULL;
                MMALLOC(td[i],sizeof(struct thread_data));
                td[i]->id = i;
                td[i]->jobs = jobs;
                td[i]->num_jobs = num_jobs;
                td[i]->num_threads = d->num_threads;
        }

        for(i = 0; i < d->num_threads ;i++){
                if((status = thr_pool_queue(pool, run_kalign_thread, td[i])) == -1) ERROR_MSG("Adding job to queue failed.");
        }
        thr_pool_wait(pool);
        for(i = 0; i < num_jobs;i++){
                //fprintf(stdout,"%s %s: %f %f\n",td[i]->in,td[i]->param,td[i]->SP,td[i]->TC);
                d->population[jobs[i]->param_index]->score += jobs[i]->SP;// td[i]->SP;
        }
        for(i = 0; i < d->mu;i++){
                d->population[i]->score /= (double) num_infiles;
                //fprintf(stdout,"SCORE: %d %f\n",i,d->population[i]->score);
        }

        for(i = 0; i < num_jobs;i++){
                MFREE(jobs[i]);
        }
        MFREE(jobs);
        for(i = 0; i < d->num_threads;i++){
                MFREE(td[i]);
        }
        MFREE(td);
        gfree(param_file_name_buffer);
        return OK;
ERROR:
        return FAIL;

}

void* run_kalign_thread(void *threadarg)
{
        char buffer[32];

        struct jobs** jobs = NULL;
        struct thread_data *data;
        int num_jobs;
        int id;
        int i;
        data = (struct thread_data *) threadarg;
        ASSERT(data != NULL, "No data");

        jobs = data->jobs;
        num_jobs = data->num_jobs;
        id = data->id;


        for(i = 0; i < num_jobs;i++){
                if(i % data->num_threads == id){
                        snprintf(buffer, 32, "job%d.msf",  jobs[i]->job_number);
                //fprintf(stdout,"%s %s\n", data->in,data->param);
                        RUN(run_kalign_and_score(jobs[i]->in, jobs[i]->param, &jobs[i]->SP, &jobs[i]->TC, buffer));
                }
        }
        return NULL;
ERROR:
        return NULL;
}


int write_kalign_parameter_files(struct pbil_data* d)
{
        FILE* f_ptr = NULL;
        char buffer[BUFFER_LEN];
        int i,j;
        for(i = 0; i < d->mu;i++){
                snprintf(buffer, BUFFER_LEN, "kalign_param_%d.txt", i+1);
                //            LOG_MSG("Open %s",buffer);
                RUNP( f_ptr = fopen(buffer, "w"));
                for(j =0; j < d->num_param;j++){
                        fprintf(f_ptr,"%f\n",  d->population[i]->param[j]);
                }
                fclose(f_ptr);

        }
        return OK;
ERROR:
        return FAIL;
}

int set_pbil_based_on_pop(struct pbil_data* d)
{
        int i,j,c,f;
        double sum;

        qsort(d->population, d->mu, sizeof(struct individual*), sort_pop_by_score);
        sum = 0.0;
        for(i = 0;i < d->lambda;i++){
                fprintf(stdout,"%f ", d->population[i]->score);
                sum += d->population[i]->score;
        }
        fprintf(stdout,"average:%f\n",sum / (double) d->lambda);
        f = 0;
        for(j= 0 ;j < d->num_param;j++){

                for(c = 0;c < d->bits_per_param;c++){
                        sum = 0;
                        for(i = 0; i < d->lambda;i++){
                                if(d->population[i]->code[j] & (1 << c)){


                                        sum++;
                                }
                        }
                        sum = sum / (double) d->lambda;


                        d->bit_prob[f] = sum;

                        //fprintf(stdout,"%d %f\n",f, d->bit_prob[f]);
                        f++;
                }
        }

        return OK;
}

int update_pbil(struct pbil_data* d)
{
        int i,j,c,f;
        double sum;
        qsort(d->population, d->mu, sizeof(struct individual*), sort_pop_by_score);
        sum = 0.0;

        for(i = 0;i < d->lambda;i++){
                //fprintf(stdout,"%0.2f ", d->population[i]->score);
                sum += d->population[i]->score;
        }

        //fprintf(stdout,"average:%f\n",sum / (double) d->lambda);
        if(d->population[0]->score > d->best->score){

                for(i = 0; i < d->num_param;i++){
                        d->best->param[i] = d->population[0]->param[i];
                }
                d->best->score =d->population[0]->score;
                LOG_MSG("Found new best: %f.",d->best->score);
        }
        f = 0;
        for(j= 0 ;j < d->num_param;j++){

                for(c = 0;c < d->bits_per_param;c++){
                        sum = 0.0;
                        for(i = 0; i < d->lambda;i++){
                                if(d->population[i]->code[j] & (1 << c)){


                                        sum++;
                                }
                        }
                        sum = sum / (double) d->lambda;
                        sum = sum * d->gamma;

                        d->bit_prob[f] = d->bit_prob[f] * (1.0 - d->gamma) + sum;

                        //fprintf(stdout,"%d %f\n",f, d->bit_prob[f]);
                        f++;
                }
        }
        return OK;
}


int print_best(struct pbil_data* d,char* out)
{
        char buffer[BUFFER_LEN];
        FILE* f_ptr = NULL;
        int i,j,c;


        ASSERT(d != NULL, "No data");

        snprintf(buffer, BUFFER_LEN, "%s_flat.txt", out);

        //LOG_MSG("Open %s",buffer);
        RUNP( f_ptr = fopen(buffer, "w"));
        for( i = 0; i < d->num_param;i++){

                fprintf(f_ptr,"%f\n", d->best->param[i]);
        }
        fflush(f_ptr);
        fclose(f_ptr);

        snprintf(buffer, BUFFER_LEN, "%s_arr.txt", out);

        //LOG_MSG("Open %s",buffer);
        RUNP( f_ptr = fopen(buffer, "w"));
        fprintf(f_ptr,"float balimt[]={\n");
        c =0;
        for(i = 0; i < ALPHABET_LEN;i++){
                //fprintf(stdout,"%d",i);
                for(j = 0; j <= i;j++){


                        fprintf(f_ptr," %f,",d->best->param[c]);
                        c++;
                }
                fprintf(f_ptr,"\n");
        }
        fprintf(f_ptr,"};\n");

        fprintf(f_ptr,"ap->gpo = %f;\n", d->best->param[c]);
        c++;

        fprintf(f_ptr,"ap->gpe =  %f;\n", d->best->param[c]);
        c++;

        fprintf(f_ptr,"ap->tgpe =  %f;\n", d->best->param[c]);
        c++;
        fflush(f_ptr);
        fclose(f_ptr);


        return OK;
ERROR:
        return FAIL;
}

int random_score(struct pbil_data* d)
{
        int i;
        double r;
        for(i = 0; i < d->mu;i++){
                r =  tl_random_double(d->rng);
                //RUN(drand48_r(&d->randBuffer, &r));
                d->population[i]->score = r;

        }

        return OK;
}

int sample_pop(struct pbil_data* d)
{
        int i,j;
        double r;
        for(i = 0; i < d->mu;i++){
                for(j = 0; j < d->num_param;j++){
                        d->population[i]->code[j] = 0;
                }

                for(j = 0; j < d->n_bits;j++){
                        r = tl_random_double(d->rng);

                        //fprintf(stdout,"%f",r);
                        if(r <= d->bit_prob[j]){
                                d->population[i]->code[j/d->bits_per_param] |= 1 << (j % d->bits_per_param);
                        }

                }
                for(j = 0; j < d->num_param;j++){
                        d->population[i]->code[j] = GrayToBinary32(d->population[i]->code[j]);
                        d->population[i]->param[j] =  d->min[j] + d->step[j] * d->population[i]->code[j];

                        //fprintf(stdout,"%f ", d->population[i]->param[j]);
                }

                //fprintf(stdout," %d\n",d->num_param);
        }
        return OK;
}




int mutate_prob_vector(struct pbil_data* d)
{
        int i;
        double r;
        for(i = 0; i < d->n_bits;i++){
                r = tl_random_double(d->rng);
                //drand48_r(&d->randBuffer, &r);
                if(r < d->bit_mutation_prob){
                        r = tl_random_double(d->rng);

                        d->bit_prob[i] = d->bit_prob[i] * (1.0 - d->bit_mutation_shift) + r * d->bit_mutation_shift;
                }


                if(d->bit_prob[i] <= 0.0){
                        d->bit_prob[i] = 0.00001;
                }
                if(d->bit_prob[i] >= 1.0){
                        d->bit_prob[i] = 0.99999;
                }
        }
        return OK;
}

struct pbil_data* init_pbil_data(int num_param,int num_bits,int num_gen,int sampled, int selected)
{
        struct pbil_data* d = NULL;
        int i = 0;
        MMALLOC(d, sizeof(struct pbil_data));
        d->bit_prob = NULL;
        d->gamma = 0.01;// 1.0 / (double) num_gen;
        d->lambda = selected;
        d->mu = sampled;
        d->num_gen = num_gen;
        d->num_param = num_param;
        d->population = NULL;
        d->n_bits = num_bits * num_param;
        d->bits_per_param = num_bits;
        d->step = NULL;
        d->min = NULL;
        d->max = NULL;
        d->best = NULL;
        d->bit_mutation_prob = 0.02;
        d->bit_mutation_shift = 0.05;
        MMALLOC(d->bit_prob, sizeof(double) * num_bits * num_param);
        MMALLOC(d->step, sizeof(double) * num_param);
        MMALLOC(d->min, sizeof(double) * num_param);
        MMALLOC(d->max, sizeof(double) * num_param);

        for(i =0; i < num_param;i++){

                d->min[i] = -10.0;
                d->max[i] = 10.0;

        }
        for(i =num_param-3; i < num_param;i++){

                d->min[i] = 0.0;
                d->max[i] = 10.0;

        }
        for(i =0; i < num_param;i++){
                d->step[i] = (1 << num_bits) -1;

                d->step[i] = (d->max[i] - d->min[i]) / d->step[i];

        }
        //srand48_r(time(NULL), &d->randBuffer);
        for(i = 0; i < d->n_bits;i++){
                d->bit_prob[i] = 0.5f;
        }


        MMALLOC(d->best, sizeof(struct individual));
        d->best->score = 0.0f;
        d->best->code = NULL;
        d->best->param = NULL;;
        MMALLOC(d->best->code, sizeof(unsigned int) * num_param);
        MMALLOC(d->best->param, sizeof(double) * num_param);

        MMALLOC(d->population,sizeof(struct individual*) * d->mu);
        for(i = 0; i < d->mu;i++){
                d->population[i] = NULL;
                MMALLOC(d->population[i], sizeof(struct individual));
                d->population[i]->score = 0.0f;
                d->population[i]->code = NULL;
                d->population[i]->param = NULL;;
                MMALLOC(d->population[i]->code, sizeof(unsigned int) * num_param);
                MMALLOC(d->population[i]->param, sizeof(double) * num_param);
        }

        RUNP(d->rng = init_rng(0));
        return d;
ERROR:
        free_pbil_data(d);
        return NULL;

}

void free_pbil_data(struct pbil_data* d)

{
        int i;
        if(d){
                if(d->population){
                        for (i = 0;i < d->mu;i++){
                                MFREE(d->population[i]->param);
                                MFREE(d->population[i]->code);
                                MFREE(d->population[i]);
                        }
                        MFREE(d->population);
                }
                MFREE(d->best->param);
                MFREE(d->best->code);
                MFREE(d->best);
                MFREE(d->max);
                MFREE(d->min);
                MFREE(d->step);
                MFREE(d->bit_prob);
                MFREE(d->rng);
                MFREE(d);
        }
}

int run_kalign_and_score(char* infile,char* p_file, double* SP, double* TC,char* name)
{

        FILE *pf;
        char cmd[BUFFER_LEN+32];
        char ret[BUFFER_LEN];

        snprintf(cmd, BUFFER_LEN+32, "run_kalign_bali_score.sh -a %s -i %s -n %s",p_file, infile,name);
        // Execute a process listing


        // Setup our pipe for reading and execute our command.
        RUNP(pf = popen(cmd,"r"));


        // Grab data from process execution


        while (fgets(ret, BUFFER_LEN, pf)){
                //fprintf(stdout,"%s", ret);
        }

        //fprintf(stdout,"-%s-", ret);

        sscanf(ret, "%*s %*s %lf %lf", SP,TC);
        pclose(pf);

        return OK;
ERROR:
        return FAIL;
}



int sort_pop_by_score(const void *a, const void *b)
{
        struct individual* const *one = a;
        struct individual* const *two = b;

        if((*one)->score >= (*two)->score){
                return -1;
        }else{
                return 1;
        }
}



int init_pop_from_seed(struct pbil_data*d, char* infile)
{
        int i,j;
        double r;
        double tmp;
        RUN(read_aln_param_from_file(d,infile));

        for(i = 0; i < d->num_param;i++){
                d->min[i] = d->population[0]->param[i] - 1.0;
                d->max[i] = d->population[0]->param[i] + 1.0;
                d->step[i] = (1 << d->bits_per_param) -1;

                d->step[i] = (d->max[i] - d->min[i]) / d->step[i];

        }
        for(i = d->num_param -3; i < d->num_param;i++){
                d->min[i] = d->population[0]->param[i] - 1.0;
                d->max[i] = d->population[0]->param[i] + 1.0;
                if(d->min[i] < 0.0){
                        d->min[i] = 0.0;
                }
                d->step[i] = (1 << d->bits_per_param) -1;

                d->step[i] = (d->max[i] - d->min[i]) / d->step[i];

        }

        for(i = 0; i < d->num_param;i++){
                tmp = d->population[0]->param[i] - d->min[i];
                tmp = tmp/ (d->max[i] - d->min[i]); /* should be number between 0 and 1 */
                d->population[0]->code[i] = (unsigned int) ((double) (1 << d->bits_per_param) * tmp);
                d->population[0]->code[i] = BinaryToGray(d->population[0]->code[i]);
                d->population[0]->code[i] = GrayToBinary32(d->population[0]->code[i]);
                //fprintf(stdout,"%d %f ", d->population[0]->code[i],d->population[0]->param[i]);
                d->population[0]->param[i] =  d->min[i] + d->step[i]* d->population[0]->code[i];
                //fprintf(stdout,"%f INIT\n", d->population[0]->param[i]);
        }
        //exit(0);
        for(i = 0; i < d->n_bits;i++){
                if(     d->population[0]->code[i/d->bits_per_param] & 1 << (i % d->bits_per_param)){
                        //c = i % d->bits_per_param;
                        d->bit_prob[i] = 0.99;

                }else{

                        //c = i % d->bits_per_param;
                        d->bit_prob[i] = 0.01;
                }

        }

        //

        i = 1;

        for(j = 0; j < d->num_param;j++){
                d->population[i]->code[j] = 0;
        }

        for(j = 0; j < d->n_bits;j++){
                r = tl_random_double(d->rng);
                //fprintf(stdout,"%f\n", d->bit_prob[j]);
                if(r <= d->bit_prob[j]){
                        d->population[i]->code[j/d->bits_per_param] |= 1 << (j % d->bits_per_param);
                }


        }

        for(j = 0; j < d->num_param;j++){
                d->population[i]->code[i] = BinaryToGray(d->population[i]->code[i]);
                d->population[i]->code[i] = GrayToBinary32(d->population[i]->code[i]);

                d->population[i]->param[j] =  d->min[j]+ d->step[j]* d->population[i]->code[j];

                //fprintf(stdout,"%d %f %f\n", d->population[i]->code[j], d->population[0]->param[j], d->population[1]->param[j]);
        }
        //exit(0);
        return OK;
ERROR:
        return FAIL;
}



unsigned int BinaryToGray(unsigned int num)
{
    return num ^ (num >> 1);
}

/*
 * A more efficient version for Gray codes 32 bits or fewer
 * through the use of SWAR (SIMD within a register) techniques.
 * It implements a parallel prefix XOR function.  The assignment
 * statements can be in any order.
 *
 * This function can be adapted for longer Gray codes by adding steps.
 * A 4-bit variant changes a binary number (abcd)2 to (abcd)2 ^ (00ab)2,
 * then to (abcd)2 ^ (00ab)2 ^ (0abc)2 ^ (000a)2.
 */
unsigned int GrayToBinary32(unsigned int num)
{
    num = num ^ (num >> 16);
    num = num ^ (num >> 8);
    num = num ^ (num >> 4);
    num = num ^ (num >> 2);
    num = num ^ (num >> 1);
    return num;
}
