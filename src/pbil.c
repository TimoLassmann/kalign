

#include "tldevel.h"

#include <getopt.h>

#include "thr_pool.h"




struct thread_data{
        char* in;
        char* param;
        int param_index;
        double SP;
        double TC;
};

struct individual{
        unsigned int* code;
        double* param;
        double score;
};

struct pbil_data{
        struct individual** population;
        struct drand48_data randBuffer;
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
        double gamma;           /* learning rate */
};

/* Memory functions  */
struct pbil_data* init_pbil_data(int num_param,int num_bits,int num_gen,int sampled, int selected);
void free_pbil_data(struct pbil_data* d);

/* core pbil functions */
int sample_pop(struct pbil_data* d);
int write_kalign_parameter_files(struct pbil_data* d);
int update_pbil(struct pbil_data* d);
int mutate_prob_vector(struct pbil_data* d);

/* Misc  */
unsigned int BinaryToGray(unsigned int num);
unsigned int GrayToBinary32(unsigned int num);

/* objective function */
int eval(struct pbil_data* d,char** infile, int num_infiles,struct thr_pool* pool);
void* run_kalign_thread(void *threadarg);
int run_kalign_and_score(char* infile,char* p_file,double* SP, double* TC);


/* seeding */
int set_pbil_based_on_pop(struct pbil_data* d);
int init_pop_from_seed(struct pbil_data*d, char* infile);
int read_aln_param_from_file(struct pbil_data*d, char* infile);

/* for testing */
int random_score(struct pbil_data* d);

int sort_pop_by_score(const void *a, const void *b);

int main(int argc, char *argv[])
{
        int i;
        //double SP,TC;

        struct pbil_data* d = NULL;


        int num_infiles = 0;
        char** infile = NULL;
        char* outfile = NULL;
        char* seedfile = NULL;
        int c = 0;


        //int help = 0;
        while (1){
                static struct option long_options[] ={
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

                case 's':
                        seedfile = optarg;
                        break;
                case 'o':
                        outfile = optarg;
                        break;
                case 'h':
                        //help = 1;
                        break;
                default:
                        ERROR_MSG("not recognized");
                        break;
                }
        }


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
        int num_threads = 5;
        struct thr_pool* pool;
        RUNP(pool = thr_pool_create(num_threads, num_threads, 0, 0));
        RUNP(d = init_pbil_data(234, 16, 1000, 10, 1));

        if(seedfile){

                RUN(init_pop_from_seed(d, seedfile));

        }
        //d->lambda = 5;


        for(i = 0; i < d->num_gen;i++){
                RUN(sample_pop(d));

                RUN(write_kalign_parameter_files(d));
                RUN(eval(d, infile, num_infiles, pool));
                //random_score(d);

                RUN(update_pbil(d));
                RUN(mutate_prob_vector(d));
        }


        //RUN(run_kalign_and_score("~/data/bb3_release/RV30/BB30001.xml","test_kalign_param.txt", &SP, &TC));
        //fprintf(stdout,"%f %f\n",SP,TC);
        free_pbil_data(d);

        return EXIT_SUCCESS;
ERROR:
        return EXIT_FAILURE;
}


int eval(struct pbil_data* d,char** infile, int num_infiles,struct thr_pool* pool)
{
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

        MMALLOC(td, sizeof(struct thread_data*) * num_jobs);
        for(i = 0; i < num_jobs;i++){
                td[i] = NULL;
                MMALLOC(td[i], sizeof(struct thread_data));
                td[i]->SP = 0.0;
                td[i]->TC = 0.0;
                td[i]->in = NULL;
                td[i]->param = NULL;
        }
        c = 0;
        for(i = 0; i < num_infiles;i++){

                for(j = 0; j < d->mu;j++){
                        td[c]->in = infile[i];
                        td[c]->param = param_file_name_buffer[j];
                        td[c]->param_index = j;
                        c++;
                }
        }

        for(i = 0; i < num_jobs;i++){
                if((status = thr_pool_queue(pool, run_kalign_thread, td[i])) == -1) ERROR_MSG("Adding job to queue failed.");
        }
        thr_pool_wait(pool);
        for(i = 0; i < num_jobs;i++){
                //fprintf(stdout,"%s %s: %f %f\n",td[i]->in,td[i]->param,td[i]->SP,td[i]->TC);
                d->population[td[i]->param_index]->score += td[i]->SP;
        }
        for(i = 0; i < d->mu;i++){

                d->population[i]->score /= (double) num_infiles;
                //fprintf(stdout,"SCORE: %d %f\n",i,d->population[i]->score);
        }


        for(i = 0; i < num_jobs;i++){
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


        struct thread_data *data;
        data = (struct thread_data *) threadarg;
        ASSERT(data != NULL, "No data");
        //fprintf(stdout,"%s %s\n", data->in,data->param);
        RUN(run_kalign_and_score(data->in, data->param, &data->SP, &data->TC));
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

        fprintf(stdout,"average:%f\n",sum / (double) d->lambda);
        if(d->population[0]->score > d->best->score){

                for(i = 0; i < d->num_param;i++){
                        d->best->param[i] = d->population[0]->param[i];
                }
                d->best->score =d->population[0]->score;

                fprintf(stdout,"New best: %f\t",d->best->score);
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

int random_score(struct pbil_data* d)
{
        int i;
        double r;
        for(i = 0; i < d->mu;i++){
                RUN(drand48_r(&d->randBuffer, &r));
                d->population[i]->score = r;

        }

        return OK;
ERROR:
        return FAIL;
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
                        RUN(drand48_r(&d->randBuffer, &r));
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
ERROR:
        return FAIL;
}




int mutate_prob_vector(struct pbil_data* d)
{
        int i;
        double r;
        for(i = 0; i < d->n_bits;i++){
                drand48_r(&d->randBuffer, &r);
                if(r < d->bit_mutation_prob){
                        drand48_r(&d->randBuffer, &r);

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
        srand48_r(time(NULL), &d->randBuffer);
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
                MFREE(d->bit_prob);
                MFREE(d);
        }
}

int run_kalign_and_score(char* infile,char* p_file, double* SP, double* TC)
{

        FILE *pf;
        char cmd[BUFFER_LEN];
        char ret[BUFFER_LEN];

        snprintf(cmd, BUFFER_LEN, "run_kalign_bali_score.sh -a %s -i %s ",p_file, infile);
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
        int c;
        double r;
        double tmp;
        RUN(read_aln_param_from_file(d,infile));

        for(i = 0; i < d->num_param;i++){
                d->min[i] = d->population[0]->param[i] - 1.0;
                d->max[i] = d->population[0]->param[i] + 1.0;
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
                RUN(drand48_r(&d->randBuffer, &r));
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

int read_aln_param_from_file(struct pbil_data*d, char* infile)
{
        char line_buffer[BUFFER_LEN];
        FILE* f_ptr = NULL;


        int m_pos;
        double val = 0;



        if(!my_file_exists(infile)){
                ERROR_MSG("File: %s does not exist", infile);
        }
        RUNP( f_ptr = fopen(infile, "r"));

        m_pos =0;
        while (fgets(line_buffer , BUFFER_LEN, f_ptr)){
                sscanf(line_buffer, "%lf", &val);
                d->population[0]->param[m_pos] = val;
                m_pos++;
                //fprintf(stdout,"%s", ret);
        }


        fclose(f_ptr);

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
