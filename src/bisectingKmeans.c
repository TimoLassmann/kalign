/*
  Kalign - a multiple sequence alignment program

  Copyright 2006, 2019 Timo Lassmann

  This file is part of kalign.

  Kalign is free software: you can redistribute it and/or modify
  it under the terms of the GNU General Public License as published by
  the Free Software Foundation, either version 3 of the License, or
  (at your option) any later version.

  This program is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
  GNU General Public License for more details.

  You should have received a copy of the GNU General Public License
  along with this program.  If not, see <https://www.gnu.org/licenses/>.

*/
#ifdef HAVE_OPENMP
#include <omp.h>
#endif

#ifdef HAVE_AVX2
#include <xmmintrin.h>
#endif

#include <mm_malloc.h>

#include "tlrng.h"
#include "msa.h"

#include "bisectingKmeans.h"

#include "global.h"
#include "msa.h"
#include "aln_task.h"
#include "sequence_distance.h"
#include "euclidean_dist.h"

/* #include "alignment.h" */
#include "pick_anchor.h"
#include "esl_stopwatch.h"

struct node{
        struct node* left;
        struct node* right;
        int id;
};

struct kmeans_result{
        int* sl;
        int* sr;
        int nl;
        int nr;
        float score;
};

static struct kmeans_result* alloc_kmeans_result(int num_samples);
static void free_kmeans_results(struct kmeans_result* k);

struct node* upgma(float **dm,int* samples, int numseq);
static struct node* alloc_node(void);

static int label_internal(struct node*n, int label);
static void create_tasks(struct node*n, struct aln_tasks* t);


static struct node* bisecting_kmeans_serial(struct msa* msa, struct node* n, float** dm,int* samples,int numseq, int num_anchors,int num_samples);
static struct node* bisecting_kmeans_parallel(struct msa* msa, struct node* n, float** dm,int* samples,int numseq, int num_anchors,int num_samples);

static int split(float** dm,int* samples, int num_anchors,int num_samples,int seed_pick,struct kmeans_result** ret);

int build_tree_kmeans(struct msa* msa, struct aln_param* ap,struct aln_tasks** tasks)
{
        struct aln_tasks* t = NULL;
        struct node* root = NULL;
        float** dm = NULL;
        int* samples = NULL;
        int* anchors = NULL;
        int num_anchors;
        int numseq;

        int i;

        ASSERT(msa != NULL, "No alignment.");
        //ASSERT(param != NULL, "No input parameters.");
        ASSERT(ap != NULL, "No alignment parameters.");


        t = *tasks;
        if(!t){
                RUN(alloc_tasks(&t, msa->numseq));
        }
        numseq = msa->numseq;

        DECLARE_TIMER(timer);
        /* pick anchors . */
        LOG_MSG("Calculating pairwise distances");
        START_TIMER(timer);
        RUNP(anchors = pick_anchor(msa, &num_anchors));

        RUNP(dm = d_estimation(msa, anchors, num_anchors,0));//les,int pair)

        STOP_TIMER(timer);
        GET_TIMING(timer);
        //LOG_MSG("Done in %f sec.", GET_TIMING(timer));

        MFREE(anchors);

        MMALLOC(samples, sizeof(int)* numseq);
        for(i = 0; i < numseq;i++){
                samples[i] = i;
        }
        LOG_MSG("%d anchors ", num_anchors);
        //RUNP(root = alloc_node());

        START_TIMER(timer);
        LOG_MSG("Building guide tree.");

        if(ap->nthreads == 1){
                root = bisecting_kmeans_serial(msa,root, dm, samples, numseq, num_anchors, numseq);
        }else{
                root = bisecting_kmeans_parallel(msa,root, dm, samples, numseq, num_anchors, numseq);
        }

        STOP_TIMER(timer);
        GET_TIMING(timer);

        label_internal(root, numseq);

        create_tasks(root, t);


        /*exit(0);
        ap->tree[0] = 1;
        ap->tree = readbitree(root, ap->tree);
        for (i = 0; i < (numseq*3);i++){
                tree[i] = tree[i+1];
                }*/
        MFREE(root);
        for(i =0 ; i < msa->numseq;i++){
                _mm_free(dm[i]);
        }
        MFREE(dm);
        DESTROY_TIMER(timer);
        return OK;
ERROR:
        return FAIL;
}





struct node* bisecting_kmeans_parallel(struct msa* msa, struct node* n, float** dm,int* samples,int numseq, int num_anchors,int num_samples)
{
        struct kmeans_result* res_tmp = NULL;
        struct kmeans_result* best = NULL;
        struct kmeans_result** res_ptr = NULL;
        int i,j;
        int tries = 40;
        /* int t_iter; */
        /* int r; */
        int* sl = NULL;
        int* sr = NULL;
        int num_l,num_r;
        if(num_samples < 100){
                float** dm = NULL;
                RUNP(dm = d_estimation(msa, samples, num_samples,1));// anchors, num_anchors,1));
                n = upgma(dm,samples, num_samples);

                gfree(dm);
                MFREE(samples);
                return n;
        }else if(num_samples < 1000){
                n = bisecting_kmeans_serial(msa, n, dm, samples, numseq, num_anchors, num_samples);
                return n;
        }


        best = NULL;
        res_tmp = NULL;

        MMALLOC(res_ptr, sizeof(struct kmeans_result*) * 4);
        for(i = 0; i < 4;i++){
                res_ptr[i] = NULL;
        }
        tries = MACRO_MIN(tries, num_samples);
        int step = num_samples / tries;
        int change = 0;
        for(i = 0;i < tries;i += 4){
                change = 0;
#ifdef HAVE_OPENMP
#pragma omp parallel for num_threads(4)
#endif
                for(j = 0; j < 4;j++){
                        split(dm,samples,num_anchors, num_samples, (i+ j)*step, &res_ptr[j]);
                }


                for(j = 0; j < 4;j++){
                        if(!best){
                                change++;
                                best = res_ptr[j];
                                res_ptr[j] = NULL;
                        }else{
                                if(best->score > res_ptr[j]->score){
                                        res_tmp = best;
                                        best = res_ptr[j];
                                        res_ptr[j] = res_tmp;
                                        /* LOG_MSG("Better!!! %f %f", res_tmp->score,best->score); */

                                        change++;
                                }



                        }
                }
                if(!change){
                        break;
                }

        }
        sl = best->sl;
        sr = best->sr;

        num_l = best->nl;

        num_r = best->nr;

        for(i = 0; i < 4;i++){
                free_kmeans_results(res_ptr[i]);
        }
        MFREE(res_ptr);
        MFREE(best);

        MFREE(samples);
        n = alloc_node();

#ifdef HAVE_OPENMP
#pragma omp parallel num_threads(2)
        {
#pragma omp single nowait
                {
                        //LOG_MSG("Done");

#pragma omp task shared(msa,n,dm,numseq,num_anchors)
#endif
                        {
                                n->left = bisecting_kmeans_parallel(msa,n->left, dm, sl, numseq, num_anchors,num_l);
                        }

#ifdef HAVE_OPENMP
#pragma omp task shared(msa,n,dm,numseq,num_anchors)
#endif
                        {
                                n->right = bisecting_kmeans_parallel(msa,n->right, dm, sr, numseq, num_anchors, num_r);
                        }

#ifdef HAVE_OPENMP
#pragma omp taskwait
                }
        }
#endif
        return n;
ERROR:
        return NULL;
}

struct node* bisecting_kmeans_serial(struct msa* msa, struct node* n, float** dm,int* samples,int numseq, int num_anchors,int num_samples)
{
        struct kmeans_result* res_tmp = NULL;
        struct kmeans_result* best = NULL;
        struct kmeans_result** res_ptr = NULL;
        int i,j;
        int tries = 40;
        /* int t_iter; */
        /* int r; */
        int* sl = NULL;
        int* sr = NULL;
        int num_l,num_r;
        if(num_samples < 100){
                float** dm = NULL;
                RUNP(dm = d_estimation(msa, samples, num_samples,1));// anchors, num_anchors,1));
                n = upgma(dm,samples, num_samples);

                gfree(dm);
                MFREE(samples);
                return n;
        }

        best = NULL;
        res_tmp = NULL;

        MMALLOC(res_ptr, sizeof(struct kmeans_result*) * 4);
        for(i = 0; i < 4;i++){
                res_ptr[i] = NULL;
        }
        tries = MACRO_MIN(tries, num_samples);
        int step = num_samples / tries;
        int change = 0;
        for(i = 0;i < tries;i += 4){
                change = 0;
                for(j = 0; j < 4;j++){
                        split(dm,samples,num_anchors, num_samples, (i+ j)*step, &res_ptr[j]);
                }


                for(j = 0; j < 4;j++){
                        if(!best){
                                change++;
                                best = res_ptr[j];
                                res_ptr[j] = NULL;
                        }else{
                                if(best->score > res_ptr[j]->score){
                                        res_tmp = best;
                                        best = res_ptr[j];
                                        res_ptr[j] = res_tmp;
                                        /* LOG_MSG("Better!!! %f %f", res_tmp->score,best->score); */

                                        change++;
                                }



                        }
                }
                if(!change){
                        break;
                }

        }
        sl = best->sl;
        sr = best->sr;

        num_l = best->nl;

        num_r = best->nr;

        for(i = 0; i < 4;i++){
                free_kmeans_results(res_ptr[i]);
        }
        MFREE(res_ptr);
        MFREE(best);

        MFREE(samples);
        n = alloc_node();

        n->left = bisecting_kmeans_serial(msa,n->left, dm, sl, numseq, num_anchors,num_l);
        n->right = bisecting_kmeans_serial(msa,n->right, dm, sr, numseq, num_anchors, num_r);
        return n;
ERROR:
        return NULL;
}

int split(float** dm,int* samples, int num_anchors,int num_samples,int seed_pick,struct kmeans_result** ret)
{
        struct kmeans_result* res = NULL;
        int* sl = NULL;
        int* sr = NULL;
        int num_l,num_r;
        float* w = NULL;
        float* wl = NULL;
        float* wr = NULL;
        float* cl = NULL;
        float* cr = NULL;
        float dl = 0.0F;
        float dr = 0.0F;
        float score;
        int num_var;
        int i;
        int s;
        int j;
        int stop = 0;

        num_var = num_anchors / 8;
        if( num_anchors%8){
                num_var++;
        }
        num_var = num_var << 3;


        wr = _mm_malloc(sizeof(float) * num_var,32);
        wl = _mm_malloc(sizeof(float) * num_var,32);
        cr = _mm_malloc(sizeof(float) * num_var,32);
        cl = _mm_malloc(sizeof(float) * num_var,32);
        w = _mm_malloc(sizeof(float) * num_var,32);
        if(*ret){
                res = *ret;

        }else{
                RUNP(res = alloc_kmeans_result(num_samples));
        }

        res->score = FLT_MAX;

        sl = res->sl;
        sr = res->sr;


        for(i = 0; i < num_var;i++){
                w[i] = 0.0F;
                wr[i] = 0.0F;
                wl[i] = 0.0F;
                cr[i] = 0.0F;
                cl[i] = 0.0F;
        }
        for(i = 0; i < num_samples;i++){
                s = samples[i];
                for(j = 0; j < num_anchors;j++){
                        w[j] += dm[s][j];
                }
        }

        for(j = 0; j < num_anchors;j++){
                w[j] /= (float)num_samples;
        }
        //r = tl_random_int(rng  , num_samples);
        //r = sel[t_iter];

        s = samples[seed_pick];
        /* LOG_MSG("Selected %d\n",s); */
        for(j = 0; j < num_anchors;j++){
                cl[j] = dm[s][j];
        }

        for(j = 0; j < num_anchors;j++){
                cr[j] = w[j] - (cl[j] - w[j]);
                //      fprintf(stdout,"%f %f  %f\n", cl[j],cr[j],w[j]);
        }

        _mm_free(w);

        /* check if cr == cl - we have identical sequences  */
        s = 0;
        for(j = 0; j < num_anchors;j++){
                if(fabsf(cl[j]-cr[j]) >  1.0E-6){
                        s = 1;
                        break;
                }
        }

        if(!s){
                score = 0.0F;
                num_l = 0;
                num_r = 0;
                sl[num_l] = samples[0];
                num_l++;

                for(i =1 ; i <num_samples;i++){
                        sr[num_r] = samples[i];
                        num_r++;
                }
        }else{
                w = NULL;
                while(1){
                        stop++;
                        if(stop == 10000){
                                ERROR_MSG("Failed.");
                        }
                        num_l = 0;
                        num_r = 0;

                        for(i = 0; i < num_anchors;i++){
                                wr[i] = 0.0F;
                                wl[i] = 0.0F;
                        }
                        score = 0.0f;
                        for(i = 0; i < num_samples;i++){
                                s = samples[i];
#ifdef HAVE_AVX2
                                edist_256(dm[s], cl, num_anchors, &dl);
                                edist_256(dm[s], cr, num_anchors, &dr);
#else
                                edist_serial(dm[s], cl, num_anchors, &dl);
                                edist_serial(dm[s], cr, num_anchors, &dr);
#endif
                                score += MACRO_MIN(dl,dr);

                                if(dr < dl){
                                        w = wr;
                                        sr[num_r] = s;
                                        num_r++;
                                }else if (dr > dl){
                                        w = wl;
                                        sl[num_l] = s;
                                        num_l++;
                                }else{
                                        if(i & 1){
                                                w = wr;
                                                sr[num_r] = s;
                                                num_r++;
                                        }else{
                                                w = wl;
                                                sl[num_l] = s;
                                                num_l++;
                                        }
                                }

                                for(j = 0; j < num_anchors;j++){
                                        w[j] += dm[s][j];
                                }

                        }

                        for(j = 0; j < num_anchors;j++){
                                wl[j] /= (float)num_l;
                                wr[j] /= (float)num_r;
                        }

                        s = 0;

                        for(j = 0; j < num_anchors;j++){
                                if(wl[j] != cl[j]){
                                        s = 1;
                                        break;
                                }
                                if(wr[j] != cr[j]){
                                        s = 1;
                                        break;
                                }
                        }

                        if(s){
                                w = cl;
                                cl = wl;
                                wl = w;

                                w = cr;
                                cr = wr;
                                wr = w;
                        }else{
                                break;
                        }
                }
        }
        _mm_free(wr);
        _mm_free(wl);
        _mm_free(cr);
        _mm_free(cl);

        res->nl =  num_l;
        res->nr =  num_r;
        res->score = score;

        *ret = res;
        return OK;
ERROR:
        return FAIL;
}


struct node* upgma(float **dm,int* samples, int numseq)
{
        struct node** tree = NULL;
        struct node* tmp = NULL;

        int i,j;
        int *as = NULL;

        float max;
        int node_a = 0;
        int node_b = 0;
        int cnode = numseq;
        int numprofiles;


        numprofiles = (numseq << 1) - 1;

        MMALLOC(as,sizeof(int)*numseq);
        for (i = numseq; i--;){
                as[i] = i+1;
        }


        MMALLOC(tree,sizeof(struct node*)*numseq);
        for (i = 0;i < numseq;i++){
                tree[i] = NULL;
                tree[i] = alloc_node();
                tree[i]->id = samples[i];
        }

        while (cnode != numprofiles){
                max = FLT_MAX;
                for (i = 0;i < numseq-1; i++){
                        if (as[i]){
                                for ( j = i + 1;j < numseq;j++){
                                        if (as[j]){
                                                if (dm[i][j] < max){
                                                        max = dm[i][j];
                                                        node_a = i;
                                                        node_b = j;
                                                }
                                        }
                                }
                        }
                }
                tmp = NULL;
                tmp = alloc_node();
                tmp->left = tree[node_a];
                tmp->right = tree[node_b];


                tree[node_a] = tmp;
                tree[node_b] = NULL;

                /*deactivate  sequences to be joined*/
                as[node_a] = cnode+1;
                as[node_b] = 0;


                cnode++;

                /*calculate new distances*/
                for (j = numseq;j--;){
                        if (j != node_b){
                                dm[node_a][j] = (dm[node_a][j] + dm[node_b][j])*0.5F + 0.001F;

                        }
                        //fprintf(stdout,"\n");
                }
                dm[node_a][node_a] = 0.0F;
                for (j = numseq;j--;){
                        dm[j][node_a] = dm[node_a][j];
                }
        }
        tmp = tree[node_a];
        MFREE(tree);
        MFREE(as);
        return tmp;
ERROR:
        return NULL;
}

struct node* alloc_node(void)
{
        struct node* n = NULL;
        MMALLOC(n, sizeof(struct node));
        n->left = NULL;
        n->right = NULL;
        n->id = -1;
        return n;
ERROR:
        return NULL;
}

int label_internal(struct node*n, int label)
{
        //n->d = d;
        if(n->left){
                label = label_internal(n->left, label);
        }
        if(n->right){
                label = label_internal(n->right, label);
        }
        if(n->id == -1){
                n->id = label;
                label++;
        }
        return label;

}

void create_tasks(struct node*n, struct aln_tasks* t)
{


        if(n->left && n->right){
                struct task* task;

                task = t->list[t->n_tasks];
                task->a = n->left->id;
                task->b = n->right->id;
                task->c = n->id;
                /* task->p = n->d; */
                /* task->n = n->n; */
                //fprintf(stdout,"Node %d  %d  depends on %d %d \n", n->id, n->d  , n->left->id, n->right->id);

                t->n_tasks++;
        }
        if(n->left){
                create_tasks(n->left,t);
        }
        if(n->right){
                create_tasks(n->right,t);
        }
        if(n->left){
                if(n->right){
                        MFREE(n->left);
                        MFREE(n->right);
                }
        }
}


struct kmeans_result* alloc_kmeans_result(int num_samples)
{
        struct kmeans_result* k = NULL;
        ASSERT(num_samples != 0, "No samples???");

        MMALLOC(k, sizeof(struct kmeans_result));

        k->nl = 0;
        k->nr = 0;
        k->sl = NULL;
        k->sr = NULL;
        MMALLOC(k->sl, sizeof(int) * num_samples);
        MMALLOC(k->sr, sizeof(int) * num_samples);
        k->score = FLT_MAX;
        return k;
ERROR:
        free_kmeans_results(k);
        return NULL;
}

void free_kmeans_results(struct kmeans_result* k)
{
        if(k){
                if(k->sl){
                        MFREE(k->sl);
                }
                if(k->sr){
                        MFREE(k->sr);
                }
                MFREE(k);
        }
}
