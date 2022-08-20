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

#include "kmeans.h"



double** kmeans(double** data,int* cluster_assignment, int len_a,int len_b, int k)
{
        double** means = NULL;
        double** tmp = NULL;
        double** tmp2 = NULL;
        double** tmp_ptr = NULL;
        double best_solution;
        double score;
        double old_score;
        double min;
        double d;
        int min_index;
        int i,j;

        int num_attempts = 10000;
        int a_item;

        struct rng_state* rng;
        int* sel = NULL;
        double* n_item = NULL;

        ASSERT(k > 0,"K needs to be greater than one");
        ASSERT(k < len_a,"K larger than number of items");

        RUN(rng = init_rng(0));
        //srand48_r(time(NULL), &randBuffer);


        sel = galloc(sel,len_a);

        for(i = 0; i < len_a;i++){
                sel[i] = i;
        }


        n_item = galloc(n_item,k);

        for(i = 0; i < k;i++){
                n_item[i] = 0.0;
        }
        tmp = galloc(tmp,k,len_b,0.0f);
        tmp2= galloc(tmp2,k,len_b,0.0f);
        means= galloc(means,k,len_b,0.0f);

        best_solution = DBL_MAX;

        for(a_item = 0; a_item < num_attempts;a_item++){
                /* reset tmp */
                for(i = 0; i < k;i++){
                        for(j = 0; j < len_b;j++){
                                tmp[i][j] = 0.0;
                        }
                }
                /* shuffle to select first k means */
                shuffle_arr_r(sel, len_a,rng);


                /* initial selection  */

                for(i = 0; i < k;i++){
                        for(j = 0; j < len_b;j++){
                                tmp[i][j] = data[sel[i]][j];
                                tmp2[i][j] = 0.0;
                        }
                }

                score =1.0;
                old_score = 0.0;
                while(score != old_score){
                        for(i = 0; i < k;i++){
                                n_item[i] = 0.0;
                        }
                        for(i = 0; i < k;i++){
                                for(j = 0; j < len_b;j++){

                                        tmp2[i][j] = 0.0;
                                }
                        }

                        score = 0.0;
                        for(i= 0; i < len_a;i++){
                                min = DBL_MAX;
                                min_index = -1;
                                for(j = 0; j < k;j++){
                                        edist_serial_d(data[i], tmp[j], len_b, &d);
                                        if(d < min){
                                                min = d;
                                                min_index = j;
                                        }
                                }
                                n_item[min_index]++;
                                score += min;
                                for(j = 0; j < len_b;j++){
                                        tmp2[min_index][j] += data[i][j];
                                }
                                if(cluster_assignment){
                                        cluster_assignment[i] = min_index;
                                }
                        }
                        /* calculate new means  */
                        for(i = 0; i < k;i++){
                                for(j = 0; j < len_b;j++){

                                        tmp2[i][j] /= n_item[i];
                                }
                        }
                        /* switch */

                        tmp_ptr = tmp;
                        tmp= tmp2;
                        tmp2 = tmp_ptr;
                        if(old_score == score){
                                break;
                        }
                        old_score = score;
                }

                if(score < best_solution){
                        fprintf(stdout,"%f better than %f\n", score,best_solution);
                        best_solution = score;
                        for(i = 0; i < k;i++){
                                for(j = 0; j < len_b;j++){
                                        means[i][j] = tmp[i][j];
                                }
                        }
                }
        }
        for(i = 0; i < k;i++){
                n_item[i] = 0.0;
        }

        for(i= 0; i < len_a;i++){
                if(cluster_assignment){
                        n_item[cluster_assignment[i]]++;
                }
        }
        for(i = 0; i < k;i++){
                fprintf(stdout,"%d %f\n",i,n_item[i]);
        }

        gfree(tmp);
        gfree(tmp2);

        gfree(n_item);
        gfree(sel);

        return means;
ERROR:
        return NULL;
}
