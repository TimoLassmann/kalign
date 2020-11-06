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

#ifndef PARAMETERS_H
#define PARAMETERS_H


#define KALIGNDIST_ALN 0
#define KALIGNDIST_BPM 1
#define KALIGNDIST_WU 2

struct parameters{
        char **infile;
        char *input;
        char *outfile;
        char* format;
        char* aln_param_file;
        float gpo;
        float gpe;
        float tgpe;
        float matadd;
        int chaos;
        int out_format;
        int param_set;
        int dist_method;
        int num_infiles;
        int reformat;
        int rename;             /* rename sequences - to make bali_score swallow the alignments */
        int dump_internal;
        int nthreads;
        int clean;
        int unalign;
        int help_flag;
};

extern struct parameters* init_param(void);
extern void free_parameters(struct parameters* param);
#endif
