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


#include "global.h"
#include "parameters.h"





struct parameters*init_param(void)
{
        struct parameters* param = NULL;
        MMALLOC(param, sizeof(struct parameters));
        param->dist_method = KALIGNDIST_BPM;
        param->aln_param_file = NULL;
        param->param_set = -1;
        param->infile = NULL;
        param->num_infiles = 0;
        param->input = NULL;
        param->outfile = NULL;
        param->format = NULL;
        param->reformat = 0;
        param->rename = 0;
        param->help_flag = 0;
        param->dump_internal = 0;

        param->gpo = FLT_MAX;
        param->gpe = FLT_MAX;
        param->tgpe = FLT_MAX;
        param->matadd = 0.0F;
        param->chaos = 0;
        param->nthreads = 4;
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
