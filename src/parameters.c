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

#include "tldevel.h"

#include <string.h>
#define PARAMETERS_IMPORT
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

        param->type = -1;

        param->gpo = -1.0;
        param->gpe = -1.0;
        param->tgpe = -1.0;
        param->matadd = 0.0F;
        param->chaos = 0;
        param->nthreads = 4;
        param->clean = 0;
        param->unalign = 0;
        param->quiet = 0;
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
