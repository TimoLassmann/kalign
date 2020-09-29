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

#ifndef BISECTINGKMEANS_H
#define BISECTINGKMEANS_H

#include "alignment_parameters.h"

struct aln_tasks;

int build_tree_kmeans(struct msa* msa, struct aln_param* ap,struct  aln_tasks** task_list);

//int build_tree_kmeans(struct msa* msa, struct aln_param* ap);

//extern int build_tree_kmeans(struct msa* msa, struct aln_param* ap);
#endif
