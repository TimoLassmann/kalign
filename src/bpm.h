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

#ifndef BPM_H
#define BPM_H

#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

#include "tldevel.h"

/* Must be called before bpm_256!!!!  */
extern void set_broadcast_mask(void);

extern uint8_t bpm_256(const uint8_t* t,const uint8_t* p,int n,int m);
extern uint8_t bpm(const uint8_t* t,const uint8_t* p,int n,int m);



#endif
