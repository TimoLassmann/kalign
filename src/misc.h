#ifndef MISC_H
#define MISC_H

#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

#include <stdlib.h>
#include <stdio.h>


#include "tldevel.h"



extern int byg_detect(int* text,int n);
extern int byg_start(char* pattern,char*text);
extern int byg_end(char* pattern,char*text);
extern int byg_count(char* pattern,char*text);

#endif
