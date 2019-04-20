

#include "misc.h"

int byg_detect(int* text,int n)
{
        int Tc;
        int i  = 0;
        int s = 0;
        int T[256];
        for (i = 0;i < 256;i++){
                T[i] = 0;
        }
        int mb = 1;
        //char *unique_aa = "EFILPQXZ";//permissiv
        //ABCDEFGHIJKLMNOPQRSTUVWXYZ
        char *unique_aa = "BDEFHIJKLMNOPQRSVWYZ";//restrictive
        int aacode[26] = {0,1,2,3,4,5,6,7,8,-1,9,10,11,12,23,13,14,15,16,17,17,18,19,20,21,22};
        for (i= 0;i < 20;i++){
                T[(int)aacode[unique_aa[i]-65]] |= 1;
        }
        for (i = 0;i < n;i++){
                //	fprintf(stderr,"%d\n",text[i]);
                if(text[i] != -1){
                        s <<= 1;
                        s |= 1;
                        Tc = T[text[i]];
                        s &= Tc;
                        if(s & mb){
                                return 0;
                        }
                }
        }
        return 1;
}


int byg_count(char* pattern,char*text)
{
        int Tc;
        int count = 0;
        int i  = 0;
        int s = 0;
        int T[256];
        for (i = 0;i < 256;i++){
                T[i] = 0;
        }

        int m = strlen(pattern);
        int n = strlen (text);
        int mb = (1 << (m-1));

        for (i= 0;i < m;i++){
                T[(int)pattern[i]] |= (1 << i);
        }

        for (i = 0;i < n;i++){
                s <<= 1;
                s |= 1;
                Tc = T[(int)text[i]];
                s &= Tc;
                if(s & mb){
                        count++;
                }
        }
        return count;
}

int byg_end(char* pattern,char*text)
{
        int Tc;
        int i  = 0;
        int s = 0;
        int T[256];
        for (i = 0;i < 256;i++){
                T[i] = 0;
        }

        int m = strlen(pattern);
        int n = strlen (text);
        int mb = (1 << (m-1));

        for (i= 0;i < m;i++){
                T[(int)pattern[i]] |= (1 << i);
        }

        for (i = 0;i < n;i++){
                s <<= 1;
                s |= 1;
                if(!text[i]){
                        return -1;
                }
                Tc = T[(int)text[i]];
                s &= Tc;
                if(s & mb){
                        return i+1;
                }
        }
        return -1;
}


int byg_start(char* pattern,char*text)
{
        int Tc;
        int i  = 0;
        int s = 0;
        int T[256];
        for (i = 0;i < 256;i++){
                T[i] = 0;
        }

        int m = strlen(pattern);
        int n = strlen(text);
        int mb = (1 << (m-1));

        for (i= 0;i < m;i++){
                T[(int)pattern[i]] |= (1 << i);
        }

        for (i = 0;i < n;i++){
                s <<= 1;
                s |= 1;
                Tc = T[(int)text[i]];
                s &= Tc;
                if(s & mb){
                        return i-m+1;
                }
        }
        return -1;
}
