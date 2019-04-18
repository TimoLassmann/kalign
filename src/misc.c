

#include "misc.h"

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
