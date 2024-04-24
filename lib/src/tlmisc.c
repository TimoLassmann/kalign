

#include <sys/stat.h>

#include <stdio.h>
#include <string.h>

#include "strnlen_compat.h"

#define TLMISC_IMPORT
#include "tlmisc.h"

#include "tldevel.h"

#define MAX_CMD_LEN 16384




int my_str_cpy(char* target, char* source, int t_size,int s_size)
{
        if(s_size > t_size){
                ERROR_MSG("Target size is < than source size");
        }
        int s_len = strnlen(source, s_size);
        int t_len = strnlen(target, t_size);

        if(s_len > t_size){
                ERROR_MSG("Target len is < than source size");
        }

        for(int i = 0; i < t_len;i++){
                target[i] = source[i];
        }
        target[t_len] = 0;

        return OK;
ERROR:
        return FAIL;
}

int my_str_append(char* target, char* source, int t_size,int s_size)
{
        if(s_size > t_size){
                ERROR_MSG("Target size is < than source size");
        }
        int s_len = strnlen(source, s_size);
        int t_len = strnlen(target, t_size);

        if(s_len > t_size){
                ERROR_MSG("Target len is < than source size");
        }

        if(t_size - t_len < s_len){
                ERROR_MSG("Target has insufficient space.");
        }
        int c = t_len;
        for(int i = 0; i < t_len;i++){
                target[c] = source[i];
                c++;
        }
        target[c] = 0;

        return OK;
ERROR:
        return FAIL;
}


char* basename(const char* name)
{
        int i= 0;
        int c = 0;

        while(1){
                if(name[i] == '/'){
                        c = i+1;
                }
                if(!name[i]){
                        break;
                }
                i++;
        }
        return (char*)(name +c);
}


int my_file_exists(const char* name)
{
        struct stat buf;
        int ret,local_ret;
        ret = 0;
        local_ret= stat ( name, &buf );
        /* File found */
        if ( local_ret == 0 )
        {
                ret++;
        }
        return ret;
}

/* I don't like that both libgen and string have functions to work with
   directory / filenames. The functions below copy the input path and alloc
   a new character array to store the output. Needs to be MFREE'd...
 */
int tlfilename(char* path, char** out)
{
        char* tmp = NULL;
        int len;
        int i;
        int c;
        len = 0;

        len = (int) strlen(path);
        MMALLOC(tmp, sizeof(char) * (len+1));
        c = 0;
        for(i = 0;i < len;i++){
                tmp[c] = path[i];
                c++;
                if(path[i] == '/'){
                        c = 0;
                }

        }
        tmp[c] = 0;
        if(c == 0){
                ERROR_MSG("No filename found in: %s", path);
        }
        *out = tmp;
        return OK;
ERROR:
        if(tmp){
                MFREE(tmp);
        }
        return FAIL;
}

int tldirname(char* path, char** out)
{
        char* tmp = NULL;
        int len;
        int i;
        int c;
        int e;
        len = 0;

        len = (int) strlen(path);
        MMALLOC(tmp, sizeof(char) * (len+1));
        c = 0;
        e = 0;
        for(i = 0;i < len;i++){

                tmp[c] = path[i];
                if(path[i] == '/'){
                        e = c;
                }
                c++;
        }
        tmp[e] = 0;
        if(e == 0){
                ERROR_MSG("No dirname found in: %s", path);
        }
        *out = tmp;
        return OK;
ERROR:
        if(tmp){
                MFREE(tmp);
        }
        return FAIL;
}

int make_cmd_line(char** command, const int argc,char* const argv[])
{
        char* cmd = NULL;
        int i,j,c,g;
        int alloc_len = 16;


        int old_len;
        RUN(galloc(&cmd,alloc_len));
        for(i =0 ; i < alloc_len;i++){
                cmd[i] = 0;
        }
        c = 0;
        for(i =0 ; i < argc;i++){
                //fprintf(stdout,"%s\n", argv[i]);
                for(j = 0; j < (int) strlen(argv[i]);j++){
                        if(c == 16384-1){
                                break;
                        }
                        cmd[c] = argv[i][j];
                        c++;
                        if(c == alloc_len){
                                old_len = alloc_len;
                                alloc_len = alloc_len + alloc_len /2;
                                RUN(galloc(&cmd,alloc_len));
                                for(g = old_len; g < alloc_len;g++){
                                        cmd[g] = 0;
                                }
                        }

                }
                if(c >= MAX_CMD_LEN){
                        ERROR_MSG("Command line too long! Allocated: %d", alloc_len);
                }
                if(c == 16384-1){
                        break;
                }
                cmd[c] = ' ';
                c++;
                if(c == alloc_len){
                        old_len = alloc_len;
                        alloc_len = alloc_len + alloc_len /2;
                        RUN(galloc(&cmd,alloc_len));
                        for(g = old_len; g < alloc_len;g++){
                                cmd[g] = 0;
                        }
                }


        }
        cmd[c] = 0;
        *command = cmd;
        return OK;
ERROR:
        if(cmd){
                gfree(cmd);
        }
        return FAIL;
}
