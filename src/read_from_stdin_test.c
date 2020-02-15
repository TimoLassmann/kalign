#include <unistd.h>
#include <stdlib.h>
#include <stdio.h>

#include "tldevel.h"

int main(void)

{
        char* line_buffer = NULL;
        size_t line_len = 0;
        ssize_t nread;
        if (isatty(fileno(stdin))){

                puts("stdin is connected to a terminal");
        }else{

                while ((nread = getline(&line_buffer, &line_len, stdin)) != -1){
                        fprintf(stdout,"%s",line_buffer);
                }
                puts("stdin is NOT connected to a terminal");
                MFREE(line_buffer);
        }


        return 0;
}
