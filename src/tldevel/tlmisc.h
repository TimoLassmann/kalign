#ifndef TLMISC_H
#define TLMISC_H



kalign_extern char* basename(const char* name);
kalign_extern int my_str_cpy(char* target, char* source, int t_size,int s_size);
kalign_extern int my_str_append(char* target, char* source, int t_size,int s_size);
kalign_extern int my_file_exists(const char* name);
kalign_extern int make_cmd_line(char** command, const int argc,char* const argv[]);
kalign_extern int tlfilename(char* path, char** out);
kalign_extern int tldirname(char* path, char** out);


#endif
