#ifndef TLMISC_H
#define TLMISC_H


#ifdef TLMISC_IMPORT
#define EXTERN
#else
#ifdef __cplusplus
#define EXTERN extern "C"
#else
#define EXTERN extern
#endif
#endif

EXTERN char* basename(const char* name);
EXTERN int my_str_cpy(char* target, char* source, int t_size,int s_size);
EXTERN int my_str_append(char* target, char* source, int t_size,int s_size);
EXTERN int my_file_exists(const char* name);
EXTERN int make_cmd_line(char** command, const int argc,char* const argv[]);
EXTERN int tlfilename(char* path, char** out);
EXTERN int tldirname(char* path, char** out);

#undef TLMISC_IMPORT
#undef EXTERN

#endif
