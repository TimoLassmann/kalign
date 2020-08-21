

#include <stdlib.h>
#include <stdarg.h>
#include <stdio.h>




#define TLDEVEL_IMPORT
#include "tldevel.h"


#define TYPE_MARGIN 8

static void verror(FILE* f_ptr, const char *location, const char *format,  va_list argp);
static void vwarning(FILE* f_ptr,const char *location, const char *format,  va_list argp);
static void vlog(FILE* f_ptr,const char *format,  va_list argp);
static int get_time(char* time_ptr, int size);


typedef struct {
        int dim1;
        int dim2;
} mem_i;


const char* tldevel_version(void)
{
        return TLDEVEL_VERSION;
}


int galloc_unknown_type_error (void* p, ...)
{
        error(AT, "galloc was called with pointer of unknown type");
        return FAIL;
}

int galloc_too_few_arg_error (void* p)
{
        error(AT,"galloc was called with only one argument");
        return FAIL;
}


/* g memory stuff */

int get_dim1(void* ptr,int* d)
{
        if(ptr){
                *d = ((mem_i*)((void*) ((char*)ptr - sizeof(mem_i))))->dim1;
                return OK;
        }
        return FAIL;
}

int get_dim2(void* ptr,int* d)
{
        if(ptr){
                *d = ((mem_i*)((void*) ((char*)ptr - sizeof(mem_i))))->dim2;
                return OK;
        }
        return FAIL;
}


#define ALLOC_1D_ARRAY(type)                                            \
        int alloc_1D_array_size_ ##type (type **array, int dim1) {      \
                mem_i* h = NULL;                                        \
                void* tmp = NULL;                                       \
                ASSERT(dim1 >= 1,"DIM1 is too small: %d",dim1);         \
                if(*array == NULL){                                      \
                        MMALLOC(tmp,(dim1  * sizeof **array + sizeof(mem_i))); \
                }else{                                                  \
                        tmp = *array;                                   \
                        tmp = (void*) ((char*)tmp - sizeof(mem_i));     \
                        h = (mem_i*)(tmp);                              \
                        if(h->dim1 < dim1){                             \
                                MREALLOC(tmp,(dim1  * sizeof **array + sizeof(mem_i))); \
                        }else{                                          \
                                *array = (type*) ((char*)tmp + sizeof(mem_i)); \
                                return OK;                              \
                        }                                               \
                }                                                       \
                h = (mem_i*)(tmp);                                      \
                h->dim1  = dim1;                                        \
                h->dim2  = 0;                                           \
                *array= (type*)  ((char*)tmp + sizeof(mem_i));          \
                return OK;                                              \
        ERROR:                                                          \
                gfree(*array);                                         \
                return FAIL;                                            \
        }


ALLOC_1D_ARRAY(char)
ALLOC_1D_ARRAY(int8_t)
ALLOC_1D_ARRAY(uint8_t)
ALLOC_1D_ARRAY(int16_t)
ALLOC_1D_ARRAY(uint16_t)
ALLOC_1D_ARRAY(int32_t)
ALLOC_1D_ARRAY(uint32_t)
ALLOC_1D_ARRAY(int64_t)
ALLOC_1D_ARRAY(uint64_t)
ALLOC_1D_ARRAY(float)
ALLOC_1D_ARRAY(double)

#undef ALLOC_1D_ARRAY

#define ALLOC_2D_ARRAY(type)                                            \
        int alloc_2D_array_size_ ##type (type ***array, int dim1,int dim2) { \
                int i,j,c;                                              \
                mem_i* h = NULL;                                        \
                type** ptr_t = NULL;                                    \
                type* ptr_tt = NULL;                                    \
                void* tmp = NULL;                                       \
                int max1, max2;                                         \
                int o1, o2;                                             \
                ASSERT(dim1 >= 1,"DIM1 is too small: %d",dim1);         \
                ASSERT(dim2 >= 1,"DIM1 is too small: %d",dim2);         \
                if(*array == NULL){                                     \
                        MMALLOC(tmp,(dim1  * sizeof **array+ sizeof(mem_i))); \
                        MMALLOC(ptr_tt,((dim1 * dim2)  * sizeof ***array)); \
                        h = (mem_i*)tmp;                                \
                        h->dim1  = dim1;                                \
                        h->dim2  = dim2;                                \
                        max1 = dim1;                                    \
                        max2 = dim2;                                    \
                        ptr_t =(type**) ((char*)tmp + sizeof(mem_i));   \
                        for(i = 0;i< dim1;i++){                         \
                                ptr_t[i] = ptr_tt + i * dim2;           \
                        }                                               \
                        *array = ptr_t;                                 \
                }else{                                                  \
                        ptr_tt = *array[0];                             \
                        tmp = (void*)( (char*)*array -sizeof(mem_i));   \
                        h = (mem_i*)tmp;                                \
                        o1 = h->dim1;                                   \
                        o2 = h->dim2;                                   \
                        max1 = MACRO_MAX(dim1,o1);                      \
                        max2 = MACRO_MAX(dim2,o2);                      \
                        if(dim1 > o1){                                  \
                                MREALLOC(tmp,(dim1  * sizeof **array+ sizeof(mem_i))); \
                                MREALLOC(ptr_tt,((dim1* max2)  * sizeof ***array )); \
                        }else if(dim2 > o2){                            \
                                MREALLOC(ptr_tt,((max1 * dim2) * sizeof ***array )); \
                        }else{                                          \
                                return OK;                              \
                        }                                               \
                        if(dim2 > o2){                                  \
                                for(i = o1-1; i >= 0;i-- ){             \
                                        c =  i* max2;                   \
                                        for(j = o2-1;j >=0;j--){        \
                                                *(ptr_tt + c + j) =*(ptr_tt + i*o2 + j); \
                                        }                               \
                                }                                       \
                        }                                               \
                        h = (mem_i*)tmp;                                \
                        h->dim1 = max1;                                 \
                        h->dim2 = max2;                                 \
                        ptr_t = (type**) ((char*)tmp + sizeof(mem_i));  \
                        for(i = 0; i < max1;i++){                       \
                                ptr_t[i] = ptr_tt + i * max2;           \
                        }                                               \
                        *array = ptr_t;                                 \
                }                                                       \
                return OK;                                              \
        ERROR:                                                          \
                gfree(*array);                                          \
                return FAIL;                                            \
        }

ALLOC_2D_ARRAY(char)
ALLOC_2D_ARRAY(int8_t)
ALLOC_2D_ARRAY(uint8_t)
ALLOC_2D_ARRAY(int16_t)
ALLOC_2D_ARRAY(uint16_t)
ALLOC_2D_ARRAY(int32_t)
ALLOC_2D_ARRAY(uint32_t)
ALLOC_2D_ARRAY(int64_t)
ALLOC_2D_ARRAY(uint64_t)
ALLOC_2D_ARRAY(float)
ALLOC_2D_ARRAY(double)


#define FREE_VOID(type)                         \
        void gfree_void_ ##type(type *a){\
                error(AT, "free was called on wrong type (%p)",(void*)a); \
        }

FREE_VOID(char)
FREE_VOID(int8_t)
FREE_VOID(uint8_t)
FREE_VOID(int16_t)
FREE_VOID(uint16_t)
FREE_VOID(int32_t)
FREE_VOID(uint32_t)
FREE_VOID(int64_t)
FREE_VOID(uint64_t)
FREE_VOID(float)
FREE_VOID(double)

#undef FREE_VOID


#define FREE_1D_ARRAY(type)                               \
        void free_1d_array_ ##type(type **array){                        \
                if(*array){                                              \
                        void* ptr = (void*)((char*)*array - sizeof(mem_i)); \
                        MFREE(ptr);                                     \
                        *array = NULL;                                  \
                }                                                       \
        }


FREE_1D_ARRAY(char)
FREE_1D_ARRAY(int8_t)
FREE_1D_ARRAY(uint8_t)
FREE_1D_ARRAY(int16_t)
FREE_1D_ARRAY(uint16_t)
FREE_1D_ARRAY(int32_t)
FREE_1D_ARRAY(uint32_t)
FREE_1D_ARRAY(int64_t)
FREE_1D_ARRAY(uint64_t)
FREE_1D_ARRAY(float)
FREE_1D_ARRAY(double)

#undef FREE_1D_ARRAY


#define FREE_2D_ARRAY(type)                               \
        void free_2d_array_ ##type(type ***array){         \
        if(*array){                                                \
                if(*array[0]){                                     \
                        MFREE(*array[0]);                          \
                }                                                 \
                void* ptr = (void*)((char*)*array- sizeof(mem_i)); \
                MFREE(ptr);                                       \
                *array = NULL;                                    \
        }                                                         \
        }

FREE_2D_ARRAY(char)
FREE_2D_ARRAY(int8_t)
FREE_2D_ARRAY(uint8_t)
FREE_2D_ARRAY(int16_t)
FREE_2D_ARRAY(uint16_t)
FREE_2D_ARRAY(int32_t)
FREE_2D_ARRAY(uint32_t)
FREE_2D_ARRAY(int64_t)
FREE_2D_ARRAY(uint64_t)
FREE_2D_ARRAY(float)
FREE_2D_ARRAY(double)

#undef FREE_2D_ARRAY

void error(const char *location, const char *format, ...)
{
        va_list argp;
        va_start(argp, format);
        verror(stderr,location,format,argp);
        va_end(argp);
}

void warning(const char *location, const char *format, ...)
{
        va_list argp;
        va_start(argp, format);
        vwarning(stdout,location, format, argp);
        va_end(argp);
}

void log_message( const char *format, ...)
{
        va_list argp;
        va_start(argp, format);
        vlog(stdout,format, argp);
        va_end(argp);
}

int get_time(char* time_ptr, int size)
{
        struct tm *ptr;
        time_t current = time(NULL);
        ptr = localtime(&current);
        if(!strftime(time_ptr, size, "[%F %H:%M:%S] ", ptr))ERROR_MSG("write failed");
        return OK;
ERROR:
        return FAIL;
}

void verror(FILE* f_ptr, const char *location, const char *format,  va_list argp)
{
        char time_string[200];
        if(get_time(time_string, 200) != OK){
                fprintf(stderr,"notime");
        }
        fprintf(f_ptr,"%*s: ",MESSAGE_MARGIN,time_string);
        fprintf(f_ptr,"%*s: ",TYPE_MARGIN,"ERROR ");
        vfprintf(f_ptr, format, argp);
        fprintf(f_ptr," (%s)\n",location);
        fflush(f_ptr);
}

void vwarning(FILE* f_ptr,const char *location, const char *format,  va_list argp)
{
        char time_string[200];

        if(get_time(time_string, 200) != OK){
                fprintf(stderr,"notime");
        }
        fprintf(f_ptr,"%*s: ",MESSAGE_MARGIN,time_string);
        fprintf(f_ptr,"%*s: ",TYPE_MARGIN,"WARNING ");
        vfprintf(f_ptr, format, argp);
        fprintf(f_ptr," (%s)\n",location);
        fflush(f_ptr);
}

void vlog(FILE* f_ptr,const char *format,  va_list argp)
{

        char time_string[200];

        if(get_time(time_string, 200) != OK){
                fprintf(stderr,"notime");
        }
        fprintf(f_ptr,"%*s: ",MESSAGE_MARGIN,time_string);
        fprintf(f_ptr,"%*s: ",TYPE_MARGIN,"LOG ");
        vfprintf(f_ptr, format, argp);
        fprintf(f_ptr,"\n");
        fflush(f_ptr);
}


int nearly_equal_float(float a, float b)
{
        float absa = fabsf(a);
        float absb = fabsf(b);
        float d = fabsf(a-b);
        if(a == b){
                return 1;
        }else if (a == 0.0f || b == 0 || (absa + absb < FLT_MIN)){
                return d < (FLT_EPSILON * FLT_MIN);
        }else{
                return d / MACRO_MIN((absa+absb), FLT_MIN) < FLT_EPSILON;
        }
}

int nearly_equal_double(double a, double b)
{
        double absa = fabs(a);
        double absb = fabs(b);
        double d = fabs(a-b);
        if(a == b){
                return 1;
        }else if (a == 0.0f || b == 0 || (absa + absb < DBL_MIN)){
                return d < (DBL_EPSILON * DBL_MIN);
        }else{
                return d / MACRO_MIN((absa+absb), DBL_MIN) < DBL_EPSILON;
        }
}
