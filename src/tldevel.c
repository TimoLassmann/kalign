#ifdef HAVE_CONFIG_H
#include "config.h"
#endif
#include "tldevel.h"

static void verror(FILE* f_ptr, const char *location, const char *format,  va_list argp);
static void vwarning(FILE* f_ptr,const char *location, const char *format,  va_list argp);
static void vlog(FILE* f_ptr,const char *format,  va_list argp);
static void vmessage(FILE* f_ptr,const char *location, const char *format,  va_list argp);
static void print_program_description(char * const argv[],const char* description);
static void echo_build_config (void);


FREE_VOID(char)
FREE_VOID(int)
FREE_VOID(double)
FREE_1D_ARRAY(char)
FREE_1D_ARRAY(int)
FREE_1D_ARRAY(ulong)
FREE_1D_ARRAY(float)
FREE_1D_ARRAY(double)

FREE_2D_ARRAY(char)
FREE_2D_ARRAY(int)
FREE_2D_ARRAY(ulong)
FREE_2D_ARRAY(float)
FREE_2D_ARRAY(double)


ALLOC_1D_ARRAY(char)
ALLOC_1D_ARRAY(int)
ALLOC_1D_ARRAY(ulong)
ALLOC_1D_ARRAY(float)
ALLOC_1D_ARRAY(double)

ALLOC_2D_ARRAY(char)
ALLOC_2D_ARRAY(int)
ALLOC_2D_ARRAY(ulong)
ALLOC_2D_ARRAY(float)
ALLOC_2D_ARRAY(double)


int print_program_header(char* const argv[],const char* description)
{
        echo_build_config();
        print_program_description(argv,description);
        return OK;
}


void print_program_description(char * const argv[],const char* description)
{
        int i;
        int newline = 0;
        fprintf(stdout,"%-*s: %s\n" ,MESSAGE_MARGIN,"Running program", basename(argv[0]));
        if(description){
                newline = 1;
                fprintf(stdout,"%-*s: ",MESSAGE_MARGIN,"Description");
                for(i = 0;i < (int) strlen(description);i++){
                        if(description[i] == '\n'){
                                fprintf(stdout,"\n");
                                fprintf(stdout,"%-*s: ",MESSAGE_MARGIN,"");
                                newline = 1;

                        }else{
                                if(isspace(description[i]) &&newline){

                                }else{
                                        fprintf(stdout,"%c",description[i]);
                                        newline = 0;
                                }
                        }
                }
                fprintf(stdout,"\n");
        }
        fflush(stdout);
}


void echo_build_config (void)
{
        fprintf(stdout,"\n%s\n",build_config);
        fflush(stdout);
}

int log_command_line(const int argc,char* const argv[])
{
        char* buffer = NULL;

        RUNP(buffer = make_cmd_line(argc,argv));
        LOG_MSG("%s",buffer);
        MFREE(buffer);
        return OK;
ERROR:
        MFREE(buffer);
        return FAIL;
}

char* make_cmd_line(const int argc,char* const argv[])
{
        char* cmd = NULL;
        int i,j,c;

        MMALLOC(cmd, 16384);

        c = 0;
        for(i =0 ; i < argc;i++){
                for(j = 0; j < strlen(argv[i]);j++){
                        if(c == 16384-1){
                                break;
                        }
                        cmd[c] = argv[i][j];
                        c++;

                }
                if(c == 16384-1){
                        break;
                }
                cmd[c] = ' ';
                c++;


        }
        cmd[c] = 0;

        return cmd;
ERROR:
        return NULL;
}


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

void message(const char *location, const char *format, ...)
{
        va_list argp;
        va_start(argp, format);
        vmessage(stdout,location,format,argp);
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


void vmessage(FILE* f_ptr,const char *location, const char *format,  va_list argp)
{
        char time_string[200];

        if(get_time(time_string, 200) != OK){
                fprintf(stderr,"notime");
        }
        fprintf(f_ptr,"%*s: ",MESSAGE_MARGIN,time_string);
        fprintf(f_ptr,"%*s: ",TYPE_MARGIN,"MESSAGE ");
        vfprintf(f_ptr, format, argp);
        fprintf(f_ptr," (%s)\n",location);
        fflush(f_ptr);
}


int set_checkpoint_file(struct checkpoint* chk,char* function,char* location,char* cmd)
{
        char buffer[BUFFER_LEN];
        FILE* f_ptr = NULL;

        struct tm *ptr;

        char time_string[200];

        time_t current = time(NULL);
        ptr = localtime(&current);

        if(!strftime(time_string, 200, "%F %H:%M:%S", ptr)){
                ERROR_MSG("Write failed");
        }
        snprintf(buffer,BUFFER_LEN ,"%s/%s_%d.chk", chk->base_dir,chk->base_name,chk->test_num );
        RUNP(f_ptr = fopen(buffer , "w" ));
        fprintf(f_ptr,"%*s: %s\n",MESSAGE_MARGIN, "command", cmd);
        fprintf(f_ptr,"%*s: %d\n",MESSAGE_MARGIN, "checkpoint ID", chk->test_num);
        fprintf(f_ptr,"%*s: %s\n",MESSAGE_MARGIN, "function", function);
        fprintf(f_ptr,"%*s: %s\n",MESSAGE_MARGIN, "called in", location);
        fprintf(f_ptr,"%*s: %s\n",MESSAGE_MARGIN, "at time", time_string);

        fclose(f_ptr);

        return OK;
ERROR:
        if(f_ptr){
                fclose(f_ptr);
        }
        return FAIL;
}

int test_for_checkpoint_file(struct checkpoint* chk,char* function,char* location, char* cmd)
{
        FILE* f_ptr = NULL;
        char buffer[BUFFER_LEN];
        static int8_t found = 0;

        snprintf(buffer,BUFFER_LEN ,"%s/%s_%d.chk", chk->base_dir,chk->base_name,chk->test_num );
        if(my_file_exists(buffer) && !found){
                RUNP(f_ptr = fopen(buffer , "r" ));
                /* get first line and compare to  */
                buffer[0]= 0;
                if(fscanf(f_ptr,"%*s %99[^\n]s",buffer) != 1){
                        ERROR_MSG("fscanf failed.");
                }
                fclose(f_ptr);
                //fprintf(stdout,"%s\n%s\n",cmd,buffer);
                if(!strncmp(cmd,buffer,99)){
                        return 1;
                }

                LOG_MSG("   Re-running: %s (%s)",function,location);
                LOG_MSG("   arguments have changed from:");
                LOG_MSG("     %s",cmd);
                LOG_MSG("   to:");
                LOG_MSG("     %s",buffer);
                LOG_MSG("   will re-run everything from this point.");

                found = 1;
        }else{
                found = 1;
        }
        return 0;
ERROR:
        LOG_MSG("test_for_checkpoint file has failed.");
        return 0;
}

struct checkpoint* init_checkpoint(char* base_name,char* target_dir)
{
        struct checkpoint* chk = NULL;
        size_t i = 0;
        int j;
        MMALLOC(chk, sizeof(struct checkpoint));

        chk->test_num = 0;
        chk->base_dir = NULL;
        chk->base_name = NULL;

        i = strlen(target_dir);
        MMALLOC(chk->base_dir, sizeof(char) * (i+1));

        for(j = 0;j < i;j++){
                chk->base_dir[j] = target_dir[j];
        }
        chk->base_dir[i] = 0;

        i = strlen(base_name);
        MMALLOC(chk->base_name, sizeof(char) * (i+1));
        for(j = 0;j < i;j++){
                chk->base_name[j] = base_name[j];
        }
        chk->base_name[i] = 0;

        return chk;
ERROR:
        return NULL;
}

void free_checkpoint(struct checkpoint* chk)
{
        if(chk){
                MFREE(chk->base_dir);
                MFREE(chk->base_name);
                MFREE(chk);
                chk = NULL;
        }
}

uint32_t* make_bitvector(uint32_t num_elem)
{
        uint32_t* x = NULL;
        MMALLOC(x , sizeof(int) *((num_elem / 32) + 1));
        RUN(clear_bitvector(x, num_elem));
        return x;
ERROR:
        return NULL;
}

int clear_bitvector(uint32_t* array,uint32_t num_elem)
{
        memset(array, 0, sizeof(int) *((num_elem / 32) + 1));
        return OK;
}

void bit_set(uint32_t* array, uint32_t i)
{
        array[i >> 5] |= (1 << (i & 0x1F));
}

void bit_clr(uint32_t* array, uint32_t i)
{
        array[i >> 5] &= ~(1 << (i & 0x1F));
}

int bit_test(uint32_t* array, uint32_t i)
{
        return (array[i >> 5] & (1 << (i & 0x1F))) != 0 ;
}

int my_file_exists(char* name)
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

uint32_t adler(const void* buf, size_t len)
{
        const uint8_t* buffer = NULL;
        uint32_t s1 = 1;
        uint32_t s2 = 0;

        buffer = (const uint8_t*) buf;

        for (size_t i = 0; i < len; i++) {
                s1 = (s1 + buffer[i]) % 65521;
                s2 = (s2 + s1) % 65521;
        }
        return (s2 << 16) | s1;
}

int ulltoa(uint64_t value, char *buf, int radix)
{
        char tmp[64 + 1];
        char *p1 = tmp, *p2;
        int c;
        static const char xlat[] = "0123456789abcdefghijklmnopqrstuvwxyz";

        if(radix < 2 || radix > 36) {
                return FAIL;
        }
        c = 0;
        do {
                *p1++ = xlat[value % (unsigned)radix];
                c++;
        } while((value /= (unsigned)radix));

        while(c < 64){
                *p1++ = xlat[0];
                c++;
        }

        for(p2 = buf; p1 != tmp; *p2++ = *--p1) {
        }
        *p2 = '\0';
        return OK;
}







char** malloc_2d_char(char**m,int newdim1, int newdim2, char fill_value)
{
        int i,j;

        char** ptr_t = NULL;
        char* ptr_tt = NULL;

        int* int_ptr = NULL;

        int olddim1,olddim2;
        int max1, max2;


        ASSERT((newdim1 > 0), "Malloc 2D char failed: dim1:%d\n",newdim1);
        ASSERT((newdim2 > 0), "Malloc 2D char failed: dim2:%d\n",newdim2);


        if(m == NULL){
                MMALLOC(ptr_t, sizeof(char*) * newdim1);
                MMALLOC(ptr_tt, sizeof(char) * (newdim1*newdim2) + 3*sizeof(int) );

                int_ptr = (int*)ptr_tt;
                int_ptr[0] = 2;
                int_ptr[1] = newdim1;
                int_ptr[2] = newdim2;

                ptr_tt = (char*)(int_ptr + 3);

                for(i = 0;i< newdim1;i++){
                        ptr_t[i] = ptr_tt + i * newdim2;
                        for(j = 0; j < newdim2;j++){
                                ptr_t[i][j] = fill_value;
                        }
                }
                m = ptr_t;

        }else{
                ptr_t = m;
                int_ptr = (int*) m[0];
                int_ptr  = int_ptr -3;

                ptr_tt = (char* )int_ptr;

                olddim1 = *(int_ptr+1);
                olddim2 = *(int_ptr+2);

                DPRINTF3("%d-%d new: %d-%d", olddim1,olddim2, newdim1,newdim2);

                /* in case we want a smaller matrix don't realloc but zero out "free mem"*/
                if(olddim1 >newdim1 || olddim2 > newdim2){
                        max1 = (olddim1 > newdim1) ? olddim1:newdim1;
                        max2 = (olddim2 > newdim2) ? olddim2:newdim2;
                        for(i = 0; i <max1;i++){
                                for(j = 0; j < max2;j++){
                                        if(i >= newdim1 || j >= newdim2){
                                                m[i][j] = fill_value;
                                        }
                                }
                        }
                        if(olddim1 > newdim1){
                                newdim1 = olddim1;
                        }
                        if(olddim1 >newdim1){
                                newdim2 = olddim2;
                        }

                }

                /* case 0: old == new*/
                if(olddim1 == newdim1 && olddim2 == newdim2){
                        return m;
                }

                /*case 1 : both dimensions increase... */


                if(olddim1 < newdim1 && olddim2 < newdim2){
                        MREALLOC(ptr_t,  sizeof(char*) * newdim1);
                        MREALLOC(ptr_tt,  sizeof(char) * (newdim1*newdim2) + 3*sizeof(int) );
                        int_ptr = (int*)ptr_tt;
                        ptr_tt = (char*) (int_ptr + 3);
                        for(i = olddim1-1; i >= 0;i-- ){
                                for(j = olddim2-1;j >=0;j--){

                                        *(ptr_tt + i* newdim2 + j) =*(ptr_tt + i*olddim2 + j);
                                }
                                for(j = newdim2-1;j >= olddim2;j--){
                                        *(ptr_tt + i* newdim2 + j) = fill_value;
                                }

                        }

                        for(i = olddim1; i < newdim1;i++){
                                for(j = 0; j < newdim2;j++){
                                        *(ptr_tt + i* newdim2 + j) = fill_value;
                                }
                        }
                        for(i = 0;i< newdim1;i++){
                                ptr_t[i] = ptr_tt + i * newdim2;
                        }


                        int_ptr[0] = 2;
                        int_ptr[1] = newdim1;
                        int_ptr[2] = newdim2;
                        m = ptr_t;
                }else if(olddim1 == newdim1 && olddim2 < newdim2){
                        MREALLOC(ptr_t,  sizeof(char*) * newdim1);
                        MREALLOC(ptr_tt,  sizeof(char) * (newdim1*newdim2) + 3*sizeof(int) );
                        int_ptr = (int*)ptr_tt;
                        ptr_tt = (char*) (int_ptr + 3);
                        for(i = olddim1-1; i >= 0;i-- ){
                                for(j = olddim2-1;j >=0;j--){

                                        *(ptr_tt + i* newdim2 + j) =*(ptr_tt + i*olddim2 + j);
                                }
                                for(j = newdim2-1;j >= olddim2;j--){
                                        *(ptr_tt + i* newdim2 + j) = fill_value;
                                }

                        }

                        for(i = 0;i< newdim1;i++){
                                ptr_t[i] = ptr_tt + i * newdim2;
                        }


                        int_ptr[0] = 2;
                        int_ptr[1] = newdim1;
                        int_ptr[2] = newdim2;
                        m = ptr_t;
                }else if(olddim1 < newdim1 && olddim2 == newdim2){
                        MREALLOC(ptr_t,  sizeof(char*) * newdim1);
                        MREALLOC(ptr_tt,  sizeof(char) * (newdim1*newdim2) + 3*sizeof(int) );
                        int_ptr = (int*)ptr_tt;
                        ptr_tt = (char*) (int_ptr + 3);
                        for(i = olddim1; i < newdim1;i++){
                                for(j = 0; j < olddim2;j++){
                                        *(ptr_tt + i* olddim2 + j) = fill_value;
                                }
                        }
                        for(i = olddim1;i< newdim1;i++){
                                ptr_t[i] = ptr_tt + i * olddim2;
                        }
                        for(i = 0;i< newdim1;i++){
                                ptr_t[i] = ptr_tt + i * newdim2;
                        }
                        int_ptr[0] = 2;
                        int_ptr[1] = newdim1;
                        int_ptr[2] = newdim2;
                        m = ptr_t;

                }
        }
        return m;
ERROR:
        MFREE(ptr_t );
        MFREE(ptr_tt);
        return NULL;
}


int** malloc_2d_int(int**m,int newdim1, int newdim2, int fill_value)
{
        int i,j;
        int** ptr_t = NULL;
        int* ptr_tt = NULL;

        int olddim1,olddim2;
        int max1, max2;

        ASSERT((newdim1 > 0), "Malloc 2D int failed: dim1:%d\n",newdim1);
        ASSERT((newdim2 > 0), "Malloc 2D int failed: dim2:%d\n",newdim2);

        if(m == NULL){
                MMALLOC(ptr_t, sizeof(int*) * newdim1);
                MMALLOC(ptr_tt, sizeof(int) * (newdim1*newdim2) + 3*sizeof(int));
                ptr_tt[0] = 2;
                ptr_tt[1] = newdim1;
                ptr_tt[2] = newdim2;

                ptr_tt = ptr_tt + 3;

                for(i = 0;i< newdim1;i++){
                        ptr_t[i] = ptr_tt + i * newdim2;
                        for(j = 0; j < newdim2;j++){
                                ptr_t[i][j] = fill_value;
                        }
                }
                m = ptr_t;

        }else{
                ptr_t = m;
                ptr_tt  =m[0]-3;
                olddim1 = *(ptr_tt+1);
                olddim2 = *(ptr_tt+2);
                DPRINTF3("%d-%d new: %d-%d", olddim1,olddim2, newdim1,newdim2 );

                /* in case we want a smaller matrix don't realloc but zero out "free mem"*/
                if(olddim1 >newdim1 || olddim2 > newdim2){
                        max1 = (olddim1 > newdim1) ? olddim1:newdim1;
                        max2 = (olddim2 > newdim2) ? olddim2:newdim2;
                        for(i = 0; i <max1;i++){
                                for(j = 0; j < max2;j++){
                                        if(i >= newdim1 || j >= newdim2){
                                                m[i][j] = fill_value;
                                        }
                                }
                        }
                        if(olddim1 > newdim1){
                                newdim1 = olddim1;
                        }
                        if(olddim1 >newdim1){
                                newdim2 = olddim2;
                        }

                }

                /* case 0: old == new*/
                if(olddim1 == newdim1 && olddim2 == newdim2){
                        return m;
                }

                /*case 1 : both dimensions increase... */


                if(olddim1 < newdim1 && olddim2 < newdim2){
                        MREALLOC(ptr_t,  sizeof(int*) * newdim1);
                        MREALLOC(ptr_tt, sizeof(int) * (newdim1*newdim2) + 3*sizeof(int));
                        ptr_tt = ptr_tt + 3;
                        for(i = olddim1-1; i >= 0;i-- ){
                                for(j = olddim2-1;j >=0;j--){

                                        *(ptr_tt + i* newdim2 + j) =*(ptr_tt + i*olddim2 + j);
                                }
                                for(j = newdim2-1;j >= olddim2;j--){
                                        *(ptr_tt + i* newdim2 + j) = fill_value;
                                }

                        }

                        for(i = olddim1; i < newdim1;i++){
                                for(j = 0; j < newdim2;j++){
                                        *(ptr_tt + i* newdim2 + j) = fill_value;
                                }
                        }
                        for(i = 0;i< newdim1;i++){
                                ptr_t[i] = ptr_tt + i * newdim2;
                        }


                        ptr_tt = ptr_tt - 3;
                        ptr_tt[0] = 2;
                        ptr_tt[1] = newdim1;
                        ptr_tt[2] = newdim2;
                        m = ptr_t;
                }else if(olddim1 == newdim1 && olddim2 < newdim2){
                        MREALLOC(ptr_t,  sizeof(int*) * newdim1);
                        MREALLOC(ptr_tt, sizeof(int) * (newdim1*newdim2) + 3*sizeof(int));
                        ptr_tt = ptr_tt + 3;
                        for(i = olddim1-1; i >= 0;i-- ){
                                for(j = olddim2-1;j >=0;j--){

                                        *(ptr_tt + i* newdim2 + j) =*(ptr_tt + i*olddim2 + j);
                                }
                                for(j = newdim2-1;j >= olddim2;j--){
                                        *(ptr_tt + i* newdim2 + j) = fill_value;
                                }

                        }

                        for(i = 0;i< newdim1;i++){
                                ptr_t[i] = ptr_tt + i * newdim2;
                        }

                        ptr_tt = ptr_tt - 3;
                        ptr_tt[0] = 2;
                        ptr_tt[1] = olddim1;
                        ptr_tt[2] = newdim2;
                        m = ptr_t;
                }else if(olddim1 < newdim1 && olddim2 == newdim2){
                        MREALLOC(ptr_t,  sizeof(int*) * newdim1);
                        MREALLOC(ptr_tt, sizeof(int) * (newdim1*newdim2) + 3*sizeof(int));
                        ptr_tt = ptr_tt + 3;
                        for(i = olddim1; i < newdim1;i++){
                                for(j = 0; j < olddim2;j++){
                                        *(ptr_tt + i* olddim2 + j) = fill_value;
                                }
                        }
                        for(i = olddim1;i< newdim1;i++){
                                ptr_t[i] = ptr_tt + i * olddim2;
                        }
                        for(i = 0;i< newdim1;i++){
                                ptr_t[i] = ptr_tt + i * newdim2;
                        }


                        ptr_tt = ptr_tt - 3;
                        ptr_tt[0] = 2;
                        ptr_tt[1] = newdim1;
                        ptr_tt[2] = olddim2;
                        m = ptr_t;

                }
        }
        return m;
ERROR:
        MFREE(ptr_t );
        MFREE(ptr_tt);

        return NULL;
}

float** malloc_2d_float(float**m,int newdim1, int newdim2,float fill_value)
{
        int i,j;
        float** ptr_t = NULL;
        float* ptr_tt = NULL;

        int* int_ptr = NULL;

        int olddim1,olddim2;
        int max1, max2;

        ASSERT((newdim1 > 0), "Malloc 2D float failed: dim1:%d\n",newdim1);
        ASSERT((newdim2 > 0), "Malloc 2D float failed: dim2:%d\n",newdim2);


        if(m == NULL){
                MMALLOC(ptr_t, sizeof(float*) * newdim1);
                MMALLOC(ptr_tt, sizeof(float) * (newdim1*newdim2) + 3*sizeof(int) );

                int_ptr = (int*)ptr_tt;
                int_ptr[0] = 2;
                int_ptr[1] = newdim1;
                int_ptr[2] = newdim2;

                ptr_tt = (float*)(int_ptr + 3);

                for(i = 0;i< newdim1;i++){
                        ptr_t[i] = ptr_tt + i * newdim2;
                        for(j = 0; j < newdim2;j++){
                                ptr_t[i][j] = fill_value;
                        }
                }
                m = ptr_t;

        }else{
                ptr_t = m;
                int_ptr = (int*) m[0];
                int_ptr  = int_ptr -3;

                ptr_tt = (float* )int_ptr;

                olddim1 = *(int_ptr+1);
                olddim2 = *(int_ptr+2);
                DPRINTF3("%d-%d new: %d-%d", olddim1,olddim2, newdim1,newdim2 );

                /* in case we want a smaller matrix don't realloc but zero out "free mem"*/
                if(olddim1 >newdim1 || olddim2 > newdim2){
                        max1 = (olddim1 > newdim1) ? olddim1:newdim1;
                        max2 = (olddim2 > newdim2) ? olddim2:newdim2;
                        for(i = 0; i <max1;i++){
                                for(j = 0; j < max2;j++){
                                        if(i >= newdim1 || j >= newdim2){
                                                m[i][j] = fill_value;
                                        }
                                }
                        }
                        if(olddim1 > newdim1){
                                newdim1 = olddim1;
                        }
                        if(olddim1 >newdim1){
                                newdim2 = olddim2;
                        }

                }

                /* case 0: old == new*/
                if(olddim1 == newdim1 && olddim2 == newdim2){
                        return m;
                }

                /*case 1 : both dimensions increase... */

                if(olddim1 < newdim1 && olddim2 < newdim2){
                        MREALLOC(ptr_t,  sizeof(float*) * newdim1);
                        MREALLOC(ptr_tt,  sizeof(float) * (newdim1*newdim2) + 3*sizeof(int) );
                        int_ptr = (int*)ptr_tt;
                        ptr_tt = (float*) (int_ptr + 3);
                        for(i = olddim1-1; i >= 0;i-- ){
                                for(j = olddim2-1;j >=0;j--){

                                        *(ptr_tt + i* newdim2 + j) =*(ptr_tt + i*olddim2 + j);
                                }
                                for(j = newdim2-1;j >= olddim2;j--){
                                        *(ptr_tt + i* newdim2 + j) = fill_value;
                                }

                        }

                        for(i = olddim1; i < newdim1;i++){
                                for(j = 0; j < newdim2;j++){
                                        *(ptr_tt + i* newdim2 + j) = fill_value;
                                }
                        }
                        for(i = 0;i< newdim1;i++){
                                ptr_t[i] = ptr_tt + i * newdim2;
                        }


                        int_ptr[0] = 2;
                        int_ptr[1] = newdim1;
                        int_ptr[2] = newdim2;
                        m = ptr_t;
                }else if(olddim1 == newdim1 && olddim2 < newdim2){
                        MREALLOC(ptr_t,  sizeof(float*) * newdim1);
                        MREALLOC(ptr_tt,  sizeof(float) * (newdim1*newdim2) + 3*sizeof(int) );
                        int_ptr = (int*)ptr_tt;
                        ptr_tt = (float*) (int_ptr + 3);
                        for(i = olddim1-1; i >= 0;i-- ){
                                for(j = olddim2-1;j >=0;j--){

                                        *(ptr_tt + i* newdim2 + j) =*(ptr_tt + i*olddim2 + j);
                                }
                                for(j = newdim2-1;j >= olddim2;j--){
                                        *(ptr_tt + i* newdim2 + j) = fill_value;
                                }

                        }

                        for(i = 0;i< newdim1;i++){
                                ptr_t[i] = ptr_tt + i * newdim2;
                        }


                        int_ptr[0] = 2;
                        int_ptr[1] = newdim1;
                        int_ptr[2] = newdim2;
                        m = ptr_t;
                }else if(olddim1 < newdim1 && olddim2 == newdim2){
                        MREALLOC(ptr_t,  sizeof(float*) * newdim1);
                        MREALLOC(ptr_tt,  sizeof(float) * (newdim1*newdim2) + 3*sizeof(int) );
                        int_ptr = (int*)ptr_tt;
                        ptr_tt = (float*) (int_ptr + 3);
                        for(i = olddim1; i < newdim1;i++){
                                for(j = 0; j < olddim2;j++){
                                        *(ptr_tt + i* olddim2 + j) = fill_value;
                                }
                        }
                        for(i = olddim1;i< newdim1;i++){
                                ptr_t[i] = ptr_tt + i * olddim2;
                        }
                        for(i = 0;i< newdim1;i++){
                                ptr_t[i] = ptr_tt + i * newdim2;
                        }
                        int_ptr[0] = 2;
                        int_ptr[1] = newdim1;
                        int_ptr[2] = newdim2;
                        m = ptr_t;
                }
        }
        return m;

ERROR:
        MFREE(ptr_t );
        MFREE(ptr_tt);
        return NULL;
}


float*** malloc_3d_float(int dim1, int dim2, int dim3, float fill_value)
{
        int i,j,c;
        float*** ptr_t = NULL;
        float** ptr_tt = NULL;
        float* ptr_ttt = NULL;

        int* int_ptr = NULL;

        ASSERT((dim1 > 0), "Malloc 3D float failed: dim1:%d\n",dim1);
        ASSERT((dim2 > 0), "Malloc 3D float failed: dim2:%d\n",dim2);
        ASSERT((dim3 > 0), "Malloc 3D float failed: dim3:%d\n",dim3);

        MMALLOC(ptr_t, sizeof(float**) * dim1);
        MMALLOC(ptr_tt, sizeof(float*) * (dim1*dim2));
        MMALLOC(ptr_ttt, sizeof(float) * (dim1 * dim2 * dim3) + 4*sizeof(int) );

        int_ptr = (int*)ptr_ttt;
        int_ptr[0] = 3;
        int_ptr[1] = dim1;
        int_ptr[2] = dim2;
        int_ptr[3] = dim3;

        ptr_ttt = (float*)(int_ptr + 4);

        for (i = 0; i < dim1; i++){
                ptr_t[i] = ptr_tt + dim2*i;
                for (j = 0; j < dim2; j++){
                        ptr_t[i][j] = ptr_ttt + dim2*dim3*i + dim3*j;
                        for(c =0 ; c < dim3;c++){
                                ptr_t[i][j][c] = fill_value;
                        }
                }
        }

        return ptr_t;
ERROR:
        MFREE(ptr_t );
        MFREE(ptr_tt);
        MFREE(ptr_ttt);
        return NULL;

}

float**** malloc_4d_float(int dim1, int dim2, int dim3,int dim4, float fill_value)
{
        int i,j,c,f;
        float**** ptr_t = NULL;
        float*** ptr_tt = NULL;
        float** ptr_ttt = NULL;
        float* ptr_tttt = NULL;

        int* int_ptr = NULL;

        ASSERT((dim1 > 0), "Malloc 3D float failed: dim1:%d\n",dim1);
        ASSERT((dim2 > 0), "Malloc 3D float failed: dim2:%d\n",dim2);
        ASSERT((dim3 > 0), "Malloc 3D float failed: dim3:%d\n",dim3);
        ASSERT((dim4 > 0), "Malloc 3D float failed: dim4:%d\n",dim4);

        MMALLOC(ptr_t, sizeof(float***) * (long int)dim1);
        MMALLOC(ptr_tt, sizeof(float**) * ((long int)dim1*(long int)dim2));
        MMALLOC(ptr_ttt, sizeof(float*) * ((long int)dim1 * (long int)dim2 * (long int)dim3));
        MMALLOC(ptr_tttt, sizeof(float) * ((long int)dim1 * (long int)dim2 * (long int)dim3 *(long int)dim4) + 5*sizeof(int) );

        int_ptr = (int*)ptr_tttt;
        int_ptr[0] = 4;
        int_ptr[1] = dim1;
        int_ptr[2] = dim2;
        int_ptr[3] = dim3;
        int_ptr[4] = dim4;
        ptr_tttt = (float*)(int_ptr + 5);

        for (i = 0; i < dim1; i++){
                ptr_t[i] = ptr_tt + dim2*i;
                for (j = 0; j < dim2; j++){
                        ptr_t[i][j] = ptr_ttt + dim2*dim3*i + dim3*j;
                        for(c =0 ; c < dim3;c++){
                                ptr_t[i][j][c] =  ptr_tttt + dim2*dim3*dim4 * i + dim3*dim4*j  + dim4 *c;
                                for(f =0 ; f < dim4;f++){
                                        ptr_t[i][j][c][f] = fill_value;
                                        fprintf(stdout,"%d %d %d %d: %f\n",i,j,c,f,fill_value );
                                }
                        }
                }
        }

        return ptr_t;
ERROR:
        MFREE(ptr_t );
        MFREE(ptr_tt);
        MFREE(ptr_ttt);
        MFREE(ptr_tttt);
        return NULL;
}

double** malloc_2d_double(double**m,int newdim1, int newdim2, double fill_value)
{
        int i,j;
        double** ptr_t = NULL;
        double* ptr_tt = NULL;

        int* int_ptr = NULL;
        int max1, max2;
        int olddim1,olddim2;

        ASSERT((newdim1 > 0), "Malloc 2D double failed: dim1:%d\n",newdim1);
        ASSERT((newdim2 > 0), "Malloc 2D double failed: dim2:%d\n",newdim2);

        if(m == NULL){
                MMALLOC(ptr_t, sizeof(double*) * newdim1);
                MMALLOC(ptr_tt, sizeof(double) * (newdim1*newdim2) + 3*sizeof(int) );

                int_ptr = (int*)ptr_tt;
                int_ptr[0] = 2;
                int_ptr[1] = newdim1;
                int_ptr[2] = newdim2;

                ptr_tt = (double*)(int_ptr + 3);

                for(i = 0;i< newdim1;i++){
                        ptr_t[i] = ptr_tt + i * newdim2;
                        for(j = 0; j < newdim2;j++){
                                ptr_t[i][j] = fill_value;
                        }
                }
                m = ptr_t;

        }else{
                ptr_t = m;
                int_ptr = (int*) m[0];
                int_ptr  = int_ptr -3;

                ptr_tt = (double* )int_ptr;

                olddim1 = *(int_ptr+1);
                olddim2 = *(int_ptr+2);

                /* in case we want a smaller matrix don't realloc but zero out "free mem"*/
                if(olddim1 >newdim1 || olddim2 > newdim2){
                        max1 = (olddim1 > newdim1) ? olddim1:newdim1;
                        max2 = (olddim2 > newdim2) ? olddim2:newdim2;
                        for(i = 0; i <max1;i++){
                                for(j = 0; j < max2;j++){
                                        if(i >= newdim1 || j >= newdim2){
                                                m[i][j] = fill_value;
                                        }
                                }
                        }
                        if(olddim1 > newdim1){
                                newdim1 = olddim1;
                        }
                        if(olddim1 >newdim1){
                                newdim2 = olddim2;
                        }

                }


                /* case 0: old == new*/
                if(olddim1 == newdim1 && olddim2 == newdim2){
                        return m;
                }

                /*case 1 : both dimensions increase... */


                if(olddim1 < newdim1 && olddim2 < newdim2){
                        MREALLOC(ptr_t,  sizeof(double*) * newdim1);
                        MREALLOC(ptr_tt,  sizeof(double) * (newdim1*newdim2) + 3*sizeof(int) );
                        int_ptr = (int*)ptr_tt;
                        ptr_tt = (double*) (int_ptr + 3);
                        for(i = olddim1-1; i >= 0;i-- ){
                                for(j = olddim2-1;j >=0;j--){

                                        *(ptr_tt + i* newdim2 + j) =*(ptr_tt + i*olddim2 + j);
                                }
                                for(j = newdim2-1;j >= olddim2;j--){
                                        *(ptr_tt + i* newdim2 + j) = fill_value;
                                }

                        }

                        for(i = olddim1; i < newdim1;i++){
                                for(j = 0; j < newdim2;j++){
                                        *(ptr_tt + i* newdim2 + j) = fill_value;
                                }
                        }
                        for(i = 0;i< newdim1;i++){
                                ptr_t[i] = ptr_tt + i * newdim2;
                        }


                        int_ptr[0] = 2;
                        int_ptr[1] = newdim1;
                        int_ptr[2] = newdim2;
                        m = ptr_t;
                }else if(olddim1 == newdim1 && olddim2 < newdim2){
                        MREALLOC(ptr_t,  sizeof(double*) * newdim1);
                        MREALLOC(ptr_tt,  sizeof(double) * (newdim1*newdim2) + 3*sizeof(int) );
                        int_ptr = (int*)ptr_tt;
                        ptr_tt = (double*) (int_ptr + 3);
                        for(i = olddim1-1; i >= 0;i-- ){
                                for(j = olddim2-1;j >=0;j--){

                                        *(ptr_tt + i* newdim2 + j) =*(ptr_tt + i*olddim2 + j);
                                }
                                for(j = newdim2-1;j >= olddim2;j--){
                                        *(ptr_tt + i* newdim2 + j) = fill_value;
                                }

                        }

                        for(i = 0;i< newdim1;i++){
                                ptr_t[i] = ptr_tt + i * newdim2;
                        }


                        int_ptr[0] = 2;
                        int_ptr[1] = newdim1;
                        int_ptr[2] = newdim2;
                        m = ptr_t;
                }else if(olddim1 < newdim1 && olddim2 == newdim2){
                        MREALLOC(ptr_t,  sizeof(double*) * newdim1);
                        MREALLOC(ptr_tt,  sizeof(double) * (newdim1*newdim2) + 3*sizeof(int) );
                        int_ptr = (int*)ptr_tt;
                        ptr_tt = (double*) (int_ptr + 3);
                        for(i = olddim1; i < newdim1;i++){
                                for(j = 0; j < olddim2;j++){
                                        *(ptr_tt + i* olddim2 + j) = fill_value;
                                }
                        }
                        for(i = olddim1;i< newdim1;i++){
                                ptr_t[i] = ptr_tt + i * olddim2;
                        }
                        for(i = 0;i< newdim1;i++){
                                ptr_t[i] = ptr_tt + i * newdim2;
                        }
                        int_ptr[0] = 2;
                        int_ptr[1] = newdim1;
                        int_ptr[2] = newdim2;
                        m = ptr_t;

                }
        }
        return m;
ERROR:
        MFREE(ptr_t );
        MFREE(ptr_tt);

        return NULL;
}



void free_2d(void** m)
{
        int* ptr_tt  = (int*)m[0];
        ptr_tt = ptr_tt -3;
        MFREE(ptr_tt);
        MFREE(m);
}

void free_3d(void*** m)
{
        int* ptr_ttt  = (int*)m[0][0];
        ptr_ttt = ptr_ttt -4;
        MFREE(ptr_ttt);
        MFREE(m[0]);
        MFREE(m);
}

void free_4d(void**** m)
{
        int* ptr_tttt  = (int*)m[0][0][0];
        ptr_tttt = ptr_tttt -5;
        MFREE(ptr_tttt);
        MFREE(m[0][0]);
        MFREE(m[0]);
        MFREE(m);
}



uint16_t prob_to_uint16(float x)
{
        uint16_t a;

        x = scaledprob2prob(x);
        if(x > 1.0){
                x = 1.0;
        }
        if(x < 0.0){
                x = 0.0;
        }

        a = (uint16_t) roundf(x *UINT16_MAX);

        return a;
}

float uint16_to_prob(uint16_t a)
{
        float x;

        x = prob2scaledprob( (double)a  / UINT16_MAX);
        return x;
}



uint32_t prob_to_uint32(float x)
{
        uint32_t a;

        x = scaledprob2prob(x);
        if(x > 1.0){
                x = 1.0;
        }
        if(x < 0.0){
                x = 0.0;
        }

        a = (uint32_t) roundf(x *UINT32_MAX);

        return a;
}

float uint32_to_prob(uint32_t a)
{
        float x;

        x = prob2scaledprob( (double)a  / UINT32_MAX);
        return x;
}






void init_logsum()
{
        static int called = 0;
        int i;
        if(!called){
                called = 1;
                for(i = 0; i < LOGSUM_SIZE;i++){
                        logsum_lookup[i] = log(1.0 +exp((double) -i / SCALE));
                }
        }
}

float logsum(const float a,const float b)
{
        register const float max = MACRO_MAX(a, b);
        register const float min = MACRO_MIN(a, b);

        if(min == -INFINITY){
                return max;
        }
        if( (max-min) >= 15.7f){
                return max;
        }
        return  max+ logsum_lookup[(int)((max-min)*SCALE)];
}

float prob2scaledprob(float p)
{
        if(p == 0.0){
                return -INFINITY;
        }else{
                return  log(p);
        }
}


float scaledprob2prob(float p)
{
        if(p == -INFINITY){
                return 0.0;
        }else{
                return exp(p);
        }
}



char* shorten_pathname(char* p)
{
        int i;
        char* tmp = p;
        int len = 0;
        len = (int) strlen(p);
        if(len){
                for(i = 0; i< len;i++){
                        if(p[i] == '/'){
                                tmp = p+i +1;
                        }
                }
        }
        return tmp;
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



int replace_punctuation_with_underscore(char* p)
{
        int i;
        int c;
        int len = 0;
        len = (int) strlen(p);
        if(len){
                for(i = 0; i < len;i++){
                        c = (int) p[i];
                        if(ispunct(c)){
                                p[i] = '_';
                        }
                }
        }
        return OK;
}


#ifdef ITEST

int dummy_broken_func(int i);

int char_test(void);
int int_test(void);
int float_test(void);
int double_test(void);


int float_3d_test(void);

int float_4d_test(void);

int float_4d_test(void)
{
        float**** m = NULL;
        int i,j,c,f;
        int dim1,dim2,dim3,dim4;
        dim1 = 2;

        dim2 = 3;
        dim3 = 5;
        dim4 = 5;

        RUNP(m = malloc_4d_float(dim1,dim2,dim3,dim4,0.0f));

        for(i =0;i < dim1;i++){
                for(j = 0; j < dim2;j++){
                        for(c = 0; c < dim3;c++){
                                for(f = 0; f < dim4;f++){
                                        m[i][j][c][f] =(float)i/1000.0 + (float)j/10.0+  (float)c +  (float)f;
                                }
                        }
                }
        }

        for(i =0;i < dim1;i++){
                fprintf(stdout,"LEVEL: %d\n",i);

                for(j = 0; j < dim2;j++){
                        fprintf(stdout,"SUBLEVEL: %d\n",j);
                        for(c = 0; c < dim3;c++){
                                for(f = 0; f < dim4;f++){
                                        fprintf(stdout," %0.1f",m[i][j][c][f]);
                                }
                                fprintf(stdout,"\n");
                        }
                        fprintf(stdout,"\n");
                }
                fprintf(stdout,"\n");
        }
        fprintf(stdout,"\n");


        free_4d((void****)m);
        return OK;
ERROR:
        return FAIL;
}


int float_3d_test(void)
{
        float*** m = NULL;
        int i,j,c;
        int dim1,dim2,dim3;
        dim1 = 3;

        dim2 = 5;
        dim3 = 5;

        RUNP(m = malloc_3d_float(dim1,dim2,dim3,0.0f));
        for(i =0;i < dim1;i++){
                for(j = 0; j < dim2;j++){
                        for(c = 0; c < dim3;c++){
                                m[i][j][c] = (float)i/10.0+  (float)j +  (float)c;
                        }
                }
        }

        for(i =0;i < dim1;i++){
                fprintf(stdout,"LEVEL: %d\n",i);

                for(j = 0; j < dim2;j++){
                        for(c = 0; c < dim3;c++){
                                fprintf(stdout," %0.1f",m[i][j][c]);
                        }
                        fprintf(stdout,"\n");
                }
                fprintf(stdout,"\n");
        }
        fprintf(stdout,"\n");


        free_3d((void***)m);
        return OK;
ERROR:
        return FAIL;
}

int char_test(void)
{
        int dim1,dim2,i,j;
        char** m = NULL;

        DPRINTF1("Testing char 2D");

        dim1 = 5;
        dim2 = 5;

        RUNP(m = malloc_2d_char(m,dim1,dim2,0));

        if(!m){
                ERROR_MSG("malloc_2d_char failed");
        }
        for(i =0;i < dim1;i++){
                for(j = 0; j < dim2;j++){
                        m[i][j] = i+j + 65;
                }
        }
        for(i =0;i < dim1;i++){
                for(j = 0; j < dim2;j++){
                        if(m[i][j] == 0){
                                fprintf(stdout," %d",m[i][j]);
                        }else{
                                fprintf(stdout," %c",m[i][j]);
                        }
                }
                fprintf(stdout,"\n");
        }
        fprintf(stdout,"\n");
        dim1 += 1;
        dim2 += 1;
        RUNP(m = malloc_2d_char(m,dim1,dim2,0));
        for(i =0;i < dim1;i++){
                for(j = 0; j < dim2;j++){
                        if(m[i][j] == 0){
                                fprintf(stdout," %d",m[i][j]);
                        }else{
                                fprintf(stdout," %c",m[i][j]);
                        }
                }
                fprintf(stdout,"\n");
        }
        fprintf(stdout,"\n");
        dim1 += 1;
        dim2 += 0;
        RUNP(m = malloc_2d_char(m,dim1,dim2,0));
        for(i =0;i < dim1;i++){
                for(j = 0; j < dim2;j++){
                        if(m[i][j] == 0){
                                fprintf(stdout," %d",m[i][j]);
                        }else{
                                fprintf(stdout," %c",m[i][j]);
                        }
                }
                fprintf(stdout,"\n");
        }
        fprintf(stdout,"\n");

        dim1 += 0;
        dim2 += 1;
        RUNP(m = malloc_2d_char(m,dim1,dim2,0));
        for(i =0;i < dim1;i++){
                for(j = 0; j < dim2;j++){
                        if(m[i][j] == 0){
                                fprintf(stdout," %d",m[i][j]);
                        }else{
                                fprintf(stdout," %c",m[i][j]);
                        }
                }
                fprintf(stdout,"\n");
        }
        fprintf(stdout,"\n");

        dim1 = 3;
        dim2 = 3;
        RUNP(m = malloc_2d_char(m,dim1,dim2,0));
        for(i =0;i < dim1;i++){
                for(j = 0; j < dim2;j++){
                        if(m[i][j] == 0){
                                fprintf(stdout," %d",m[i][j]);
                        }else{
                                fprintf(stdout," %c",m[i][j]);
                        }
                }
                fprintf(stdout,"\n");
        }
        fprintf(stdout,"\n");


        dim1 = 10;
        dim2 = 10;
        RUNP(m = malloc_2d_char(m,dim1,dim2,0));
        for(i =0;i < dim1;i++){
                for(j = 0; j < dim2;j++){
                        if(m[i][j] == 0){
                                fprintf(stdout," %d",m[i][j]);
                        }else{
                                fprintf(stdout," %c",m[i][j]);
                        }
                }
                fprintf(stdout,"\n");
        }
        fprintf(stdout,"\n");

        free_2d((void**) m);
        return OK;
ERROR:
        return FAIL;
}


int int_test(void)
{
        int dim1,dim2,i,j;
        int** m = NULL;

        fprintf(stdout,"Testing int 2D\n");

        dim1 = 5;
        dim2 = 5;
        RUNP(m = malloc_2d_int(m,dim1,dim2,0));
        for(i =0;i < dim1;i++){
                for(j = 0; j < dim2;j++){
                        m[i][j] = i+j;
                }
        }
        for(i =0;i < dim1;i++){
                for(j = 0; j < dim2;j++){
                        fprintf(stdout," %d",m[i][j]);
                }
                fprintf(stdout,"\n");
        }
        fprintf(stdout,"\n");
        dim1 += 1;
        dim2 += 1;
        RUNP(m = malloc_2d_int(m,dim1,dim2,0));
        for(i =0;i < dim1;i++){
                for(j = 0; j < dim2;j++){
                        fprintf(stdout," %d",m[i][j]);
                }
                fprintf(stdout,"\n");
        }
        fprintf(stdout,"\n");
        dim1 += 1;
        dim2 += 0;
        RUNP(m = malloc_2d_int(m,dim1,dim2,0));
        for(i =0;i < dim1;i++){
                for(j = 0; j < dim2;j++){
                        fprintf(stdout," %d",m[i][j]);
                }
                fprintf(stdout,"\n");
        }
        fprintf(stdout,"\n");

        dim1 += 0;
        dim2 += 1;
        RUNP(m = malloc_2d_int(m,dim1,dim2,0));
        for(i =0;i < dim1;i++){
                for(j = 0; j < dim2;j++){
                        fprintf(stdout," %d",m[i][j]);
                }
                fprintf(stdout,"\n");
        }
        fprintf(stdout,"\n");


        dim1 = 3;
        dim2 = 3;
        RUNP(m = malloc_2d_int(m,dim1,dim2,0));
        for(i =0;i < dim1;i++){
                for(j = 0; j < dim2;j++){
                        fprintf(stdout," %d",m[i][j]);
                }
                fprintf(stdout,"\n");
        }
        fprintf(stdout,"\n");

        dim1 = 10;
        dim2 = 10;
        RUNP(m = malloc_2d_int(m,dim1,dim2,0));
        for(i =0;i < dim1;i++){
                for(j = 0; j < dim2;j++){
                        fprintf(stdout," %d",m[i][j]);
                }
                fprintf(stdout,"\n");
        }
        fprintf(stdout,"\n");

        free_2d((void**) m);
        return OK;
ERROR:
        return FAIL;
}


int float_test(void)
{
        int dim1,dim2,i,j;
        float** m = NULL;

        DPRINTF1("Testing float 2D");

        dim1 = 5;
        dim2 = 5;
        RUNP(m = malloc_2d_float(m,dim1,dim2,0.0f));
        for(i =0;i < dim1;i++){
                for(j = 0; j < dim2;j++){
                        m[i][j] = i+j + 0.1;
                }
        }
        for(i =0;i < dim1;i++){
                for(j = 0; j < dim2;j++){

                        fprintf(stdout," %0.1f",m[i][j]);

                }
                fprintf(stdout,"\n");
        }
        fprintf(stdout,"\n");
        dim1 += 1;
        dim2 += 1;
        RUNP(m = malloc_2d_float(m,dim1,dim2,0.0f));
        for(i =0;i < dim1;i++){
                for(j = 0; j < dim2;j++){
                        fprintf(stdout," %0.1f",m[i][j]);
                }
                fprintf(stdout,"\n");
        }
        fprintf(stdout,"\n");
        dim1 += 1;
        dim2 += 0;
        RUNP(m = malloc_2d_float(m,dim1,dim2,0.0f));
        for(i =0;i < dim1;i++){
                for(j = 0; j < dim2;j++){
                        fprintf(stdout," %0.1f",m[i][j]);
                }
                fprintf(stdout,"\n");
        }
        fprintf(stdout,"\n");

        dim1 += 0;
        dim2 += 1;
        RUNP(m = malloc_2d_float(m,dim1,dim2,0.0f));
        for(i =0;i < dim1;i++){
                for(j = 0; j < dim2;j++){
                        fprintf(stdout," %0.1f",m[i][j]);
                }
                fprintf(stdout,"\n");
        }
        fprintf(stdout,"\n");

        dim1 = 3;
        dim2 = 3;
        RUNP(m = malloc_2d_float(m,dim1,dim2,0.0f));
        for(i =0;i < dim1;i++){
                for(j = 0; j < dim2;j++){
                        fprintf(stdout," %0.1f",m[i][j]);
                }
                fprintf(stdout,"\n");
        }
        fprintf(stdout,"\n");

        dim1 = 10;
        dim2 = 10;
        RUNP(m = malloc_2d_float(m,dim1,dim2,0.0f));
        for(i =0;i < dim1;i++){
                for(j = 0; j < dim2;j++){
                        fprintf(stdout," %0.1f",m[i][j]);
                }
                fprintf(stdout,"\n");
        }
        fprintf(stdout,"\n");

        free_2d((void**) m);
        return OK;
ERROR:
        return FAIL;
}


int double_test(void)
{
        int dim1,dim2,i,j;
        double** m = NULL;

        DPRINTF1("Testing double 2D");

        dim1 = 5;
        dim2 = 5;

        RUNP(m = malloc_2d_double(m,dim1,dim2,0.0));
        for(i =0;i < dim1;i++){
                for(j = 0; j < dim2;j++){
                        m[i][j] = i+j + 0.1;
                }
        }
        for(i =0;i < dim1;i++){
                for(j = 0; j < dim2;j++){

                        fprintf(stdout," %0.1f",m[i][j]);

                }
                fprintf(stdout,"\n");
        }
        fprintf(stdout,"\n");
        dim1 += 1;
        dim2 += 1;
        RUNP(m = malloc_2d_double(m,dim1,dim2,0.0));
        for(i =0;i < dim1;i++){
                for(j = 0; j < dim2;j++){
                        fprintf(stdout," %0.1f",m[i][j]);
                }
                fprintf(stdout,"\n");
        }
        fprintf(stdout,"\n");
        dim1 += 1;
        dim2 += 0;
        RUNP(m = malloc_2d_double(m,dim1,dim2,0.0));
        for(i =0;i < dim1;i++){
                for(j = 0; j < dim2;j++){
                        fprintf(stdout," %0.1f",m[i][j]);
                }
                fprintf(stdout,"\n");
        }
        fprintf(stdout,"\n");

        dim1 += 0;
        dim2 += 1;
        RUNP(m = malloc_2d_double(m,dim1,dim2,0.0));
        for(i =0;i < dim1;i++){
                for(j = 0; j < dim2;j++){
                        fprintf(stdout," %0.1f",m[i][j]);
                }
                fprintf(stdout,"\n");
        }
        fprintf(stdout,"\n");

        dim1 = 3;
        dim2 = 3;
        RUNP(m = malloc_2d_double(m,dim1,dim2,0.0));
        for(i =0;i < dim1;i++){
                for(j = 0; j < dim2;j++){
                        fprintf(stdout," %0.1f",m[i][j]);
                }
                fprintf(stdout,"\n");
        }
        fprintf(stdout,"\n");

        dim1 = 10;
        dim2 = 10;
        RUNP(m = malloc_2d_double(m,dim1,dim2,0.0));
        for(i =0;i < dim1;i++){
                for(j = 0; j < dim2;j++){
                        fprintf(stdout," %0.1f",m[i][j]);
                }
                fprintf(stdout,"\n");
        }
        fprintf(stdout,"\n");


        free_2d((void**) m);
        return OK;
ERROR:
        return FAIL;
}

struct char_struct{
        int len;
        char* seq;
};


int main (int argc,char * const argv[])
{
        int* p = NULL;
        MFREE(p);
        MFREE(p);
        char* my_str_seq  = NULL;
        struct char_struct* my_str = NULL;
        print_program_header(argv, "MEMORY TESTING PROGRAM");

        fprintf(stdout,"Testing logging functions\n");
        LOG_MSG("Testing %s","the here and now.");
        log_message("TESTING");
        error(AT,"test error %d in stderr %s", 1,"this should go to stderr" );
        message(AT, "Here is a message %d", 100);
        warning(AT, "Oh dear I have a warning: %s","big warning to stdout");

        //set_logfile("testxt");

        //et_logfile("00");
        int** gaps = NULL;

        gaps = galloc(gaps,10,10,0);

        gfree(gaps);
        message(AT, "Here is a message %d", 666);
        log_message("TESTING TEE");

        error(AT,"test error %d in stderr %s", 2,"this should be tee'd" );
        //error(AT,"test TEE ERROR " );


        warning(AT, "Oh dear I have a warning: %s","big warning to be tee'd");

        //set_logfile(NULL);

        warning(AT, "Oh dear I have a warning (not to be written to log...) : %s","big warning to be tee'd");
        ASSERT(1 == 1,"Of dear 1 in NOT equal to %d", 1);



        fprintf(stdout,"Running libks sanity tests\n");

        RUN(char_test() );
        char_test();
        int_test();
        float_test();
        double_test();

        float_3d_test();
        float_4d_test();

        p= NULL;

        MCALLOC(p, 10,int);

        MFREE(p);
        float* fffffflf = NULL;
        fffffflf= galloc(fffffflf,10);
        gfree(fffffflf);
        MCALLOC(my_str, 1,struct char_struct);

        fprintf(stderr,"%d %s\n",	my_str->len,my_str->seq);

        MMALLOC(my_str_seq, sizeof(char) * 10);
        my_str->seq =my_str_seq;

        log_message("All is good");
        MFREE(my_str->seq);
        MFREE(my_str);


        log_message("Yes %s is ","it");

        char* cmd = NULL;
        RUNP(cmd =  make_cmd_line(argc,argv));
        DECLARE_CHK(MAIN,".");

        RUN_CHECKPOINT(MAIN,float_4d_test(),cmd);

        DESTROY_CHK(MAIN);
        MFREE(cmd);

        return EXIT_SUCCESS;
ERROR:
        MFREE(my_str_seq);
        MFREE(my_str);
        MFREE(p);
        return EXIT_FAILURE;
}

int dummy_broken_func(int i)
{
        return FAIL;
}

#endif
