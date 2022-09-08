/* Tracking cpu/system/elapsed time used by a process.
 *
 * SRE, Wed Feb 22 19:30:36 2006 [St. Louis] [moved to Easel]
 * SRE, Thu Aug  3 08:00:35 2000 [St. Louis] [moved to SQUID]
 * SRE, Fri Nov 26 14:54:21 1999 [St. Louis] [HMMER]
 */
#ifndef eslSTOPWATCH_INCLUDED
#define eslSTOPWATCH_INCLUDED

#include <stdio.h>
#ifdef ESL_STOPWATCH_IMPORT
   #define EXTERN
#else
   #ifndef EXTERN
      #ifdef __cplusplus
         #define EXTERN extern "C"
      #else
         #define EXTERN extern
      #endif
   #endif
#endif

#include <time.h>
#ifdef HAVE_TIMES
#include <sys/times.h>
#endif
#ifdef HAVE_UNISTD_H
#include <unistd.h>		/* need for sysconf() */
#endif

typedef struct {
#ifdef eslSTOPWATCH_HIGHRES
        double     t0;                /* baseline wall time from Nadeau routine */
#elif  HAVE_TIMES
        clock_t    t0;		/* baseline wall time, POSIX times()      */
#else
        time_t     t0;                /* baseline wall time from ANSI time()    */
#endif

#ifdef HAVE_TIMES
        struct tms cpu0;		/* baseline CPU/system time, POSIX times()      */
#else
        clock_t cpu0;			/* baseline CPU time, fallback to ANSI clock()  */
#endif

        /* elapsed/user/sys are t-t0 results for the last time the
         * watch was Stop()'ed.
         */
        double elapsed;               /* elapsed wall time, seconds */
        double user;                  /* CPU time, seconds          */
        double sys;                   /* system time, seconds       */
} ESL_STOPWATCH;


EXTERN ESL_STOPWATCH *esl_stopwatch_Create(void);
EXTERN void           esl_stopwatch_Destroy(ESL_STOPWATCH *w);

EXTERN int esl_stopwatch_Start(ESL_STOPWATCH *w);
EXTERN int esl_stopwatch_Stop(ESL_STOPWATCH *w);
EXTERN int esl_stopwatch_Display(FILE *fp, ESL_STOPWATCH *w, char *prefix);
EXTERN int tl_stopwatch_Display(ESL_STOPWATCH *w);
EXTERN double esl_stopwatch_GetElapsed(ESL_STOPWATCH *w);
EXTERN double  tl_stopwatch_utime(ESL_STOPWATCH *w);
EXTERN int esl_stopwatch_Include(ESL_STOPWATCH *master, ESL_STOPWATCH *w);


#define DECLARE_TIMER(n)  ESL_STOPWATCH* timer_##n = esl_stopwatch_Create();

#define START_TIMER(n) esl_stopwatch_Start(timer_##n);
#define STOP_TIMER(n)  esl_stopwatch_Stop(timer_##n);
#define GET_TIMING(n) tl_stopwatch_Display(timer_##n);

#define GET_ELAPSED(n) esl_stopwatch_GetElapsed(timer_##n);
#define GET_USERTIME(n) tl_stopwatch_utime(timer_##n);

#define DESTROY_TIMER(n)  esl_stopwatch_Destroy(timer_##n);

#undef ESL_STOPWATCH_IMPORT
#undef EXTERN

#endif /*eslSTOPWATCH_INCLUDED*/
