
#include "rng.h"

/* code here was adopted from:  */

/* http://vigna.di.unimi.it/xorshift/xoshiro256starstar.c */
/* http://xoshiro.di.unimi.it/splitmix64.c */

/* written by: */
/* Sebastiano Vigna (vigna@acm.org) */
/* David Blackman */

#ifdef ITESTRNG
int main(int argc, char *argv[])
{
        struct rng_state* rng = NULL;
        int i;
        RUNP(rng = init_rng(0));
        for(i = 0; i < 10;i++){
                fprintf(stdout,"%f\t%d\n", tl_random_double(rng), tl_random_int(rng,10));
        }
        free_rng(rng);
        return EXIT_SUCCESS;
ERROR:
        return EXIT_FAILURE;
}
#endif

static inline uint64_t rotl(const uint64_t x, int k);
uint64_t next(struct rng_state* s);
void jump(struct rng_state* s);
void long_jump(struct rng_state* s);

static uint64_t choose_arbitrary_seed(void);
static uint32_t jenkins_mix3(uint32_t a, uint32_t b, uint32_t c);

double tl_random_double(struct rng_state* rng)
{
        uint64_t x;
        double y;
        do{
                x = next(rng);
                y = ((double) x / 18446744073709551616.0);
        }while (y == 0.0);
        return y;
}

int tl_random_int(struct rng_state* rng,int a)
{
        return (int) (tl_random_double(rng) * a);
}

struct rng_state* init_rng(uint64_t seed)
{
        struct rng_state* s = NULL;
        uint64_t z;
        uint64_t sanity;

        MMALLOC(s, sizeof(struct rng_state));
        if(!seed){
                seed = choose_arbitrary_seed();
        }
        sanity = 0;
        while(!sanity){
                sanity = 0;
                z = (seed += 0x9e3779b97f4a7c15);
                z = (z ^ (z >> 30)) * 0xbf58476d1ce4e5b9;
                z = (z ^ (z >> 27)) * 0x94d049bb133111eb;
                s->s[0] = z ^ (z >> 31);
                if(s->s[0]){
                        sanity++;
                }
                z = (seed += 0x9e3779b97f4a7c15);
                z = (z ^ (z >> 30)) * 0xbf58476d1ce4e5b9;
                z = (z ^ (z >> 27)) * 0x94d049bb133111eb;
                s->s[1] = z ^ (z >> 31);
                if(s->s[1]){
                        sanity++;
                }

                z = (seed += 0x9e3779b97f4a7c15);
                z = (z ^ (z >> 30)) * 0xbf58476d1ce4e5b9;
                z = (z ^ (z >> 27)) * 0x94d049bb133111eb;
                s->s[2] = z ^ (z >> 31);
                if(s->s[2]){
                        sanity++;
                }

                z = (seed += 0x9e3779b97f4a7c15);
                z = (z ^ (z >> 30)) * 0xbf58476d1ce4e5b9;
                z = (z ^ (z >> 27)) * 0x94d049bb133111eb;
                s->s[3] = z ^ (z >> 31);
                if(s->s[3]){
                        sanity++;
                }
        }
        return s;
ERROR:
        return NULL;
}

void free_rng(struct rng_state* rng)
{
        if(rng){
                MFREE(rng);
        }
}
/* Taken from easel library (by Sean Eddy) */
static uint64_t choose_arbitrary_seed(void)
{
        uint32_t a = (uint32_t) time ((time_t *) NULL);
        uint32_t b = 87654321;	 // we'll use getpid() below, if we can
        uint32_t c = (uint32_t) clock();  // clock() gives time since process invocation, in msec at least, if not usec
        uint64_t seed;
#ifdef HAVE_GETPID
        b  = (uint32_t) getpid();	 // preferable b choice, if we have POSIX getpid()
#endif
        seed = jenkins_mix3(a,b,c);    // try to decorrelate closely spaced choices of pid/times
        return (seed == 0) ? 42 : seed; /* 42 is entirely arbitrary, just to avoid seed==0. */
}

/* jenkins_mix3()
 *
 * from Bob Jenkins: given a,b,c, generate a number that's distributed
 * reasonably uniformly on the interval 0..2^32-1 even for closely
 * spaced choices of a,b,c.
 */
static uint32_t jenkins_mix3(uint32_t a, uint32_t b, uint32_t c)
{
        a -= b; a -= c; a ^= (c>>13);
        b -= c; b -= a; b ^= (a<<8);
        c -= a; c -= b; c ^= (b>>13);
        a -= b; a -= c; a ^= (c>>12);
        b -= c; b -= a; b ^= (a<<16);
        c -= a; c -= b; c ^= (b>>5);
        a -= b; a -= c; a ^= (c>>3);
        b -= c; b -= a; b ^= (a<<10);
        c -= a; c -= b; c ^= (b>>15);
        return c;
}


static inline uint64_t rotl(const uint64_t x, int k) {
        return (x << k) | (x >> (64 - k));
}




uint64_t next(struct rng_state* s)
{
        const uint64_t result_starstar = rotl(s->s[1] * 5, 7) * 9;
        const uint64_t t = s->s[1] << 17;

        s->s[2] ^= s->s[0];
        s->s[3] ^= s->s[1];
        s->s[1] ^= s->s[2];
        s->s[0] ^= s->s[3];

        s->s[2] ^= t;

        s->s[3] = rotl(s->s[3], 45);

        return result_starstar;
}


/* This is the jump function for the generator. It is equivalent
   to 2^128 calls to next(); it can be used to generate 2^128
   non-overlapping subsequences for parallel computations. */

void jump(struct rng_state* s)
{
        static const uint64_t JUMP[] = { 0x180ec6d33cfd0aba,
                                         0xd5a61266f0c9392c,
                                         0xa9582618e03fc9aa,
                                         0x39abdc4529b1661c };

        uint64_t s0 = 0;
        uint64_t s1 = 0;
        uint64_t s2 = 0;
        uint64_t s3 = 0;
        for(int i = 0; i < sizeof JUMP / sizeof *JUMP; i++)
                for(int b = 0; b < 64; b++) {
                        if (JUMP[i] & UINT64_C(1) << b) {
                                s0 ^= s->s[0];
                                s1 ^= s->s[1];
                                s2 ^= s->s[2];
                                s3 ^= s->s[3];
                        }
                        next(s);
                }

        s->s[0] = s0;
        s->s[1] = s1;
        s->s[2] = s2;
        s->s[3] = s3;
}



/* This is the long-jump function for the generator. It is equivalent to
   2^192 calls to next(); it can be used to generate 2^64 starting points,
   from each of which jump() will generate 2^64 non-overlapping
   subsequences for parallel distributed computations. */

void long_jump(struct rng_state* s)
{
        static const uint64_t LONG_JUMP[] = { 0x76e15d3efefdcbbf,
                                              0xc5004e441c522fb3,
                                              0x77710069854ee241,
                                              0x39109bb02acbe635 };

        uint64_t s0 = 0;
        uint64_t s1 = 0;
        uint64_t s2 = 0;
        uint64_t s3 = 0;
        for(int i = 0; i < sizeof LONG_JUMP / sizeof *LONG_JUMP; i++)
                for(int b = 0; b < 64; b++) {
                        if (LONG_JUMP[i] & UINT64_C(1) << b) {
                                s0 ^= s->s[0];
                                s1 ^= s->s[1];
                                s2 ^= s->s[2];
                                s3 ^= s->s[3];
                        }
                        next(s);
                }

        s->s[0] = s0;
        s->s[1] = s1;
        s->s[2] = s2;
        s->s[3] = s3;
}
