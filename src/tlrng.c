
#include <time.h>
#include <sys/types.h>
#include <unistd.h>
/* #include <math.h> */

#define M_PI        3.14159265358979323846

#define M_2PI 2.0*M_PI

#include "tldevel.h"

#define TLRNG_IMPORT
#include "tlrng.h"

/* code here was adopted from:  */

/*  Written in 2018 by David Blackman and Sebastiano Vigna (vigna@acm.org)

To the extent possible under law, the author has dedicated all copyright
and related and neighboring rights to this software to the public domain
worldwide. This software is distributed without any warranty.

See <http://creativecommons.org/publicdomain/zero/1.0/>. */

/* http://vigna.di.unimi.it/xorshift/xoshiro256starstar.c */
/* http://xoshiro.di.unimi.it/splitmix64.c */

/* written by: */
/* Sebastiano Vigna (vigna@acm.org) */
/* David Blackman */

/* code for sampling from varius distributions was taken from  */

/* Copyright 2005 Robert Kern (robert.kern@gmail.com)
 *
 * Permission is hereby granted, free of charge, to any person obtaining a
 * copy of this software and associated documentation files (the
 * "Software"), to deal in the Software without restriction, including
 * without limitation the rights to use, copy, modify, merge, publish,
 * distribute, sublicense, and/or sell copies of the Software, and to
 * permit persons to whom the Software is furnished to do so, subject to
 * the following conditions:
 *
 * The above copyright notice and this permission notice shall be included
 * in all copies or substantial portions of the Software.
 *
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS
 * OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF
 * MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT.
 * IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY
 * CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT,
 * TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE
 * SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.
 */


/* And code bits from  */
/* the easel library (by Sean Eddy) */


struct rng_state{
        uint64_t s[4];
        uint8_t gen;
        double z1;
        int has_gauss; /* !=0: gauss contains a gaussian deviate */
        double gauss;
};


static inline uint64_t rotl(const uint64_t x, int k);
static uint64_t next(struct rng_state* s);
static void jump(struct rng_state* s);
void long_jump(struct rng_state* s);

static uint64_t choose_arbitrary_seed(void);
static uint32_t jenkins_mix3(uint32_t a, uint32_t b, uint32_t c);


static double tl_standard_gamma(struct rng_state* rng, double shape);
static double tl_standard_exponential(struct rng_state* rng);
static double tl_gauss(struct rng_state* rng);


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

/* from:  */
//https://en.wikipedia.org/wiki/Box%E2%80%93Muller_transform
double tl_random_gaussian(struct rng_state* rng, double mu, double sigma)
{
        rng->gen = !rng->gen;

        if (!rng->gen){
                return rng->z1 * sigma + mu;
        }

        double u1, u2;
        do{
                u1 = tl_random_double(rng);//  rand() * (1.0 / RAND_MAX);
                u2 = tl_random_double(rng);//rand() * (1.0 / RAND_MAX);
        } while (u1 <= DBL_EPSILON);

        double z0;
        z0 = sqrt(-2.0 * log(u1)) * cos(M_2PI * u2);
        rng->z1 = sqrt(-2.0 * log(u1)) * sin(M_2PI * u2);

        return z0 * sigma + mu;
}

double tl_random_gamma(struct rng_state* rng, double shape, double scale)
{
        return scale * tl_standard_gamma(rng, shape);
}

double tl_standard_gamma(struct rng_state* rng, double shape)
{
        double b, c;
        double U, V, X, Y;

        if (shape == 1.0)
        {
                return tl_standard_exponential(rng);
        }
        else if (shape < 1.0)
        {
                for (;;)
                {
                        U = tl_random_double(rng);
                        V = tl_standard_exponential(rng);
                        if (U <= 1.0 - shape)
                        {
                                X = pow(U, 1./shape);
                                if (X <= V)
                                {
                                        return X;
                                }
                        }
                        else
                        {
                                Y = -log((1-U)/shape);
                                X = pow(1.0 - shape + shape*Y, 1./shape);
                                if (X <= (V + Y))
                                {
                                        return X;
                                }
                        }
                }
        }
        else
        {
                b = shape - 1./3.;
                c = 1./sqrt(9*b);
                for (;;)
                {
                        do
                        {
                                X = tl_gauss(rng);
                                V = 1.0 + c*X;
                        } while (V <= 0.0);

                        V = V*V*V;
                        U = tl_random_double(rng);
                        if (U < 1.0 - 0.0331*(X*X)*(X*X)) return (b*V);
                        if (log(U) < 0.5*X*X + b*(1. - V + log(V))) return (b*V);
                }
        }
}

double tl_standard_exponential(struct rng_state* rng)
{
    /* We use -log(1-U) since U is [0, 1) */
        return -log(1.0 - tl_random_double(rng));
}

double tl_gauss(struct rng_state* rng)
{
        if (rng->has_gauss) {
                const double tmp = rng->gauss;
                rng->gauss = 0;
                rng->has_gauss = 0;
                return tmp;
        } else {
                double f, x1, x2, r2;

                do {
                        x1 = 2.0*tl_random_double(rng) - 1.0;
                        x2 = 2.0*tl_random_double(rng) - 1.0;
                        r2 = x1*x1 + x2*x2;
                }
                while (r2 >= 1.0 || r2 == 0.0);

                /* Box-Muller transform */
                f = sqrt(-2.0*log(r2)/r2);
                /* Keep for next call */
                rng->gauss = f*x1;
                rng->has_gauss = 1;
                return f*x2;
        }
}


struct rng_state* init_rng(uint64_t seed)
{
        struct rng_state* s = NULL;
        uint64_t z;
        uint64_t sanity;

        MMALLOC(s, sizeof(struct rng_state));
        s->gen = 0;
        s->z1 = 0.0f;
        s->gauss = 0.0;
        s->has_gauss = 0;

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

struct rng_state* init_rng_from_rng(struct rng_state* rng)
{
        struct rng_state* s = NULL;
        int i;
        MMALLOC(s, sizeof(struct rng_state));

        for(i = 0; i < 4;i++){
                s->s[i] = rng->s[i];
                s->gen = 0;
                s->z1 = 0.0;
                s->gauss = rng->gauss;
                s->has_gauss = rng->has_gauss;

        }
        jump(rng);
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


/* This is xoshiro256** 1.0, one of our all-purpose, rock-solid
   generators. It has excellent (sub-ns) speed, a state (256 bits) that is
   large enough for any parallel application, and it passes all tests we
   are aware of.

   For generating just floating-point numbers, xoshiro256+ is even faster.

   The state must be seeded so that it is not everywhere zero. If you have
   a 64-bit seed, we suggest to seed a splitmix64 generator and use its
   output to fill s. */

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
