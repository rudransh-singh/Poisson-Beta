/*% PoissonBetaCparam.c
*
* By Jong Kyoung Kim 
* jkkim@ebi.ac.uk
* Last Update: July 10 2012
*/

/* 
   A C-program for MT19937-64 (2004/9/29 version).
   Coded by Takuji Nishimura and Makoto Matsumoto.

   This is a 64-bit version of Mersenne Twister pseudorandom number
   generator.

   Before using, initialize the state by using init_genrand64(seed)  
   or init_by_array64(init_key, key_length).

   Copyright (C) 2004, Makoto Matsumoto and Takuji Nishimura,
   All rights reserved.                          

   Redistribution and use in source and binary forms, with or without
   modification, are permitted provided that the following conditions
   are met:

     1. Redistributions of source code must retain the above copyright
        notice, this list of conditions and the following disclaimer.

     2. Redistributions in binary form must reproduce the above copyright
        notice, this list of conditions and the following disclaimer in the
        documentation and/or other materials provided with the distribution.

     3. The names of its contributors may not be used to endorse or promote 
        products derived from this software without specific prior written 
        permission.

   THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
   "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
   LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR
   A PARTICULAR PURPOSE ARE DISCLAIMED.  IN NO EVENT SHALL THE COPYRIGHT OWNER OR
   CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
   EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
   PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
   PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
   LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
   NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
   SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

   References:
   T. Nishimura, ``Tables of 64-bit Mersenne Twisters''
     ACM Transactions on Modeling and 
     Computer Simulation 10. (2000) 348--357.
   M. Matsumoto and T. Nishimura,
     ``Mersenne Twister: a 623-dimensionally equidistributed
       uniform pseudorandom number generator''
     ACM Transactions on Modeling and 
     Computer Simulation 8. (Jan. 1998) 3--30.

   Any feedback is very welcome.
   http://www.math.hiroshima-u.ac.jp/~m-mat/MT/emt.html
   email: m-mat @ math.sci.hiroshima-u.ac.jp (remove spaces)
*/

#ifndef _Cparam
#define _Cparam
#include <stdlib.h>
#include <time.h>
#include <math.h>

#define NN 312
#define MM 156
#define MATRIX_A 0xB5026F5AA96619E9ULL
#define UM 0xFFFFFFFF80000000ULL /* Most significant 33 bits */
#define LM 0x7FFFFFFFULL /* Least significant 31 bits */
#define M_lnSqrt2PI 0.91893853320467274178


/* The array for the state vector */
static unsigned long long mt[NN]; 
/* mti==NN+1 means mt[NN] is not initialized */
static int mti=NN+1; 

typedef struct _tagCparam {
    int num_gene, num_replicate;
    double *Koni, *Koffi, *Si;
    double **Pij;
} Cparam;

/* initializes mt[NN] with a seed */
void init_genrand64(unsigned long long seed)
{
    mt[0] = seed;
    for (mti=1; mti<NN; mti++) 
        mt[mti] =  (6364136223846793005ULL * (mt[mti-1] ^ (mt[mti-1] >> 62)) + mti);
}

/* generates a random number on [0, 2^64-1]-interval */
unsigned long long genrand64_int64(void)
{
    int i;
    unsigned long long x;
    static unsigned long long mag01[2]={0ULL, MATRIX_A};

    if (mti >= NN) { /* generate NN words at one time */

        /* if init_genrand64() has not been called, */
        /* a default initial seed is used     */
        if (mti == NN+1) 
            init_genrand64(5489ULL); 

        for (i=0;i<NN-MM;i++) {
            x = (mt[i]&UM)|(mt[i+1]&LM);
            mt[i] = mt[i+MM] ^ (x>>1) ^ mag01[(int)(x&1ULL)];
        }
        for (;i<NN-1;i++) {
            x = (mt[i]&UM)|(mt[i+1]&LM);
            mt[i] = mt[i+(MM-NN)] ^ (x>>1) ^ mag01[(int)(x&1ULL)];
        }
        x = (mt[NN-1]&UM)|(mt[0]&LM);
        mt[NN-1] = mt[MM-1] ^ (x>>1) ^ mag01[(int)(x&1ULL)];

        mti = 0;
    }
  
    x = mt[mti++];

    x ^= (x >> 29) & 0x5555555555555555ULL;
    x ^= (x << 17) & 0x71D67FFFEDA60000ULL;
    x ^= (x << 37) & 0xFFF7EEE000000000ULL;
    x ^= (x >> 43);

    return x;
}

double loggamma(double x) {
    static double gamma_series[] = {
        76.18009172947146,
        -86.50532032941677,
        24.01409824083091,
        -1.231739572450155,
        0.1208650973866179e-2,
        -0.5395239384953e-5
    };
    int i;
    double denom, x1, series;

    /* Lanczos method */
    denom = x+1;
    x1 = x + 5.5;
    series = 1.000000000190015;
    for(i = 0; i < 6; i++) {
        series += gamma_series[i] / denom;
        denom += 1.0;
    }
    return( M_lnSqrt2PI + (x+0.5)*log(x1) - x1 + log(series/x) );
}

/* generates a random number on (0,1)-real-interval */
double my_rand(void)
{
    return ((genrand64_int64() >> 12) + 0.5) * (1.0/4503599627370496.0);
    /* divided by 2^32 */
}

double logpbeta(double Kon, double Koff, double count_data, double size_factor, double S, double x) {
    double log_Px;    
    
    log_Px = (Koff-1.0)*log(x) + (Kon-1.0)*log(1.0-x) + log(1.0-x)*count_data + size_factor*(S*(x-1.0));
    return log_Px;
}

double logpgammakoff(double alpha_koff, double beta_koff, double Kon, double *lnPij, int num_gene, int num_replicate, int ii, double x) {
    double log_Px, sum_Pij=0.0;
    int jj;
    
    for (jj=0; jj < num_replicate; jj++) {
        sum_Pij = sum_Pij + lnPij[jj*num_gene+ii];
    }     
    log_Px = -1.0*x/beta_koff + (alpha_koff-1.0)*log(x) + (double)(num_replicate*(loggamma(x+Kon) - loggamma(x))) + (x-1.0)*sum_Pij;
    return log_Px;
}

double logpgammakon(double alpha_kon, double beta_kon, double Koff, double *lnoneminusPij, int num_gene, int num_replicate, int ii, double x) {
    double log_Px, sum_Pij=0.0;
    int jj;
    
    for (jj=0; jj < num_replicate; jj++) {
        sum_Pij = sum_Pij + lnoneminusPij[jj*num_gene+ii];
    }    
    log_Px = -1.0*x/beta_kon + (alpha_kon-1.0)*log(x) + (double)(num_replicate*(loggamma(x+Koff) - loggamma(x))) + (x-1.0)*sum_Pij;
    return log_Px;
}

double logpgammas(double alpha_s, double beta_s, double *size_factor, double *Pij, double sum_count, int num_gene, int num_replicate, int ii, double x) {
    double log_Px, sum_Pij=0.0;
    int jj;
    
    for (jj=0; jj < num_replicate; jj++) {
        sum_Pij = sum_Pij + size_factor[jj*num_gene+ii]*(Pij[jj*num_gene+ii]-1.0);
    }
    log_Px = -1.0*x/beta_s + (alpha_s-1.0)*log(x) + sum_count*log(x) + x*sum_Pij;
    return log_Px;
}
        
#endif