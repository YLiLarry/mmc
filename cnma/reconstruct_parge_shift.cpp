#include "matrix.h"
#include "reconstruct_parge_shift.h"
#include <assert.h>
#include <givaro/givtimer.h>
#include <ostream>

using namespace CNMA;
using namespace std;

// Mi: precomputed array
void CNMA::garner_parge_shift(mpz_t a,         // output
                              size_t N,           // size of r[] and m[]
                              const mpz_t r[], // remainders
                              const mpz_t m[], // moduli
                              uint64_t coef)   // coefficient used to generate parge numbers
{
#if DEBUG_CNMA || TIME_CNMA
    gmp_fprintf(stderr, "########## garner_parge ##########\n");
#endif
#if TIME_CNMA
    Givaro::Timer timer;
    timer.clear();
    timer.start();
#endif
#if DEBUG_CNMA
    fprintf(stderr, " - coef: %lu\n", coef);
    for (size_t i = 0; i < N; i++)
    {
        gmp_fprintf(stderr, " - r[%d] = %Zd - m[%d] = %Zd\n", i, r[i], i, m[i]);
#if CHECK_CNMA
        assert(mpz_cmp(r[i], m[i]) < 0);
#endif
    }
#endif
    assert(N > 0);
    if (N == 1)
    {
        mpz_set(a, r[0]);
        return;
    }
    mpz_set_ui(a, 0);
    // initialize array
    mpz_t arr[N];
    mpz_init_set(arr[0], r[0]);
    // initialize temporary vars
    mpz_t t;
    mpz_init(t);
    mpz_t temp;
    mpz_init(temp);
    mpz_t Mi[N];
    for (size_t i = 0; i < N; i++)
    {
        mpz_init(Mi[i]);
    }
    // garner main body
    // starting from line 7 in Eugene's paper
    for (size_t i = 1; i < N; i++)
    {
        mpz_init(arr[i]);
        mpz_set(t, arr[i - 1]); // line 8
        for (int j = i - 2; j >= 0; j--)
        {                          // for loop line 9
            mpz_mul(t, t, m[j]);   // line 10
            mpz_add(t, t, arr[j]); // line 11
        }                          // end for
        mpz_sub(t, r[i], t);       // line 13
        // first term
        mpz_mul_2exp(arr[i], t, (coef << i) - 1);
        // second term
        mpz_mul_2exp(temp, t, coef - 1);
        mpz_sub(arr[i], arr[i], temp);
        // last term

        mpz_add(arr[i], arr[i], t);
        // done
        mpz_mod(arr[i], arr[i], m[i]);
    }                       // end for line 15
    mpz_set(a, arr[N - 1]); // line 16
    for (int i = N - 1; i >= 0; i--)
    { // line 17
        mpz_mul(a, a, m[i]);
        mpz_add(a, a, arr[i]);
        mpz_clear(arr[i]);
    } // line for line 20
    // free used temporary vars
    mpz_clear(t);
    mpz_clear(temp);
    // outputs in a
#if DEBUG_CNMA
    gmp_fprintf(stderr, " - recovered: %Zd\n", a);
#endif
#if TIME_CNMA
    timer.stop();
    cerr << "Timer: " << timer << endl;
#endif
#if DEBUG_CNMA || TIME_CNMA
    gmp_fprintf(stderr, "########## garner_parge ends ##########\n");
#endif
}

// see header for documentation
void CNMA::garner_parge_shift_mixed(mpz_t a,                  
                                    size_t N,                    
                                    const mpz_t r[],          
                                    const uint64_t bitsize[], 
                                    const mpz_t m[],          
                                    const mpz_t Mi[],         
                                    uint64_t coef,            
                                    mpz_t work[])             
{
#if DEBUG_CNMA || TIME_MMA
    gmp_fprintf(stderr, "########## garner_parge_shift_mixed ##########\n");
#endif
#if TIME_CNMA
    Givaro::Timer timer;
    timer.clear();
    timer.start();
#endif
#if DEBUG_CNMA
    for (size_t i = 0; i < N; i++)
    {
        gmp_fprintf(stderr, " - r[%d] = %Zd - m[%d] = %Zd - bitsize[%d] = %d\n", i, r[i], i, m[i], i, bitsize[i]);
#if CHECK_CNMA
        assert(mpz_cmp(r[i], m[i]) < 0);
#endif
    }
#endif
    assert(N > 0);
    if (N == 1)
    {
        mpz_set(a, r[0]);
        return;
    }
    mpz_set_ui(a, 0);
    mpz_set(work[0], r[0]);
    // initialize temporary vars
    mpz_t t;
    mpz_init(t);
    mpz_t temp;
    mpz_init(temp);
    // garner_parge_shift_mixed main body
    // starting from line 7 in Eugene's paper
    for (size_t i = 1; i < N; i++)
    {
        mpz_set(t, work[i - 1]); // line 8
        int j = i - 2;
        while (j >= 3) { 
            // for loop line 9
            // line 10
            mpz_add(temp, t, work[j]);
            mpz_mul_2exp(t, t, bitsize[j]);
            mpz_add(t, t, temp);
            // gmp_fprintf(stderr, "What is temp = %Zd\n", temp);
            j--;
        }   
        // m[2] is a random prime                     
        if (j >= 2) 
        {
            mpz_mul(t, t, m[2]);
            mpz_add(t, t, work[2]);
            j--;
        }
        // m[1] is 2^n+3
        if (j >= 1) {
            mpz_set(temp, work[1]);
            mpz_addmul_ui(temp, t, 3); // temp = 3*t + a[1]
            mpz_mul_2exp(t, t, bitsize[1]);
            mpz_add(t, t, temp);
            j--;
        }
        // m[0] is 2^n
        if (j >= 0) {
            mpz_mul_2exp(t, t, bitsize[0]); 
            mpz_add(t, t, work[0]); // t += t*2^n
        }
        mpz_sub(t, r[i], t);     // line 13
        mpz_mul(temp, t, Mi[i]); // line 14
        mpz_mod(work[i], temp, m[i]);
    }                       // end for line 15
    mpz_set(a, work[N - 1]); // line 16
    int i = N - 1;
    while (i >= 3)
    { // line 17
        mpz_add(temp, a, work[i]);
        mpz_mul_2exp(a, a, bitsize[i]);
        mpz_add(a, a, temp);
        i--;
    } // line for line 20
    // m[2] is a random prime
    if (i >= 2) {
        mpz_mul(a, a, m[2]);
        mpz_add(a, a, work[2]);
    }
    // m[1] is 2^n+3
    if (i >= 1) {
        mpz_set(temp, work[1]);
        mpz_addmul_ui(temp, a, 3); 
        // temp = 3*a + a[1]
        mpz_mul_2exp(a, a, bitsize[1]);
        mpz_add(a, a, temp); 
    }
    // m[0] is 2^n
    if (i >= 0) {
        mpz_mul_2exp(a, a, bitsize[0]);
        mpz_add(a, a, work[0]);
    }
    // free used temporary vars
    mpz_clear(t);
    mpz_clear(temp);
    // outputs in a
#if DEBUG_CNMA
    gmp_fprintf(stderr, " - a: %Zd\n", a);
#endif
#if TIME_CNMA
    timer.stop();
    cerr << "Timer: " << timer << endl;
#endif
#if DEBUG_CNMA || TIME_CNMA
    gmp_fprintf(stderr, "########## garner_parge_shift_mixed ends ##########\n");
#endif
}



// void prod(mpz_t M, mpz_t Mi, mpz_t *m, int N)
// {
//     mpz_set_ui(M, 1);
//     mpz_t plchldr;
//     mpz_init(plchldr);
//     for (int i = 1; i <= N - 1; i++)
//     {
//         mpz_mul(M, M, m[i - 1]);
//         mpz_gcdext(plchldr, Mi, plchldr, M, m[i]);
//     }
//     mpz_clear(plchldr);
// }
