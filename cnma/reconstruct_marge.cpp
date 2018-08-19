#include "matrix.h"
#include "reconstruct_marge.h"
#include "marge_num.h"
#include <assert.h>
#include <givaro/givtimer.h>

using namespace CNMA;
using namespace std;

void CNMA::precompute_Mi_marge(mpz_t Mi[], const mpz_t m[], const size_t N)
{
#if DEBUG_CNMA || TIME_CNMA
    gmp_fprintf(stderr, ".......... precompute_Mi_marge ..........\n");
#endif
    // line 1
    mpz_t M;
    mpz_init(M);
    mpz_set_ui(M, 1);
    // line 2
    for (size_t i = 1; i < N; i++)
    {
        mpz_mul(M, M, m[i - 1]);    // line 3
        mpz_invert(Mi[i], M, m[i]); // line 4
#if DEBUG_CNMA
        gmp_fprintf(stderr, "precompute_Mi_marge iteration#%d: M=%Zd, Mi[%d]=%Zd\n", i, M, i, Mi[i]);
#endif
#if DEBUG_CNMA
        gmp_fprintf(stderr, " - M: %Zd\n", M);
        gmp_fprintf(stderr, " - Mi: [");
        gmp_fprintf(stderr, " %Zd ", Mi[1]);
        for (size_t i = 2; i < N; i++)
        {
            gmp_fprintf(stderr, ", %Zd ", Mi[i]);
        }
        gmp_fprintf(stderr, "]\n");
#endif
    }
    mpz_clear(M);
    // outputs Mi[]
#if DEBUG_CNMA || TIME_CNMA
    gmp_fprintf(stderr, ".......... precompute_Mi_marge ends ..........\n");
#endif
}

void CNMA::garner_marge(mpz_t a,               // output
                        int N,                 // size of r[], expo[], m[], Mi[], work[]
                        const mpz_t r[],       // remainders
                        const uint64_t expo[], // array of exponents n as in moduli 2^n - 1
                        const mpz_t m[],       // moduli of the form 2^n - 1
                        const mpz_t Mi[],      // precomputed Mi (see paper)
                        mpz_t work[])          // a work array, caller is responsible for initializing and freeing this for efficiency reason
{
#if DEBUG_CNMA || TIME_CNMA
    gmp_fprintf(stderr, "########## garner_marge ##########\n");
#endif
#if TIME_CNMA
    Givaro::Timer timer;
    timer.clear();
    timer.start();
#endif
#if DEBUG_CNMA
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
    mpz_set(work[0], r[0]);
    // initialize temporary vars
    mpz_t t;
    mpz_init(t);
    mpz_t temp;
    mpz_init(temp);
    // garner_marge main body
    // starting from line 7 in Eugene's paper
    for (int i = 1; i < N; i++)
    {
        mpz_set(t, work[i - 1]); // line 8
        for (int j = i - 2; j >= 0; j--)
        { // for loop line 9
            // mpz_mul(t, t, m[j]);
            // line 10
            mpz_sub(temp, work[j], t);
            mpz_mul_2exp(t, t, expo[j]);
            mpz_add(t, t, temp); // line 11
            // gmp_fprintf(stderr, "What is temp = %Zd\n", temp);
        }                        // end for
        mpz_sub(t, r[i], t);     // line 13
        mpz_mul(work[i], t, Mi[i]); // line 14
        // mpz_mod(work[i], temp, m[i]);
        dc_reduce_minus(work[i], expo[i]);
    }                        // end for line 15
    mpz_set(a, work[N - 1]); // line 16
    for (int i = N - 1; i >= 0; i--)
    { // line 17
        mpz_sub(temp, work[i], a);
        mpz_mul_2exp(a, a, expo[i]);
        mpz_add(a, a, temp);
    } // line for line 20
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
    gmp_fprintf(stderr, "########## garner_marge ends ##########\n");
#endif
}

// Mi: precomputed array
void garner_simple_marge(mpz_t a, int N, const mpz_t r[], const mpz_t m[], const mpz_t Mi[])
{
#if DEBUG_CNMA
    gmp_fprintf(stderr, "########## garner_simple_marge ##########\n");
#endif
#if DEBUG_CNMA
    for (size_t i = 0; i < N; i++)
    {
        gmp_fprintf(stderr, " - r[%d] = %Zd - m[%d] = %Zd\n", i, r[i], i, m[i]);
#if CHECK_CNMA
        assert(mpz_cmp(r[i], m[i]) < 0);
#endif
    }
#endif
    assert(N > 0);
    mpz_set(a, r[0]); // line 6
    if (N == 1)
    {
        return;
    }
    // line 7
    mpz_t M;
    mpz_init_set(M, m[0]);
    mpz_t t;
    mpz_init(t);
    // line 8
    for (int i = 1; i < N; i++)
    {
        mpz_sub(t, r[i], a);  // line 8
        mpz_mul(t, t, Mi[i]); // line 9
        mpz_mod(t, t, m[i]);  // line 10
        mpz_addmul(a, t, M);  // line 11
        mpz_mul(M, M, m[i]);  // line 12
    }                         // end for line 15
    mpz_clear(t);
    mpz_clear(M);
    // outputs in a
#if DEBUG_CNMA
    gmp_fprintf(stderr, " - a: %Zd\n", a);
#endif
#if DEBUG_CNMA
    gmp_fprintf(stderr, "########## garner_simple_marge ends ##########\n");
#endif
}

void prod(mpz_t M, mpz_t Mi, mpz_t *m, int N)
{
    mpz_set_ui(M, 1);
    mpz_t plchldr;
    mpz_init(plchldr);
    for (int i = 1; i <= N - 1; i++)
    {
        mpz_mul(M, M, m[i - 1]);
        mpz_gcdext(plchldr, Mi, plchldr, M, m[i]);
    }
    mpz_clear(plchldr);
}
