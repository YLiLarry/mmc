#include "matrix.h"
#include "reconstruct.h"
#include <assert.h>

void precompute_Mi(mpz_t Mi[], const mpz_t m[], const size_t N)
{
#if DEBUG_MMC || TIME_MMC
    gmp_fprintf(stderr, ".......... precompute_Mi ..........\n");
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
#if DEBUG_MMC
        gmp_fprintf(stderr, "precompute_Mi iteration#%d: M=%Zd, Mi[%d]=%Zd\n", i, M, i, Mi[i]);
#endif
#if DEBUG_MMC
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
#if DEBUG_MMC || TIME_MMC
    gmp_fprintf(stderr, ".......... precompute_Mi ends ..........\n");
#endif
}

// Mi: precomputed array
void garner(mpz_t a, int N, const mpz_t r[], const mpz_t m[], const mpz_t Mi[])
{
#if DEBUG_MMC
    gmp_fprintf(stderr, "########## garner ##########\n");
#endif
#if DEBUG_MMC
    for (size_t i = 0; i < N; i++)
    {
        gmp_fprintf(stderr, " - r[%d] = %Zd - m[%d] = %Zd\n", i, r[i], i, m[i]);
#if TEST_MMC
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
    for (int i = 0; i < N; i++)
    {
        mpz_init(arr[i]);
    }
    mpz_set(arr[0], r[0]);
    // initialize temporary vars
    mpz_t t;
    mpz_init(t);
    mpz_t temp;
    mpz_init(temp);
    // garner main body
    // starting from line 7 in Eugene's paper
    for (int i = 1; i < N; i++)
    {
        mpz_set(t, arr[i - 1]); // line 8
        for (int j = i - 2; j >= 0; j--)
        {                          // for loop line 9
            mpz_mul(t, t, m[j]);   // line 10
            mpz_add(t, t, arr[j]); // line 11
        }                          // end for
        mpz_sub(t, r[i], t);       // line 13
        mpz_mod(temp, Mi[i], m[i]);
        mpz_mul(arr[i], t, temp); // line 14
    }                             // end for line 15
    mpz_set(a, arr[N - 1]);       // line 16
    for (int i = N - 1; i >= 0; i--)
    { // line 17
        mpz_mul(a, a, m[i]);
        mpz_add(a, a, arr[i]);
    } // line for line 20
    // free used temporary vars
    mpz_clear(t);
    mpz_clear(temp);
    for (int i = 0; i < N; i++)
    {
        mpz_clear(arr[i]);
    }
    // outputs in a
#if DEBUG_MMC
    gmp_fprintf(stderr, " - a: %Zd\n", a);
#endif
#if DEBUG_MMC
    gmp_fprintf(stderr, "########## garner ends ##########\n");
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
