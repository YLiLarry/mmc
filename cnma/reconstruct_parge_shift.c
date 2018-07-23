#include "matrix.h"
#include "reconstruct_marge.h"
#include <assert.h>

// Mi: precomputed array
void garner_parge_shift(mpz_t a,         // output
                        int N,           // size of r[] and m[]
                        const mpz_t r[], // remainders
                        const mpz_t m[], // moduli
                        uint64_t coef)   // coefficient used to generate parge numbers
{
#if DEBUG_MMC
    gmp_fprintf(stderr, "########## garner_parge ##########\n");
#endif
#if DEBUG_MMC
    fprintf(stderr, " - coef: %lu\n", coef);
    for (size_t i = 0; i < N; i++)
    {
        gmp_fprintf(stderr, " - r[%d] = %Zd - m[%d] = %Zd\n", i, r[i], i, m[i]);
#if CHECK_MMC
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
    for (int i = 1; i < N; i++)
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
#if DEBUG_MMC
    gmp_fprintf(stderr, " - recovered: %Zd\n", a);
#endif
#if DEBUG_MMC
    gmp_fprintf(stderr, "########## garner_parge ends ##########\n");
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
