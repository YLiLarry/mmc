/*
 * MargeNum.c
 * Uses moduli of the form 2- and DC-reduce-minus to reduce
 *
 *  Created on: Mar 4, 2018
 *      Author: Lidiia Zhitova
 */

#include "matrix.h"
#include "marge_num.h"

using namespace CNMA;

void CNMA::dc_reduce_minus(mpz_t a, unsigned long int n)
{
    mpz_t t, p;
    int s;
    int b;
    uint64_t k;
    mpz_init(t);
    mpz_init_set_ui(p, 1);
    mpz_mul_2exp(p, p, n);
    int cmp_val = mpz_cmp(a, p);

    while (cmp_val >= 0)
    {
        s = mpz_sizeinbase(a, 2);
        b = (s - 1) / n + 1;
        b >>= 1;
        k = b * n;
        mpz_tdiv_q_2exp(t, a, k);
        mpz_tdiv_r_2exp(a, a, k);
        mpz_add(a, a, t);
        cmp_val = mpz_cmp(a, p);
    }
    mpz_clear(t);
    mpz_clear(p);
}

void CNMA::minadd(mpz_t res, mpz_t a, mpz_t b, int n)
{
    int szm, szres;
    mpz_t m;
    mpz_init(m);
    get_mod(m, n);
    szm = mpz_sizeinbase(m, 2);
    mpz_add(res, a, b);
    szres = mpz_sizeinbase(res, 2);

    if (szres >= szm)
    {
        mpz_t t;
        mpz_init(t);
        bits(t, res, n, szres);
        bits(res, res, 0, n);
        mpz_add(res, res, t);
    }
    dc_reduce_minus(res, n);
}

void CNMA::minsub(mpz_t res, mpz_t a, mpz_t b, int n)
{
    int sza, szb;
    sza = mpz_sizeinbase(a, 2);
    szb = mpz_sizeinbase(b, 2);
    mpz_sub(res, a, b);

    if (sza < szb)
    {
        mpz_t t;
        mpz_init(t);
        bits(t, res, n, szb);
        bits(res, res, 0, n);
        mpz_sub_ui(res, res, 1);
    }
    dc_reduce_minus(res, n);
}

void CNMA::minmul(mpz_t res, mpz_t a, mpz_t b, int n)
{
    mpz_t mod;
    mpz_init(mod);
    get_mod(mod, n);
    mpz_mul(res, a, b);
    dc_reduce_minus(res, n);
    if (mpz_cmp(res, mod) >= 0)
    {
        mpz_t low, high;
        mpz_init(low);
        mpz_init(high);
        int sz = mpz_sizeinbase(res, 2);
        bits(low, res, 0, n);
        bits(high, res, n, sz);
        minadd(res, low, high, n);
    }
}

void CNMA::get_mod(mpz_t m, int n)
{
    mpz_ui_pow_ui(m, 2, n);
    mpz_sub_ui(m, m, 1);
}
