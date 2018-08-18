/*
 * PargeNum.c
 * Uses moduli of the form 2+ and DC-reduce-plus to reduce.
 *
 *  Created on: Mar 4, 2018
 *      Author: Lidiia Zhitova
 */

#include "matrix.h"
#include "parge_num.h"
#include "marge_num.h"

using namespace CNMA;

void CNMA::dc_reduce_plus(mpz_t a, long unsigned int n)
{
    mpz_t t;
    mpz_init(t);
    dc_reduce_minus(a, 2 * n);
    mpz_tdiv_q_2exp(t, a, n); // right shift
    mpz_tdiv_r_2exp(a, a, n); // 
    mpz_sub(a, a, t);
    if (mpz_cmp_ui(a, 0) < 0)
    {
        mpz_t p;
        mpz_init(p);
        mpz_ui_pow_ui(p, 2, n);
        mpz_add(a, a, p);
        mpz_add_ui(a, a, 1);
        mpz_clear(p);
    }
    mpz_clear(t);
}

void CNMA::plusadd(mpz_t res, mpz_t a, mpz_t b, int n)
{
    int szm, szres;
    mpz_t m;
    mpz_init(m);
    get_mod_plus(m, n);
    szm = mpz_sizeinbase(m, 2);
    mpz_add(res, a, b);
    szres = mpz_sizeinbase(res, 2);

    if (szres >= szm)
    {
        mpz_t t;
        mpz_init(t);
        bits(t, res, n, szres);
        bits(res, res, 0, n);
        mpz_sub(res, res, t);
    }

    dc_reduce_plus(res, n);
}

void CNMA::plussub(mpz_t res, mpz_t a, mpz_t b, int n)
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
        mpz_add_ui(res, res, 1);
    }
    dc_reduce_plus(res, n);
}

void CNMA::plusmul(mpz_t res, mpz_t a, mpz_t b, int n)
{
    mpz_t mod;
    mpz_init(mod);
    get_mod_plus(mod, n);
    mpz_mul(res, a, b);
    dc_reduce_plus(res, n);
    if (mpz_cmp(res, mod) >= 0)
    {
        mpz_t low, high;
        mpz_init(low);
        mpz_init(high);
        // double p = pow(2,n);
        int sz = mpz_sizeinbase(res, (long int)2);
        bits(low, res, 0, n);
        bits(high, res, n, sz);
        plussub(res, low, high, n);
    }
}

void CNMA::get_mod_plus(mpz_t m, int n)
{
    mpz_ui_pow_ui(m, 2, n);
    mpz_add_ui(m, m, 1);
}
