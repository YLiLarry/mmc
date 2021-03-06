/* -*- mode: C++; tab-width: 8; indent-tabs-mode: t; c-basic-offset: 8 -*- */
// vim:sts=8:sw=8:ts=8:noet:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s

/*
 * Copyright (C) FFLAS-FFPACK
 * Written by Pascal Giorgi <pascal.giorgi@lirmm.fr>
 *
 * This file is Free Software and part of FFLAS-FFPACK.
 *
 * ========LICENCE========
 * This file is part of the library FFLAS-FFPACK.
 *
 * FFLAS-FFPACK is free software: you can redistribute it and/or modify
 * it under the terms of the  GNU Lesser General Public
 * License as published by the Free Software Foundation; either
 * version 2.1 of the License, or (at your option) any later version.
 *
 * This library is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 * Lesser General Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public
 * License along with this library; if not, write to the Free Software
 * Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301  USA
 * ========LICENCE========
 *.
 */

#define __FFLASFFPACK_SEQUENTIAL

#include "fflas-ffpack/fflas-ffpack-config.h"
#include <chrono>
#include <iostream>
#include <string>
#include <vector>
using namespace std;

#include "fflas-ffpack/fflas/fflas.h"
#include "fflas-ffpack/field/rns-double.h"
#include "fflas-ffpack/field/rns-integer.h"
#include "fflas-ffpack/utils/args-parser.h"
#include "fflas-ffpack/utils/timer.h"
#include "givaro/givinteger.h"
#include <gmp++/gmp++.h>
#include "givaro/zring.h"
#include <givaro/givrns.h>

#define __
{
#include "flint/flint.h"
#include "flint/fmpz.h"
#include "flint/fmpz_mat.h"
#include "flint/long_extras.h"
#include "flint/longlong.h"
}

#include "algebramix/crt_integer.hpp"
#include "algebramix/matrix.hpp"
#include "algebramix/matrix_crt.hpp"
#include "algebramix/matrix_integer.hpp"

template <typename T>
void write_matrix(Givaro::Integer p, size_t m, size_t n, T *C, size_t ldc)
{

    size_t www = (p.bitsize() * log(2.)) / log(10.);
    for (size_t i = 0; i < m; ++i)
    {
        cout << "[ ";
        cout.width(www + 1);
        cout << std::right << C[i * ldc];
        for (size_t j = 1; j < n; ++j)
        {
            cout << " ";
            cout.width(www);
            cout << std::right << C[i * ldc + j];
        }
        cout << "]" << endl;
    }
    cout << endl;
}

template <typename Field>
void run_bench(size_t n, size_t primes_bits, size_t b, size_t iters, size_t seed, bool matmul = false)
{
    Givaro::Integer p;
    FFLAS::Timer chrono;
    double timeFFLASToRNS = 0., timeFFLASFromRNS = 0., timeFFLASPrecomp = 0.;
    double timeFlintToRNS = 0., timeFlintFromRNS = 0., timeFlintPrecomp = 0.;
    double timeMMXToRNS = 0., timeMMXFromRNS = 0., timeMMXPrecomp = 0.;
    double timeMMXToRNS2 = 0., timeMMXFromRNS2 = 0., timeMMXPrecomp2 = 0.;

    size_t loop = 0;
    size_t bits = b;
    auto start = std::chrono::system_clock::now();
    auto end = std::chrono::system_clock::now();
    for (; (loop < iters || std::chrono::duration<double>(end - start) < std::chrono::seconds(10)); loop++)
    {
        //cout<<loop<<" "<<end-start<<endl;
        Givaro::Integer::random_exact_2exp(p, bits);
        //nextprime(p,p);
        //cout<<p<<endl;
        Field F(p);

        typename Field::RandIter Rand(F, seed);
        typename Field::Element_ptr A;
        A = FFLAS::fflas_new(F, n, n);

        for (size_t i = 0; i < n; ++i)
            for (size_t j = 0; j < n; ++j)
                Rand.random(A[i * n + j]);

        // A is a random array of size n^2 in field F(p)

        if (matmul)
            bits = 2 * b;

        // FLINT RNS CODE //
        {
            fmpz_mat_t AA;
            fmpz_mat_init(AA, n, n);
            for (size_t i = 0; i < n; ++i)
                for (size_t j = 0; j < n; ++j)
                    fmpz_set_mpz(fmpz_mat_entry(AA, i, j), *(reinterpret_cast<const mpz_t *>(A + i * n + j)));

            fmpz_comb_t comb;
            fmpz_comb_temp_t comb_temp;
            nmod_mat_t *mod_A;
            mp_limb_t *primes;
            mp_limb_t *residues;
            size_t num_primes = (bits + primes_bits - 1) / primes_bits;

            /* FLINT RNS initialization */
            residues = (mp_limb_t *)flint_malloc(sizeof(mp_limb_t) * num_primes);
            mod_A = (nmod_mat_t *)flint_malloc(sizeof(nmod_mat_t) * num_primes);
            primes = (mp_limb_t *)flint_malloc(sizeof(mp_limb_t) * num_primes);
            primes[0] = n_nextprime(UWORD(1) << primes_bits, 0);
            nmod_mat_init(mod_A[0], AA->r, AA->c, primes[0]);
            for (size_t i = 1; i < num_primes; i++)
            {
                primes[i] = n_nextprime(primes[i - 1], 0);
                nmod_mat_init(mod_A[i], AA->r, AA->c, primes[i]);
            }

            chrono.clear();
            chrono.start();
            fmpz_comb_init(comb, primes, num_primes);
            fmpz_comb_temp_init(comb_temp, comb);
            chrono.stop();
            timeFlintPrecomp += chrono.usertime();
            chrono.clear();
            chrono.start();

            /* Calculate residues of AA */
            for (long i = 0; i < AA->r; i++)
            {
                for (long j = 0; j < AA->c; j++)
                {
                    fmpz_multi_mod_ui(residues, &(AA->rows[i][j]), comb, comb_temp);
                    for (size_t k = 0; k < num_primes; k++)
                        mod_A[k]->rows[i][j] = residues[k];
                }
            }

            chrono.stop();
            timeFlintToRNS += chrono.usertime();
            chrono.clear();
            chrono.start();

            /* Chinese remaindering */
            for (long i = 0; i < AA->r; i++)
            {
                for (long j = 0; j < AA->c; j++)
                {
                    for (size_t k = 0; k < num_primes; k++)
                        residues[k] = mod_A[k]->rows[i][j];
                    fmpz_multi_CRT_ui(&AA->rows[i][j], residues, comb, comb_temp, 1);
                }
            }
            chrono.stop();
            timeFlintFromRNS += chrono.usertime();

            fmpz_mat_clear(AA);
            flint_free(mod_A);
            fmpz_comb_temp_clear(comb_temp);
            fmpz_comb_clear(comb);
            flint_free(residues);
            flint_free(primes);
        }
        //END FLINT CODE //

        // MMX Code
        {
            //using namespace mmx;
            typedef mmx::integer mT;
            typedef mmx::matrix_integer mM;
            mmx::threads_set_number(1);
            mmx::matrix<mT, mM> AA(1, n, n);
            mmx::matrix<mT, mM> BB(1, n, n);
            for (size_t i = 0; i < n; ++i)
                for (size_t j = 0; j < n; ++j)
                    mpz_set(*AA(i, j), *(reinterpret_cast<const mpz_t *>((A + i * n + j))));

            /* compute RNS structure*/

            const mmx::nat mbit = 20; //8 * sizeof(modulus_base) - 12;
            struct mat_crt_naive
            {
                typedef mmx::integer base;
                typedef int modulus_base;
                typedef mmx::modulus_int_preinverse<mbit> modulus_base_variant;
                typedef mmx::Modulus_variant(mmx::integer) modulus_variant;
            };
            struct mat_crt_dicho
            {
                typedef mmx::integer base;
                typedef mmx::integer modulus_base;
                typedef mmx::modulus_integer_naive modulus_base_variant;
                typedef mmx::Modulus_variant(mmx::integer) modulus_variant;
            };

            typedef mmx::crt_naive_transformer<mmx::integer, mat_crt_naive> CRTer;
            //typedef mmx::crt_dicho_transformer<mmx::integer, mat_crt_dicho> CRTer;

            typedef typename CRTer::modulus_base I;
            typedef mmx::modulus<I, typename CRTer::modulus_base_variant> Modulus;
            typedef mmx::modular<Modulus, mmx::modular_matrix_crt<mmx::integer>> Modular;
            typedef mmx::matrix<Modular> Matrix_modular;
            mmx::vector<Modulus> primes;
            chrono.clear();
            chrono.start();
            typedef mmx::moduli_helper<mmx::integer, Modulus,
                                       mmx::pr_prime_sequence_int<mbit>>
                Sequence;

            if (!Sequence::covering(primes, bits))
            {
                std::cerr << "MMX prime for slow CRT error" << std::endl;
                return;
            }

            CRTer crter(primes);
            mmx::nat nn = N(crter);
            chrono.stop();
            timeMMXPrecomp += chrono.usertime();

            // do not count memory allocation
            Matrix_modular *mm1 = mmx::mmx_new<Matrix_modular>(nn);
            for (mmx::nat k = 0; k < nn; k++)
                mm1[k] = Matrix_modular(Modular(), n, n);
            I *aux = mmx::mmx_new<I>(nn);

            chrono.clear();
            chrono.start();

            /* Calculate residues of A */
            for (mmx::nat i = 0; i < n; i++)
                for (mmx::nat j = 0; j < n; j++)
                {
                    mmx::direct_crt(aux, AA(i, j), crter);
                    for (mmx::nat k = 0; k < nn; k++)
                        mm1[k](i, j) = aux[k];
                }

            chrono.stop();
            timeMMXToRNS += chrono.usertime();
            chrono.clear();
            chrono.start();
            /* Chinese remaindering */
            for (mmx::nat i = 0; i < n; i++)
                for (mmx::nat j = 0; j < n; j++)
                {
                    for (mmx::nat k = 0; k < nn; k++)
                        aux[k] = *mm1[k](i, j);
                    mmx::inverse_crt(BB(i, j), aux, crter);
                }
            chrono.stop();
            timeMMXFromRNS += chrono.usertime();
            mmx::mmx_delete<I>(aux, nn);
            mmx::mmx_delete<Matrix_modular>(mm1, nn);
        }

        // MMX Code
        {

            //using namespace mmx;
            typedef mmx::integer mT;
            typedef mmx::matrix_integer mM;
            mmx::threads_set_number(1);
            mmx::matrix<mT, mM> AA(1, n, n);
            mmx::matrix<mT, mM> BB(1, n, n);
            for (size_t i = 0; i < n; ++i)
                for (size_t j = 0; j < n; ++j)
                    mpz_set(*AA(i, j), *(reinterpret_cast<const mpz_t *>((A + i * n + j))));

            /* compute RNS structure*/

            const mmx::nat mbit = 20; //8 * sizeof(modulus_base) - 12;
            struct mat_crt_naive
            {
                typedef mmx::integer base;
                typedef int modulus_base;
                typedef mmx::modulus_int_preinverse<mbit> modulus_base_variant;
                typedef mmx::Modulus_variant(mmx::integer) modulus_variant;
            };
            struct mat_crt_dicho
            {
                typedef mmx::integer base;
                typedef mmx::integer modulus_base;
                typedef mmx::modulus_integer_naive modulus_base_variant;
                typedef mmx::Modulus_variant(mmx::integer) modulus_variant;
            };

            //typedef mmx::crt_naive_transformer<mmx::integer, mat_crt_naive> CRTer;
            typedef mmx::crt_dicho_transformer<mmx::integer, mat_crt_naive> CRTer;

            typedef typename CRTer::modulus_base I;
            typedef mmx::modulus<I, typename CRTer::modulus_base_variant> Modulus;
            typedef mmx::modular<Modulus, mmx::modular_matrix_crt<mmx::integer>> Modular;
            typedef mmx::matrix<Modular> Matrix_modular;
            mmx::vector<Modulus> primes;
            chrono.clear();
            chrono.start();

            typedef mmx::moduli_helper<mmx::integer, Modulus,
                                       mmx::pr_prime_sequence_int<mbit>>
                Sequence;

            if (!Sequence::covering(primes, bits))
            {
                std::cerr << "MMX prime for slow CRT error" << std::endl;
                return;
            }

            CRTer crter(primes);
            mmx::nat nn = N(crter);

            chrono.stop();
            timeMMXPrecomp2 += chrono.usertime();

            // do not count memory allocation
            Matrix_modular *mm1 = mmx::mmx_new<Matrix_modular>(nn);
            for (mmx::nat k = 0; k < nn; k++)
                mm1[k] = Matrix_modular(Modular(), n, n);
            I *aux = mmx::mmx_new<I>(nn);

            chrono.clear();
            chrono.start();

            /* Calculate residues of A */
            for (mmx::nat i = 0; i < n; i++)
                for (mmx::nat j = 0; j < n; j++)
                {
                    mmx::direct_crt(aux, AA(i, j), crter);
                    for (mmx::nat k = 0; k < nn; k++)
                        mm1[k](i, j) = aux[k];
                }

            chrono.stop();
            timeMMXToRNS2 += chrono.usertime();
            chrono.clear();
            chrono.start();
            /* Chinese remaindering */
            for (mmx::nat i = 0; i < n; i++)
                for (mmx::nat j = 0; j < n; j++)
                {
                    for (mmx::nat k = 0; k < nn; k++)
                        aux[k] = *mm1[k](i, j);
                    mmx::inverse_crt(BB(i, j), aux, crter);
                }
            chrono.stop();
            timeMMXFromRNS2 += chrono.usertime();
            mmx::mmx_delete<I>(aux, nn);
            mmx::mmx_delete<Matrix_modular>(mm1, nn);
        }

        // FFLAS RNS CODE
        {
            chrono.clear();
            chrono.start();
            // RNS: contains an array of primes whose product is >= p, each of primes_bits long
            FFPACK::rns_double RNS(p, primes_bits);
            // RNSInteger is a decorator (wrapper) on FFPACK::rns_double that adds some functionalities
            typedef FFPACK::RNSInteger<FFPACK::rns_double> RnsDomain;
            RnsDomain Zrns(RNS);
            typename RnsDomain::Element_ptr mod_A = FFLAS::fflas_new(Zrns, n, n); // allocate array of n^2 * RNS.size() doubles
            chrono.stop();
            timeFFLASPrecomp += chrono.usertime();
            chrono.clear();
            chrono.start();
            // from integer to rns
            // A: integer array
            // mod_A: output
            // k: ceil(|p|/16) how many 16-bits chunks
            // lda: =n? what is lda?
            FFLAS::finit_rns(Zrns, n, n, (p.bitsize() / 16) + ((p.bitsize() % 16) ? 1 : 0), A, n, mod_A);
            chrono.stop();
            timeFFLASToRNS += chrono.usertime();
            chrono.clear();
            chrono.start();
            // Zrns: primes
            // A: random matrix
            // n:
            // mod_A: n^2 dimensional matrix
            // from rns to integer
            FFLAS::fconvert_rns(Zrns, n, n, Givaro::Integer(0), A, n, mod_A);
            chrono.stop();
            timeFFLASFromRNS += chrono.usertime();
            FFLAS::fflas_delete(mod_A);
        }
        end = std::chrono::system_clock::now();
    }
#define SPC1 11
#define SPC2 12
#define SPC3 30
#define PREC 2

    cout << setw(SPC1 + 5)
         << loop << " |"
         << setw(SPC1)
         << n << " |"
         << setw(SPC1)
         << bits << " |"
         << setw(SPC1)
         << primes_bits << " |";

    //#define TIMING(x,y,z) cout<<setw(SPC2)<<std::fixed << std::setprecision(PREC)<<(x+y+z)/loop<<setw(SPC3)<<"("+to_string(x/loop)+","+to_string(y/loop)+","+to_string(z/loop)+")"<<" |";
#define TIMING(x, y, z) cout << setw(SPC2) << std::fixed << std::setprecision(PREC)          \
                             << setw(SPC3 / 3) << ((x / loop) * 1000000.) << ","             \
                             << setw(SPC3 / 3) << ((y / loop) * (1000000. / (n * n))) << "," \
                             << setw(SPC3 / 3) << ((z / loop) * (1000000. / (n * n))) << " |";
    TIMING(timeFFLASPrecomp, timeFFLASFromRNS, timeFFLASToRNS);
    TIMING(timeFlintPrecomp, timeFlintFromRNS, timeFlintToRNS);
    TIMING(timeMMXPrecomp, timeMMXFromRNS, timeMMXToRNS);
    TIMING(timeMMXPrecomp2, timeMMXFromRNS2, timeMMXToRNS2);
    cout << endl;
}

int main(int argc, char **argv)
{
    //srand((int)time(NULL));
    //srand48(time(NULL));
    static size_t iters = 10;
    static unsigned long b = 512;
    static unsigned long p = 20;
    static size_t m = 512;
    static bool longbench = false;
    static bool matmul = false;

    //static size_t k = 512 ;
    //static size_t n = 512 ;
    static Argument as[] = {
        {'b', "-b B", "Set the bitsize of the matrix entries.", TYPE_INT, &b},
        {'p', "-p B", "Set the bitsize of the RNS Prime.", TYPE_INT, &p},
        {'m', "-m M", "Set the dimension m of the matrix.", TYPE_INT, &m},
        {'i', "-i R", "Set minimal number of repetitions.", TYPE_INT, &iters},
        {'l', "-l Y", "Running a long benchmark .", TYPE_BOOL, &longbench},
        {'t', "-t Y", "Running a matmul size benchmark .", TYPE_BOOL, &matmul},
        END_OF_ARGUMENTS};
    FFLAS::parseArguments(argc, argv, as);

    size_t seed = time(NULL);
    typedef Givaro::Modular<Givaro::Integer> Field;

    cout << "### running RNS conversions benchmark ###" << endl;
    cout << "### timing in microseconds (from/to are given per element)" << endl;
    cout << "### |"
         << setw(SPC1) << "loop"
         << " |"
         << setw(SPC1) << "matrix dim"
         << " |"
         << setw(SPC1) << "bitsize"
         << " |"
         << setw(SPC1) << "mod bitsize"
         << " |"
         << setw(SPC3 + 2) << "FFLAS (precomp,from, to)"
         << " |"
         << setw(SPC3 + 2) << "FLINT (precomp,from, to)"
         << " |"
         << setw(SPC3 + 2) << "MMX   (precomp,from, to)"
         << " |"
         << setw(SPC3 + 2) << "MMX   (precomp,from, to)"
         << " |" << endl;

    if (!longbench)
        run_bench<Field>(m, p, b, iters, seed, matmul);
    else
    {
        size_t matdim = 128;
        std::vector<std::pair<size_t, size_t>> data = {
            {matdim, (1 << 7)},
            {matdim, (1 << 8)},
            {matdim, (1 << 9)},
            {matdim, (1 << 10)},
            {matdim, (1 << 11)},
            {matdim, (1 << 12)},
            {matdim, (1 << 13)},
            {matdim, (1 << 14)},
            {matdim, (1 << 15)},
            {matdim, (1 << 16)},
            {matdim, (1 << 17)},
            {matdim, (1 << 18)},
            { matdim, (1 << 18) },
        };

        for (size_t i = 0; i < data.size(); i++)
            run_bench<Field>(data[i].first, p, data[i].second, iters, seed++, matmul);
    }
    return 0;
}
