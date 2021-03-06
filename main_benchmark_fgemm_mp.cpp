/* -*- mode: C++; tab-width: 4; indent-tabs-mode: t; c-basic-offset: 4 -*- */
// vim:sts=4:sw=4:ts=4:noet:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s

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

#if not defined(MG_DEFAULT)
#define MG_DEFAULT MG_ACTIVE
#endif
#if not defined(STD_RECINT_SIZE)
#define STD_RECINT_SIZE 8
#endif

#include "fflas-ffpack/fflas-ffpack-config.h"
#include <iostream>
#include <typeinfo>
#include <vector>
#include <string>
using namespace std;

#include "fflas-ffpack/fflas/fflas.h"
#include "fflas-ffpack/utils/fflas_io.h"
#include "fflas-ffpack/utils/timer.h"
#include "fflas-ffpack/utils/args-parser.h"
#include <gmp++/gmp++.h>
#include "givaro/givcaster.h"
#include "fflas-ffpack/paladin/parallel.h"
#ifdef BENCH_RECINT
#include "recint/recint.h"
#endif
#include "two_phase_marge_least.h"
#include "two_phase_marge_most.h"
#include "two_phase_parge_block.h"
#include "two_phase_parge_shift.h"
#include <ostream>

#include "containers.h"
#include "sim_rns.h"
#include <algorithm>

#ifdef BENCH_FLINT
#define __GMP_BITS_PER_MP_LIMB 64
extern "C"
{
#include "flint/longlong.h"
#include "flint/long_extras.h"
#include "flint/fmpz_mat.h"
#include "flint/fmpz.h"
#include "flint/flint.h"
}
#endif


#define BENCH_FFLAS_PPACK 0
#define BENCH_TWO_PHASE_PARGE_SHIFT 0
#define BENCH_TWO_PHASE_PARGE_BLOCK 0
#define BENCH_TWO_PHASE_MARGE_LEAST 1
#define BENCH_TWO_PHASE_MARGE_MOST 1

static size_t iters = 1;
static Givaro::Integer q = -1;
static unsigned long b = 17;
static size_t d = 64;
static size_t e = 13;
static int nbw = -1;
static size_t seed = time(NULL);
static Argument as[] = {
    {'q', "-q Q", "Set the field characteristic (-1 for random).", TYPE_INTEGER, &q},
    {'b', "-b B", "Set the bitsize of the random characteristic.", TYPE_INT, &b},
    {'e', "-e E", "Set the bitsize of the level 1 moduli.", TYPE_INT, &e},
    {'d', "-d D", "Set the dimension of matrices.", TYPE_INT, &d},
    // {'m', "-m M", "Set the dimension m of the matrix.", TYPE_INT, &m},
    // {'k', "-k K", "Set the dimension k of the matrix.", TYPE_INT, &k},
    // {'n', "-n N", "Set the dimension n of the matrix.", TYPE_INT, &n},
    {'w', "-w N", "Set the number of winograd levels (-1 for random).", TYPE_INT, &nbw},
    {'i', "-i R", "Set number of repetitions.", TYPE_INT, &iters},
    {'s', "-s S", "Sets seed.", TYPE_INT, &seed},
    END_OF_ARGUMENTS};

template <typename Ints>
int tmain()
{

    static size_t m = d;
    static size_t k = d;
    static size_t n = d;

    srand((int)seed);
    srand48(seed);
    Givaro::Integer::seeding(seed);

    typedef Givaro::ZRing<Ints> Field;
    typedef Givaro::Modular<Ints> ModularField;
    Givaro::Integer p = 1;
    p <<= (1 << b);
    FFLAS::Timer chrono, TimFreivalds;
    double time = 0., timev = 0.;
#ifdef BENCH_FLINT
    double timeFlint = 0.;
#endif
#if BENCH_TWO_PHASE_MARGE_LEAST || BENCH_TWO_PHASE_MARGE_MOST || BENCH_TWO_PHASE_PARGE_BLOCK || BENCH_TWO_PHASE_PARGE_SHIFT
    double time_two_phase_marge_least = 0.;
    double time_two_phase_marge_most = 0.;
    double time_two_phase_parge_block = 0.;
    double time_two_phase_parge_shift = 0.;
#endif
    double time_fflas_ppack = 0.;
    for (size_t loop = 0; loop < iters; loop++)
    {
        // Givaro::Integer::random_exact_2exp(p, 1 << b);

        // Givaro::IntPrimeDom IPD;
        // IPD.prevprimein(p);
        // Ints ip;
        // Givaro::Caster<Ints, Givaro::Integer>(ip, p);
        // Givaro::Caster<Givaro::Integer, Ints>(p, ip); // to check consistency

        Field F;
        ModularField MF(p);
        size_t lda, ldb, ldc;
        lda = k;
        ldb = n;
        ldc = n;

        typename ModularField::RandIter Rand(MF, seed);
        typename Field::Element_ptr A, B, C;
        A = FFLAS::fflas_new(F, m, lda);
        B = FFLAS::fflas_new(F, k, ldb);
        C = FFLAS::fflas_new(F, m, ldc);

        // 		for (size_t i=0;i<m;++i)
        // 			for (size_t j=0;j<k;++j)
        // 				Rand.random(A[i*lda+j]);
        // 		for (size_t i=0;i<k;++i)
        // 			for (size_t j=0;j<n;++j)
        // 				Rand.random(B[i*ldb+j]);
        // 		for (size_t i=0;i<m;++i)
        // 			for (size_t j=0;j<n;++j)
        // 				Rand.random(C[i*ldc+j]);

        PAR_BLOCK { FFLAS::pfrand(F, Rand, m, k, A, m / size_t(MAX_THREADS)); }
        PAR_BLOCK { FFLAS::pfrand(F, Rand, k, n, B, k / MAX_THREADS); }
        PAR_BLOCK { FFLAS::pfzero(F, m, n, C, m / MAX_THREADS); }

        Ints alpha, beta;
        alpha = F.one;
        beta = F.zero;

        vector<Givaro::Integer> A_(A, A + m * k);
        vector<Givaro::Integer> B_(B, B + k * n);

        vector<Givaro::Integer> M_;
        M_.insert(M_.end(), A_.begin(), A_.end());
        M_.insert(M_.end(), B_.begin(), B_.end());
        
        vector<size_t> M_s;
        M_s.push_back(m);
        M_s.push_back(k);
        M_s.push_back(n);

        Givaro::Integer input_bitsize = max(max_element(A_.begin(), A_.end()), max_element(B_.begin(), B_.end()))->bitsize();
        cerr << "benchmarking with random inputs in Givaro::Modular<Ints>(integers of 2^" << input_bitsize.bitsize() - 1<< " bits)" << endl;

        // cerr << A_ << endl;
        // cerr << B_ << endl;

#ifdef BENCH_FLINT
        // FLINT MUL //
        fmpz_t modp, tmp;
        fmpz_init(modp);
        fmpz_init(tmp);
        fmpz_set_mpz(modp, *(reinterpret_cast<const mpz_t *>(&p)));
        fmpz_mat_t AA, BB, CC, DD;
        fmpz_mat_init(AA, m, k);
        fmpz_mat_init(BB, k, n);
        fmpz_mat_init(CC, m, n);
        fmpz_mat_init(DD, m, n);
        fmpz_t aalpha, bbeta;
        fmpz_set_mpz(aalpha, *(reinterpret_cast<const mpz_t *>(&alpha)));
        fmpz_set_mpz(bbeta, *(reinterpret_cast<const mpz_t *>(&beta)));

        for (size_t i = 0; i < m; ++i)
            for (size_t j = 0; j < k; ++j)
                fmpz_set_mpz(fmpz_mat_entry(AA, i, j), *(reinterpret_cast<const mpz_t *>(A + i * lda + j)));
        for (size_t i = 0; i < k; ++i)
            for (size_t j = 0; j < n; ++j)
                fmpz_set_mpz(fmpz_mat_entry(BB, i, j), *(reinterpret_cast<const mpz_t *>(B + i * ldb + j)));
        for (size_t i = 0; i < m; ++i)
            for (size_t j = 0; j < n; ++j)
                fmpz_set_mpz(fmpz_mat_entry(CC, i, j), *(reinterpret_cast<const mpz_t *>(C + i * ldc + j)));
        chrono.clear();
        chrono.start();
        // DD= A.B
        fmpz_mat_mul(DD, AA, BB);
        // CC = beta.C
        fmpz_mat_scalar_mul_fmpz(CC, CC, bbeta);
        // CC = CC + DD.alpha
        fmpz_mat_scalar_addmul_fmpz(CC, DD, aalpha);
        // CC = CC mod p
        for (size_t i = 0; i < m; ++i)
            for (size_t j = 0; j < n; ++j)
                fmpz_mod(fmpz_mat_entry(CC, i, j), fmpz_mat_entry(CC, i, j), modp);

        chrono.stop();
        timeFlint += chrono.usertime();
        fmpz_mat_clear(AA);
        fmpz_mat_clear(BB);
#endif

#if BENCH_FFLAS_PPACK
        cerr << "===========================================" << endl;
        cerr << "========== Benchmark FFLAS-PPACK ==========" << endl;
        cerr << "===========================================" << endl;
        chrono.clear();
        chrono.start();
        auto C__ = SIM_RNS::fflas_mult_integer(A_, B_, m, k, n);
        chrono.stop();
        time_fflas_ppack += chrono.usertime();
#endif

#if BENCH_TWO_PHASE_PARGE_SHIFT
        {
            cerr << "===========================================" << endl;
            cerr << "======= Benchmark TwoPhasePargeShift =======" << endl;
            cerr << "===========================================" << endl;
            chrono.clear();
            chrono.start();
            TwoPhasePargeShift algo(2 * input_bitsize, input_bitsize / 2, 1);
            auto a = algo.matrix_reduce(A_, m, k);
            auto b = algo.matrix_reduce(B_, k, n);
            auto c = algo.phase2_mult(a, b);
            auto C_ = algo.matrix_recover(c);
            chrono.stop();
            time_two_phase_parge_shift += chrono.usertime();
            // this line asserts our result is the same as FFLAS::fgemm
            #if BENCH_FFLAS_PPACK
            if (! equals(C_, C__)) {
                cerr << "ERROR! RESULT IS INCORRECT" << endl
                     << " - expected: " << C__ << endl
                     << " - got: " << C_ << endl;
            }
            #endif
        }
#endif

#if BENCH_TWO_PHASE_PARGE_BLOCK
        {
            cerr << "===========================================" << endl;
            cerr << "====== Benchmark TwoPhasePargeBlock =======" << endl;
            cerr << "===========================================" << endl;
            chrono.clear();
            chrono.start();
            TwoPhasePargeBlock algo(2 * input_bitsize, 1 << e, 4);
            auto a = algo.matrix_reduce(A_, m, k);
            auto b = algo.matrix_reduce(B_, k, n);
            auto c = algo.phase2_mult(a, b);
            auto C_ = algo.matrix_recover(c);
            chrono.stop();
            time_two_phase_parge_block += chrono.usertime();
            // this line asserts our result is the same as FFLAS::fgemm
            #if BENCH_FFLAS_PPACK
            if (! equals(C_, C__)) {
                cerr << "ERROR! RESULT IS INCORRECT" << endl
                     << " - expected: " << C__ << endl
                     << " - got: " << C_ << endl;
            };
            #endif
        }
#endif

#if BENCH_TWO_PHASE_MARGE_LEAST
        {
            cerr << "===========================================" << endl;
            cerr << "======= Benchmark TwoPhaseMargeLeast ======" << endl;
            cerr << "===========================================" << endl;
            chrono.clear();
            chrono.start();
            TwoPhaseMargeLeast algo(2 * input_bitsize, 1 << e);
            auto matrices = algo.matrix_reduce(M_, M_s);
            auto a = matrices[0];
            auto b = matrices[1];
            auto c = algo.phase2_mult(a, b);
            auto C_ = algo.matrix_recover(c);
            chrono.stop();
            time_two_phase_marge_least += chrono.usertime();
            // this line asserts our result is the same as FFLAS::fgemm
            #if BENCH_FFLAS_PPACK
            if (! equals(C_, C__)) {
                cerr << "ERROR! RESULT IS INCORRECT" << endl
                     << " - expected: " << C__ << endl
                     << " - got: " << C_ << endl;
            }
            #endif
        }
#endif

#if BENCH_TWO_PHASE_MARGE_MOST
        {
            cerr << "===========================================" << endl;
            cerr << "======= Benchmark TwoPhaseMargeMost =======" << endl;
            cerr << "===========================================" << endl;
            chrono.clear();
            chrono.start();
            TwoPhaseMargeMost algo(input_bitsize << 1, 1 << e);
            auto a = algo.matrix_reduce(A_, m, k);
            auto b = algo.matrix_reduce(B_, k, n);
            auto c = algo.phase2_mult(a, b);
            auto C_ = algo.matrix_recover(c);
            chrono.stop();
            time_two_phase_marge_most += chrono.usertime();
            // this line asserts our result is the same as FFLAS::fgemm
            #if BENCH_FFLAS_PPACK
            if (! equals(C_, C__)) {
                cerr << "ERROR! RESULT IS INCORRECT" << endl
                     << " - expected: " << C__ << endl
                     << " - got: " << C_ << endl;
            }
            #endif
        }
#endif


/*
        //END FLINT CODE //
        using FFLAS::CuttingStrategy::Recursive;
        using FFLAS::StrategyParameter::TwoDAdaptive;
        // RNS MUL_LA
        Givaro::ZRing<Givaro::Integer> intF;
        chrono.clear();
        chrono.start();
        // 		PAR_BLOCK{
        //             FFLAS::fgemm(F,FFLAS::FflasNoTrans,FFLAS::FflasNoTrans,m,n,k,alpha,A,lda,B,ldb,beta,C,ldc, SPLITTER(NUM_THREADS,Recursive,TwoDAdaptive) );
        // 		}
        {
            FFLAS::fgemm(F, FFLAS::FflasNoTrans, FFLAS::FflasNoTrans, m, n, k, alpha, A, lda, B, ldb, beta, C, ldc, FFLAS::ParSeqHelper::Sequential());
        }

        chrono.stop();
        time += chrono.usertime();

        // FFLAS::WriteMatrix(cerr, F, m, k, C, k);
        // cerr << C_ << endl;

        TimFreivalds.start();
        bool pass = FFLAS::freivalds(F, FFLAS::FflasNoTrans, FFLAS::FflasNoTrans, m, n, k, alpha, A, k, B, n, C, n);
        TimFreivalds.stop();
        timev += TimFreivalds.usertime();
        if (!pass)
        {
            std::cout << "FAILED" << std::endl;
            std::cout << "p:=" << p << ';' << std::endl;
            FFLAS::WriteMatrix(std::cout << "A:=", F, m, k, A, lda) << ';' << std::endl;
            FFLAS::WriteMatrix(std::cout << "B:=", F, k, n, B, ldb) << ';' << std::endl;
            FFLAS::WriteMatrix(std::cout << "C:=", F, m, n, C, ldc) << ';' << std::endl;
        }

        FFLAS::fflas_delete(A);
        FFLAS::fflas_delete(B);
        FFLAS::fflas_delete(C);
*/
    }
    double Gflops = (2. * double(m) / 1000. * double(n) / 1000. * double(k) / 1000.0) / time * double(iters);
    // 	Gflops*=p.bitsize()/16.;
    cout << "Time fgemm: " << (time / double(iters))
         << " Gfops: " << Gflops
         << " (total:" << time << ") "
         << typeid(Ints).name()
         << " | perword: " << (Gflops * double(p.bitsize())) / 64.;

    std::cout << " | "
              << "benchmarked with random inputs in Givaro::Modular<Ints>(integer of " << p.bitsize() << " bits) | ";
    FFLAS::writeCommandString(std::cout, as);
    std::cout << "  | Freivalds: " << timev / double(iters) << std::endl;

#ifdef BENCH_FLINT
    cout << "Time FLINT: " << timeFlint << endl;
#endif

#if BENCH_TWO_PHASE_MARGE_LEAST || BENCH_TWO_PHASE_MARGE_MOST || BENCH_TWO_PHASE_PARGE_BLOCK
    cout << "Time TwoPhaseMargeLeast: " << time_two_phase_marge_least << endl;
    cout << "Time TwoPhaseMargeMost: " << time_two_phase_marge_most << endl;
    cout << "Time TwoPhasePargeBlock: " << time_two_phase_parge_block << endl;
    cout << "Time TwoPhasePargeShift: " << time_two_phase_parge_shift << endl;
#endif
    cout << "Time FFLAS-PPACK: " << time_fflas_ppack << endl;

    return 0;
}

int main(int argc, char **argv)
{
    FFLAS::parseArguments(argc, argv, as);

    int r1 = tmain<Givaro::Integer>();

#ifdef BENCH_RECINT
    r1 += tmain<RecInt::rint<STD_RECINT_SIZE>>();
#endif
    return r1;
}
