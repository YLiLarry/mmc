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

#include "fflas-ffpack/fflas-ffpack-config.h"
#include <iostream>
#include <vector>
#include <string>
#include <chrono>
using namespace std;

#include "fflas-ffpack/utils/timer.h"
#include "fflas-ffpack/fflas/fflas.h"
#include "fflas-ffpack/utils/args-parser.h"
#include <gmp++/gmp++.h>givaro/zring.h"

#define BENCH_KRONECKER
#define BENCH_FLINT
#define BENCH_MMX

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

#ifdef BENCH_MMX
#include "algebramix/matrix.hpp"
#include "algebramix/matrix_integer.hpp"
#endif

#ifdef BENCH_KRONECKER
#include "kroneckerFFT.h"
#endif
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
void run_bench(size_t m, size_t n, size_t k, size_t b, size_t iters, size_t seed)
{
	Givaro::Integer p;
	FFLAS::Timer chrono;
	double time = 0.;
#ifdef BENCH_FLINT
	double timeFlint1 = 0.;
	double timeFlint2 = 0.;
#endif
#ifdef BENCH_MMX
	double timeMMX1 = 0;
	double timeMMX2 = 0;
#endif
#ifdef BENCH_KRONECKER
	double timeKronecker = 0.;
#endif

	size_t loop = 0;
	auto start = std::chrono::system_clock::now();
	auto end = std::chrono::system_clock::now();
	for (; (loop < iters || std::chrono::duration<double>(end - start) < std::chrono::seconds(1)); loop++)
	{
		//cout<<loop<<" "<<end-start<<endl;
		Givaro::Integer::random_exact_2exp(p, b);
		//nextprime(p,p);
		//cout<<p<<endl;

		Field F(p);
		size_t lda, ldb, ldc;
		lda = k;
		ldb = n;
		ldc = n;

		typename Field::RandIter Rand(F, seed);
		typename Field::Element_ptr A, B, C, C2;
		A = FFLAS::fflas_new(F, m, lda);
		B = FFLAS::fflas_new(F, k, ldb);
		C = FFLAS::fflas_new(F, m, ldc);
		C2 = FFLAS::fflas_new(F, m, n);

		for (size_t i = 0; i < m; ++i)
			for (size_t j = 0; j < k; ++j)
				Rand.random(A[i * lda + j]);
		for (size_t i = 0; i < k; ++i)
			for (size_t j = 0; j < n; ++j)
				Rand.random(B[i * ldb + j]);
		for (size_t i = 0; i < m; ++i)
			for (size_t j = 0; j < n; ++j)
				Rand.random(C[i * ldc + j]);

		Givaro::Integer alpha, beta;
		alpha = 1;
		beta = 0;

#ifdef BENCH_FLINT
		{
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
			fmpz_mat_mul_classical_inline(DD, AA, BB);
			chrono.stop();
			timeFlint1 += chrono.usertime();
			chrono.clear();
			chrono.start();
			fmpz_mat_clear(AA);
			fmpz_mat_clear(BB);
		}
		{
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
			fmpz_mat_mul_multi_mod(DD, AA, BB);
			chrono.stop();
			timeFlint2 += chrono.usertime();
			chrono.clear();
			chrono.start();
			fmpz_mat_clear(AA);
			fmpz_mat_clear(BB);
		}

#endif
		//END FLINT CODE //

#ifdef BENCH_MMX
		{
			using namespace mmx;
			typedef mmx::integer mT;
			typedef matrix_integer mM;
			typedef implementation<matrix_multiply, matrix_naive> Mat;
			//  MMX MUL //

			threads_set_number(1);
			matrix<mT, mM> AAA(1, m, k), BBB(1, k, n);
			for (size_t i = 0; i < m; ++i)
				for (size_t j = 0; j < k; ++j)
					mpz_set(*AAA(i, j), *(reinterpret_cast<const mpz_t *>((A + i * lda + j))));
			for (size_t i = 0; i < k; ++i)
				for (size_t j = 0; j < n; ++j)
					mpz_set(*BBB(i, j), *(reinterpret_cast<const mpz_t *>((B + i * lda + j))));
			typedef implementation<matrix_multiply, mM> MUL;
			nat l = aligned_size<mT, mM>(m * n);
			mT *d = mmx_formatted_new<mT>(l, format<mT>());

			xnat sz1 = matrix_crt_multiply_helper<mT>::size(tab(AAA), Mat::index(1, 0, m, k), Mat::index(0, 1, m, k), m, k);
			xnat sz2 = matrix_crt_multiply_helper<mT>::size(tab(BBB), Mat::index(1, 0, k, n), Mat::index(0, 1, k, n), k, n);
			xnat sz = matrix_crt_multiply_helper<mT>::size(tab(AAA), Mat::index(1, 0, m, k), Mat::index(0, 1, m, k),
														   tab(BBB), Mat::index(1, 0, k, n), Mat::index(0, 1, k, n),
														   m, k, n);

			//cout<<"MMX...";
			chrono.clear();
			chrono.start();
			(void)MUL::mul_fft_3_int(d, tab(AAA), tab(BBB), m, k, n, sz, sz1, sz2);
			chrono.stop();
			//cout<<"DONE"<<endl;

			timeMMX1 += chrono.usertime();
			end = std::chrono::system_clock::now();
			mmx_delete<mT>(d, l);
		}

		{
			using namespace mmx;
			typedef mmx::integer mT;
			typedef matrix_integer_crt mM;
			//  MMX MUL //
			threads_set_number(1);
			matrix<mT, mM> AAA(1, m, k), BBB(1, k, n);
			for (size_t i = 0; i < m; ++i)
				for (size_t j = 0; j < k; ++j)
					mpz_set(*AAA(i, j), *(reinterpret_cast<const mpz_t *>((A + i * lda + j))));
			for (size_t i = 0; i < k; ++i)
				for (size_t j = 0; j < n; ++j)
					mpz_set(*BBB(i, j), *(reinterpret_cast<const mpz_t *>((B + i * lda + j))));
			typedef implementation<matrix_multiply, mM> MUL;
			nat l = aligned_size<mT, mM>(m * n);
			mT *d = mmx_formatted_new<mT>(l, format<mT>());
			//cout<<"MMX...";
			chrono.clear();
			chrono.start();
			(void)MUL::mul(d, tab(AAA), tab(BBB), m, k, n);
			chrono.stop();
			//cout<<"DONE"<<endl;

			timeMMX2 += chrono.usertime();
			end = std::chrono::system_clock::now();
			mmx_delete<mT>(d, l);
		}

#endif //END MMX CODE //

		Givaro::ZRing<Givaro::Integer> Z;

#ifndef DISABLE_FGEMM
		// RNS MUL_LA
		chrono.clear();
		chrono.start();
		FFLAS::fgemm(Z, FFLAS::FflasNoTrans, FFLAS::FflasNoTrans, m, n, k, alpha, A, lda, B, ldb, beta, C, ldc);
		chrono.stop();
		time += chrono.usertime();
#endif
#ifdef BENCH_KRONECKER
		chrono.clear();
		chrono.start();
		intmatmul(m, n, k, C2, A, B);
		chrono.stop();
		timeKronecker += chrono.usertime();

		// if (FFLAS::fequal(Z,m,n,C2,n,C,ldc))
		// 	cout<<"Le Gros Necker est correct"<<endl;
		// else{
		// 	cout<<"Le Gros Necker est FAUX !!!"<<endl;
		//cout<<A[0]<<"x"<<B[0]<<endl;
		//cout<<C2[0]<<endl;
		//cout<<"diff taille: "<<C[0].size()<<"!="<<C2[0].size()<<endl;
		//cout<<"diff valeur: "<<C[0]-C2[0]<<endl;
		//}
#endif

		FFLAS::fflas_delete(A);
		FFLAS::fflas_delete(B);
		FFLAS::fflas_delete(C);
		FFLAS::fflas_delete(C2);
	}
#define SPC1 11
#define SPC2 22
#define PREC 6

	cout << setw(SPC1 + 5)
		 << loop << " |"
		 << setw(SPC1)
		 << m << " |"
		 << setw(SPC1)
		 << b << " |"
		 << setw(SPC1);

#define TIMING(x) cout << setw(SPC2) << std::fixed << std::setprecision(PREC) << x << " |";

	TIMING(time / loop);
	TIMING(timeFlint1 / loop);
	TIMING(timeFlint2 / loop);
	TIMING(timeMMX1 / loop);
	TIMING(timeMMX2 / loop);
	TIMING(timeKronecker / loop);
	cout << endl;
}

int main(int argc, char **argv)
{
	//srand((int)time(NULL));
	//srand48(time(NULL));
	static size_t iters = 3;
	static Givaro::Integer q = -1;
	static unsigned long b = 512;
	static size_t m = 512;
	static bool longbench = false;
	//static size_t k = 512 ;
	//static size_t n = 512 ;
	static Argument as[] = {
		{'q', "-q Q", "Set the field characteristic (-1 for random).", TYPE_INTEGER, &q},
		{'b', "-b B", "Set the bitsize of the random characteristic.", TYPE_INT, &b},
		{'m', "-m M", "Set the dimension m of the matrix.", TYPE_INT, &m},
		//{ 'k', "-k K", "Set the dimension k of the matrix.",                    TYPE_INT , &k },
		//{ 'n', "-n N", "Set the dimension n of the matrix.",                    TYPE_INT , &n },
		{'i', "-i R", "Set minimal number of repetitions.", TYPE_INT, &iters},
		{'l', "-l Y", "Running a long benchmark .", TYPE_BOOL, &longbench},
		END_OF_ARGUMENTS};
	FFLAS::parseArguments(argc, argv, as);

	size_t k, n;
	k = n = m;
	size_t seed = time(NULL);
	typedef Givaro::Modular<Givaro::Integer> Field;

	cout << "### running fgemm multiprecision benchmark ###" << endl;
	cout << "### |"
		 << setw(SPC1) << "loop"
		 << " |"
		 << setw(SPC1) << "matrix dim"
		 << " |"
		 << setw(SPC1) << "bitsize"
		 << " |"
		 << setw(SPC2) << "time(FFLAS)"
		 << " |"
		 << setw(SPC2) << "time(FLINT classic)"
		 << " |"
		 << setw(SPC2) << "time(FLINT crt)"
		 << " |"
		 << setw(SPC2) << "time(MMX kron)"
		 << " |"
		 << setw(SPC2) << "time(MMX crt)"
		 << " |"
		 << setw(SPC2) << "time(LinBox)"
		 << " |" << endl;
	if (!longbench)
		run_bench<Field>(m, n, k, b, iters, seed);
	else
	{

		std::vector<std::pair<size_t, size_t>> data = {
			{16, (1 << 10)},
			{16, (1 << 11)},
			{16, (1 << 12)},
			{16, (1 << 13)},
			{16, (1 << 14)},
			{16, (1 << 15)},
			{16, (1 << 16)},
			{16, (1 << 17)},
			{16, (1 << 18)},

			{64, (1 << 10)},
			{64, (1 << 11)},
			{64, (1 << 12)},
			{16, (1 << 13)},
			{64, (1 << 14)},
			{64, (1 << 15)},
			{64, (1 << 16)},
			{64, (1 << 17)},
			{64, (1 << 18)},

			{256, (1 << 10)},
			{256, (1 << 11)},
			{256, (1 << 12)},
			{256, (1 << 13)},
			{256, (1 << 14)},
			{256, (1 << 15)},
			{256, (1 << 16)},
			{256, (1 << 17)},

			{512, (1 << 10)},
			{512, (1 << 11)},
			{512, (1 << 12)},
			{16, (1 << 13)},
			{512, (1 << 14)},
			{512, (1 << 15)},
			{512, (1 << 16)},

			{1024, (1 << 10)},
			{1024, (1 << 11)},
			{1024, (1 << 12)},
			{1024, (1 << 13)},
			{1024, (1 << 14)},
			{1024, (1 << 15)},
		};

		for (size_t i = 0; i < data.size(); i++)
			run_bench<Field>(data[i].first, data[i].first, data[i].first, data[i].second, iters, seed++);
			run_bench<Field>(data[i].first,data[i].first,data[i].first,data[i].second,iters,seed++);
	}
	return 0;
}
