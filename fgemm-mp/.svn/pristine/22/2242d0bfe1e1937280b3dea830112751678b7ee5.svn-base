/* -*- mode: C++; tab-width: 8; indent-tabs-mode: t; c-basic-offset: 8 -*- */
/*
 * Copyright (C) 2013  Pascal Giorgi
 *
 * Written by Pascal Giorgi <pascal.giorgi@lirmm.fr>
 *
 * ========LICENCE========
 * This file is part of the library LinBox.
 *
 * LinBox is free software: you can redistribute it and/or modify
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
 */
#define SPC1 11
#define SPC2 22
#define PREC 6
#define BENCH_FLINT
#define BENCH_MMX
#define gettime usertime

#include <iostream>
#include <vector>
using namespace std;
#include <linbox/ring/modular.h>
#include <linbox/randiter/random-prime.h>
#include <linbox/randiter/random-fftprime.h>
//#include <linbox/field/unparametric.h>
#include <givaro/zring.h>
#include <linbox/matrix/matrix-domain.h>
#include <linbox/util/commentator.h>
#include <linbox/util/timer.h>
#include <linbox/matrix/polynomial-matrix.h>
#include <linbox/algorithms/polynomial-matrix/polynomial-matrix-domain.h>


#ifdef BENCH_FLINT
#define __GMP_BITS_PER_MP_LIMB 64
extern "C" {
#include "flint/longlong.h"
#include "flint/ulong_extras.h"
#include "flint/nmod_poly_mat.h"
#include "flint/flint.h"
}
#endif
#ifdef BENCH_MMX
#include <numerix/modular_int.hpp>
#include <algebramix/polynomial.hpp>
#include <algebramix/polynomial_modular_int.hpp>
#include <algebramix/polynomial_tft.hpp>
#include <algebramix/matrix.hpp>
#include <algebramix/matrix_modular_int.hpp>
#include <algebramix/matrix_tft.hpp>

#endif

using namespace LinBox;




template <typename Rand, typename Vect>
void randomVect (Rand& r, Vect& v) {
	size_t s = v.size();
	for (size_t i = 0; i < s; ++i)
		r.random(v[i]);
}

template <typename Rand, typename Mat>
void randomMat (Rand& r, Mat& m) {
	for (size_t i = 0; i < m.rowdim(); ++i)
		for (size_t j = 0; j < m.coldim(); ++j)
			r.random(m.refEntry(i,j));
}


template<typename Field, typename RandIter>
void profile_matpol_mulfft(const Field& fld,  RandIter& Gen, size_t n, size_t d, size_t iters) {
	typedef PolynomialMatrix<PMType::polfirst,PMStorage::plain,Field> MatrixP;

	double timeLinBox=0., timeFlint=0., timeMMX=0.;
	Givaro::Timer chrono;
	size_t loop=0;
	auto start = std::chrono::system_clock::now();
	auto end = std::chrono::system_clock::now();
	for (;(loop<iters ||std::chrono::duration<double>(end-start) < std::chrono::seconds(1)) ;loop++){

		MatrixP A(fld,n,n,d),B(fld,n,n,d),C(fld,n,n,2*d-1);
		// Generate random matrix of polynomial
		for (size_t i=0;i<n*n;i++){
			randomVect(Gen,A(i));
			randomVect(Gen,B(i));
		}	
		typedef PolynomialMatrixFFTMulDomain<Field>    FFT;
		FFT  PMFFT(fld);
		chrono.clear();
		chrono.start();
		PMFFT.mul(C,A,B);
		chrono.stop();
		timeLinBox+=chrono.gettime();

#ifdef BENCH_FLINT
		nmod_poly_mat_t AA,BB,CC;
		nmod_poly_mat_init(AA,n,n,(uint64_t)fld.cardinality());
		nmod_poly_mat_init(BB,n,n,(uint64_t)fld.cardinality());
		nmod_poly_mat_init(CC,n,n,(uint64_t)fld.cardinality());
		flint_rand_t state;
		flint_randinit(state);
		nmod_poly_mat_randtest(AA,state,d);
		nmod_poly_mat_randtest(BB,state,d);
		chrono.start();
		nmod_poly_mat_mul(CC,AA,BB);
		chrono.stop();
		timeFlint+=chrono.gettime();
#endif
	
#ifdef BENCH_MMX
		mmx::threads_set_number (1);
		srand(time(NULL));
		//typedef mmx::modulus<uint32_t,mmx::modulus_int_preinverse<29> > MOD;
		// typedef mmx::modulus<uint64_t > MOD;
		// typedef mmx::modular<MOD> COEFF;
		// MOD M((uint64_t)fld.cardinality());
		// COEFF::set_modulus(M);
		typedef mmx::modular<mmx::modulus<uint32_t,mmx::modulus_int_preinverse<29> >, mmx::modular_fixed<mmx::nat,469762049>> COEFF;
		//mmx::mmout <<"\n MMX mod :"<< COEFF::get_modulus()<<"\n";
		typedef mmx::polynomial_tft<mmx::polynomial_naive> PV;
		typedef mmx::matrix_tft<mmx::matrix_naive> MV;
		typedef mmx::polynomial<COEFF,PV> POLY;
		typedef mmx::matrix<POLY,MV> MATRIX;
		MATRIX AAA(1,n,n), BBB(1,n,n), CCC(1,n,n);
		for (mmx::nat i=0; i<n; i++)
			for (mmx::nat j=0; j<n; j++) {
				mmx::vector<COEFF> vA,vB;
				for(size_t h=0;h<d;h++){
					vA<<COEFF(rand()%(uint64_t)fld.cardinality());
					vB<<COEFF(rand()%(uint64_t)fld.cardinality());
				}
				AAA(i,j)=POLY(vA);
				BBB(i,j)=POLY(vB);
			}
		chrono.clear();
		chrono.start();
		CCC=AAA*BBB;
		//mmx::mmout<<CCC;
		chrono.stop();
		timeMMX+=chrono.gettime();
#endif
		end = std::chrono::system_clock::now();
	}
	
	cout<<setw(SPC1+5)
	    <<loop<<" |"
	    <<setw(SPC1)		
 	    <<fld.cardinality()<<" |"
	    <<setw(SPC1)
	    <<n<<" |"
	    <<setw(SPC1)
	    <<d-1<<" |"
	    <<setw(SPC1);
#define TIMING(x) cout<<setw(SPC2)<<std::fixed << std::setprecision(PREC)<<x<<" |";	
	TIMING(timeFlint/loop);
	TIMING(timeMMX/loop);
	TIMING(timeLinBox/loop);
	cout<<endl;
}




template<typename Field>
void runTest(const Field& F, size_t n, long b, long d, long seed, bool longtest, size_t iters){
	//typename Field::RandIter G(F,b,seed);
	typename Field::RandIter G(F,seed);
	if (!longtest)
		profile_matpol_mulfft(F,G,n,d,iters);
	else {
		size_t N[15]={16,16,16,16,16,16,16,16,             64,128,256,512,512,1024,2048};
		size_t D[15]={64,128,256,512,1024,2048,4096,8192,1024,512,256,512,128,  64,32 };
		for (size_t i=0;i<15;i++)
			profile_matpol_mulfft(F,G,N[i],D[i],iters);
	}
}

int main(int argc, char** argv){
	static size_t iters = 3 ;
	static size_t  n = 32; // matrix dimension
	static long    b = 20; // entries bitsize
	static uint64_t d = 32;  // matrix degree
	static long    seed = time(NULL);
	static bool    longbench=false;
	static bool    fftp=false;

	static Argument args[] = {
		{ 'n', "-n N", "Set dimension of test matrices to NxN." , TYPE_INT,  &n },
		{ 'd', "-d D", "Set degree of test matrices to D-1."    , TYPE_INT,  &d },
		{ 'b', "-b B", "Set bitsize of the matrix entries"      , TYPE_INT,  &b },		
		{ 'f', "-f Y", "Using a FFT Prime <2^29"                , TYPE_BOOL, &fftp },
		{ 's', "-s s", "Set the random seed to a specific value", TYPE_INT,  &seed},
		{ 'i', "-i R", "Set minimal number of repetitions.",                    TYPE_INT , &iters },
		{ 'l', "-l Y", "Running a long benchmark ."             , TYPE_BOOL, &longbench },
		END_OF_ARGUMENTS
	};
	parseArguments (argc, argv, args);


	cout<<"### running polynomial matrix multiplication benchmark ###"<<endl;
	cout<<"### |"
	    <<setw(SPC1)<<"loop"<<" |"
	    <<setw(b/3.32)<<"prime"<<" |"
	    <<setw(SPC1)<<"matrix dim"<<" |"
	    <<setw(SPC1)<<"degree"<<" |"
	    <<setw(SPC2)<<"time(FLINT)"<<" |"
	    <<setw(SPC2)<<"time(MMX)"<<" |"
	    <<setw(SPC2)<<"time(LinBox)"<<" |"<<endl;
	
	if (!fftp){
		RandomPrimeIter Rd(b,seed);
		integer p= Rd.random();			
		if (b>26) {			
			Givaro::Modular<integer> F(p);
			FFT_PROF_LEVEL=2;
			runTest (F,n,b,d,seed,longbench,iters);
		}
		else {
			Givaro::Modular<double> F(p);
			//Givaro::Modular<int32_t,int64_t> F((int32_t)p);
			runTest (F,n,b,d,seed,longbench,iters);		
		}
	}
	else {
		RandomFFTPrime Rd(1<<b,seed);
		integer p = Rd.randomPrime(integer(d).bitsize()+1);
		//Givaro::Modular<int32_t,int64_t> F((int32_t)p);
		Givaro::Modular<double> F((int32_t)p);
 		runTest (F,n,b,d,seed,longbench,iters);
	}

	return 0;
} 




