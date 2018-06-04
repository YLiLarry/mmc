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

#include <functional>
#include <iostream>
#include <vector>

#include <givaro/modular.h>
#include <givaro/givranditer.h>

using namespace std; 


#include "linbox/algorithms/polynomial-matrix/polynomial-fft-transform.h"
#include "linbox/randiter/random-fftprime.h"
#include "linbox/ring/modular.h"

using namespace LinBox;



template <typename Rand, typename Vect>
void randomVect (Rand& r, Vect& v) {
	size_t s = v.size();
	for (size_t i = 0; i < s; ++i)
		r.random(v[i]);
}


/**************************************
 ****** DFT PERFORMANCE FUNCTION ******
 **************************************/
template<typename Funct, typename FFT, typename Vect>
void DFT_performance(FFT& FFTDom, Funct f, size_t lpts, const Vect& x, string msg){
	Vect z(x);
	auto Functor = bind(f, &FFTDom, &z[0]);
	Timer chrono;
        double time;
        double Miops;
        size_t ct,minct=4;
	ct = 0;
	chrono.start();
	while (chrono.realElapsedTime() < 1){
		for (size_t k=0;k<minct;k++)
			Functor();
		ct+=minct;
	}
	time = chrono.userElapsedTime()/ct;
	Miops = 17 * (lpts<<(lpts-1)) /(1e6 * time); // 3/2 n log n
	msg+="  ";
	msg.resize(45,'.');
	cout << "Timings ... " << msg <<" : ";
	cout.precision(2);
	cout.width(10);
	cout<<scientific<<time << " s, ";
	cout.precision(2);
	cout.width(10);
	cout<<fixed<<Miops << " Miops\n";
}



template<typename Field>
void bench_DIF(const Field& fld, size_t kmax, long seed) { 
	typedef typename Field::Element Element;
	for (size_t lpts = 5; lpts < kmax ; lpts++){
		size_t pts = 1 << lpts;
		cout<<"********************************************************"<<endl;
		cout<<"*** Testing polynomials of size 2^" << lpts <<endl;
		cout<<"********************************************************"<<endl;
		vector<Element> x(pts);

		// Generate random inputs
		typename Field::RandIter Gen(fld,seed);
		randomVect(Gen,x);
		FFT_transform<Field> MulDom(fld,lpts);
		typedef FFT_transform<Field> FFT_t; 

		// check 2x2 		
		//DFT_performance(MulDom,&FFT_t::FFT_DIF_Harvey_mod2p_iterative2x2,lpts, x,     "DIF_Harvey_mod2p_iterative2x2");
		// check 3x3 		
		//DFT_performance(MulDom,&FFT_t::FFT_DIF_Harvey_mod2p_iterative3x3,lpts, x,     "DIF_Harvey_mod2p_iterative3x3");
		// check 4x1 SSE		
		//DFT_performance(MulDom,&FFT_t::FFT_DIF_Harvey_mod2p_iterative4x1_SSE,lpts, x, "DIF_Harvey_mod2p_iterative4x1_SSE");
		// check 4x2 SSE		
		//DFT_performance(MulDom,&FFT_t::FFT_DIF_Harvey_mod2p_iterative4x2_SSE,lpts, x, "DIF_Harvey_mod2p_iterative4x2_SSE");
#ifdef __AVX2__
		// check 8x1 AVX		
		//DFT_performance(MulDom,&FFT_t::FFT_DIF_Harvey_mod2p_iterative8x1_AVX,lpts, x, "DIF_Harvey_mod2p_iterative8x1_AVX");
#endif
		// check Harvey SSE		
		DFT_performance(MulDom,&FFT_t::template FFT_DIF<Element>,lpts, x, "DIF_Harvey_SSE");


		//cout<<"---------------------------------------------------------------"<<endl;
		// check 2x2 		
		//DFT_performance(MulDom,&FFT_t::FFT_DIT_Harvey_mod4p_iterative2x2,lpts, x,     "DIT_Harvey_mod4p_iterative2x2");
		// check 3x3 		
		//DFT_performance(MulDom,&FFT_t::FFT_DIT_Harvey_mod4p_iterative3x3,lpts, x,     "DIT_Harvey_mod4p_iterative3x3");
		// check 4x1 SSE		
		//DFT_performance(MulDom,&FFT_t::FFT_DIT_Harvey_mod4p_iterative4x1_SSE,lpts, x, "DIT_Harvey_mod4p_iterative4x1_SSE");
#ifdef __AVX2__
		// check 8x1 AVX		
		//DFT_performance(MulDom,&FFT_t::FFT_DIT_Harvey_mod4p_iterative8x1_AVX,lpts, x, "DIT_Harvey_mod4p_iterative8x1_AVX");
#endif
		// check Harvey SSE		
		DFT_performance(MulDom,&FFT_t::template FFT_DIT<Element>,lpts, x, "DIT_Harvey_SSE");


		cout<<endl;
	}
}


int main(int argc, char** argv){
	if (argc < 2 || argc >3){
		cerr<<"usage : prime_bitsize , (seed)"<<endl;
		exit(0);
	}
	uint64_t bits =atoi(argv[1]);
	long seed=((argc>2)?atoi(argv[2]):time(NULL));	
	RandomFFTPrime Rd(1<<bits,seed);
	uint32_t p = (uint32_t)Rd.randomPrime(bits-3);
	size_t k = 0;
	size_t _p=p,_mask=2;
	while(!((_p-1)&(_mask-1))) {_mask<<=1;k++;}
	cout<<"p       = "<<p<<endl;
	cout<<"val2(p) = "<<k<<endl;
	cout<<endl;
	
	//Givaro::Modular<uint32_t,uint64_t> F(p);
	Givaro::Modular<double> F(p);
	bench_DIF(F,k,seed);


	return 0;
}
 
 
