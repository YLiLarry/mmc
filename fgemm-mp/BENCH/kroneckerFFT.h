/* -*- mode: C++; tab-width: 8; indent-tabs-mode: t; c-basic-offset: 8 -*- */
/* 
 * Copyright (C) 2013  Romain Lebreton  <lebreton@lirmm.fr>
 *               2015  Pascal Giorgi <giorgi@lirmm.fr> 
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
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.	 See the GNU
 * Lesser General Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public
 * License along with this library; if not, write to the Free Software
 * Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301  USA
 * ========LICENCE========
 */

#ifndef KRONECKERFFT_H
#define KRONECKERFFT_H
#include <cstddef>
#include "linbox/matrix/polynomial-matrix.h"
#include "linbox/matrix/dense-matrix.h"
#include "linbox/randiter/random-fftprime.h"
#include "gmp++/gmp++.h"
#include "linbox/util/timer.h"
#include "givaro/modular.h"
#include <vector>
#include <algorithm>
#include "linbox/algorithms/polynomial-matrix/matpoly-mult-fft.h"
#include "fflas-ffpack/utils/Matio.h"

typedef LinBox::Integer Integer;

typedef Givaro::Modular<double> ModField;
typedef LinBox::PolynomialMatrix<LinBox::PMType::polfirst, LinBox::PMStorage::plain, ModField> MatrixP;

// input:
//    - lb :  2^lb is a bound on the Integer matrix to be computed
//    - n  : dimension of the matrix product
// output:
//    - d : the size of the kronecker transform in base 2^16
//    - primes: a list of fourier prime that allow FFT for dxd polynomial product
void kronecker_init(uint64_t &d, vector<Integer> &primes, uint64_t n, uint64_t lb)
{
	uint64_t ln = Integer((uint64_t)n).bitsize();
	uint64_t lfflas = (53 - ln) / 2; // TODO Check !
	// We will take 2^(lfflas-1) <= p_i <x 2^lfflas

	uint64_t lp = 16;
	d = ceil((double)lb / lp);

	Integer bound = Integer(n) * Integer(d) * Integer(uint64_t(1) << 32);
	uint64_t prime_max = std::sqrt((1UL << 53) / n) + 1;
	if (prime_max > (1 << 26))
		prime_max = (1 << 26) - 1;
	if (bound.bitsize() < 64)
	{
#ifdef VERBOSE
		cout << "bound bitsize: " << bound.bitsize() << endl;
		cout << "fflas pmax: " << prime_max << endl;
		cout << "fflas pbit: " << Integer(prime_max).bitsize() << endl;
#endif
		if (bound.bitsize() >= 2 * lfflas)
		{
			lfflas = 21;
#ifdef VERBOSE
			cout << "fflas bit: " << lfflas << endl;
#endif
		}
	}

	// Find FFT primes pi that allow 2^lpts transform and s.t. pi<2^lfflas and prod(pi)<bound
	uint64_t lpts = Integer(2 * (d - 1)).bitsize();

	LinBox::RandomFFTPrime RdFFT(prime_max);
	//LinBox::RandomFFTPrime RdFFT (lfflas);

	if (!RdFFT.generatePrimes(lpts, bound, primes))
	{
		std::cout << "COULD NOT FIND ENOUGH FFT PRIME for IntegerMul kronecker" << std::endl;
		exit(1);
	}
}

void to_kronecker_16(double *out, const Integer &in, uint64_t d)
{
	const mpz_t *m0 = reinterpret_cast<const mpz_t *>(&in);
	const uint16_t *m0_ptr = reinterpret_cast<const uint16_t *>(m0[0]->_mp_d);
	uint64_t tmp = (in.size()) * sizeof(mp_limb_t) / 2;
	uint64_t maxs = std::min(d, tmp); // to ensure 32 bits portability
	uint64_t l = 0;
	for (; l < maxs; l++)
		out[l] = m0_ptr[l];
	for (; l < d; l++)
		out[l] = 0.;
}

void from_kronecker_16(Integer &out, const Integer *in, uint64_t d)
{
	out = 0;
	for (uint64_t i = 0; i < d; i++)
		out += in[i] << ((uint64_t)i << 2);
}

void from_kronecker_16(uint64_t N, Integer *out, const uint64_t *in, uint64_t d, uint64_t stride, uint64_t rns_bitsize)
{
	uint64_t k = d;
	uint64_t cst = 0; //(rns_bitsize>>4) -1;//3;
	uint64_t k4 = ((k + cst) >> 2) + (((k + cst) % 4 == 0) ? 0 : 1);
	//cout<<"k4: "<<k4<<endl;
	//uint64_t outbit = (d-1)*16+64;
	//cout<<"outbit:"<<outbit<<endl;
	//k4= outbit/64  + (outbit%64==0?0:1);

	std::vector<uint16_t> A0(k4 << 2, 0), A1(k4 << 2, 0), A2(k4 << 2, 0), A3(k4 << 2, 0);
	Integer a0, a1, a2, a3, res;
	mpz_t *m0, *m1, *m2, *m3;
	//cout<<"k4: "<<k4<<endl;
	//cout<<"mpz alloc: "<<(int) (k4*8/sizeof(mp_limb_t))<<endl;
	m0 = reinterpret_cast<mpz_t *>(&a0);
	m1 = reinterpret_cast<mpz_t *>(&a1);
	m2 = reinterpret_cast<mpz_t *>(&a2);
	m3 = reinterpret_cast<mpz_t *>(&a3);
	mp_limb_t *m0_d, *m1_d, *m2_d, *m3_d;
	m0_d = m0[0]->_mp_d;
	m1_d = m1[0]->_mp_d;
	m2_d = m2[0]->_mp_d;
	m3_d = m3[0]->_mp_d;
	m0[0]->_mp_alloc = m1[0]->_mp_alloc = m2[0]->_mp_alloc = m3[0]->_mp_alloc = (int)(k4 * 8 / sizeof(mp_limb_t)); // to ensure 32 bits portability
	m0[0]->_mp_size = m1[0]->_mp_size = m2[0]->_mp_size = m3[0]->_mp_size = (int)(k4 * 8 / sizeof(mp_limb_t));	 // to ensure 32 bits portability
	for (uint64_t i = 0; i < N; i++)
	{
		for (uint64_t l = 0; l < k; l++)
		{
			uint64_t tmp = (uint64_t)in[l + i * stride];
			uint16_t *tptr = reinterpret_cast<uint16_t *>(&tmp);
			A0[l] = tptr[0];
			A1[l + 1] = tptr[1];
			A2[l + 2] = tptr[2];
			A3[l + 3] = tptr[3];
		}
		// see A0,A1,A2,A3 as a the gmp integers a0,a1,a2,a3
		m0[0]->_mp_d = reinterpret_cast<mp_limb_t *>(&A0[0]);
		m1[0]->_mp_d = reinterpret_cast<mp_limb_t *>(&A1[0]);
		m2[0]->_mp_d = reinterpret_cast<mp_limb_t *>(&A2[0]);
		m3[0]->_mp_d = reinterpret_cast<mp_limb_t *>(&A3[0]);
		res = a0;
		res += a1;
		res += a2;
		res += a3;
		out[i] = res;
	}
	m0[0]->_mp_d = m0_d;
	m1[0]->_mp_d = m1_d;
	m2[0]->_mp_d = m2_d;
	m3[0]->_mp_d = m3_d;
	m0[0]->_mp_alloc = m1[0]->_mp_alloc = m2[0]->_mp_alloc = m3[0]->_mp_alloc = 1;
	m0[0]->_mp_size = m1[0]->_mp_size = m2[0]->_mp_size = m3[0]->_mp_size = 0;
}

void intmatmul(uint64_t m, uint64_t n, uint64_t k, Integer *c, const Integer *a, const Integer *b)
{
	LinBox::Timer chrono;
	chrono.start();
	uint64_t dc; // real degree of cx
	uint64_t pts, lpts, num_primes, bits_a, bits_b, bits_c;

	bits_a = a[0].bitsize();
	for (uint64_t i = 1; i < m * k; i++)
		bits_a = std::max(bits_a, (uint64_t)a[i].bitsize());
	bits_b = b[0].bitsize();
	for (uint64_t i = 1; i < k * n; i++)
		bits_b = std::max(bits_b, (uint64_t)b[i].bitsize());

	bits_c = bits_a + bits_b + Integer(k).bitsize();
#ifdef VERBOSE
	cout << "nbr bit of a: " << bits_a << endl;
	cout << "nbr bit of b: " << bits_b << endl;
	cout << "nbr bit of result: " << bits_c << endl;
#endif
	vector<Integer> primes;

	kronecker_init(dc, primes, k, bits_c);

	num_primes = primes.size();
	lpts = Integer(dc).bitsize();
	pts = uint64_t(1) << lpts;
	std::vector<double> basis(primes.size());
	std::copy(primes.begin(), primes.end(), basis.begin());
	FFPACK::rns_double RNS(basis);
	Integer primesprod(1);
	for (uint64_t i = 0; i < num_primes; i++)
		primesprod *= Integer(primes[i]);
#ifdef VERBOSE
	cout << "\n--------------- ";
	cout << "  " << num_primes << "- primes FFT " << endl;
	cout << "              primes : ";
	for (uint64_t i = 0; i < num_primes; i++)
		cout << primes[i] << " ";
	cout << endl;
	cout << "   kronecker degree  : " << dc << endl;
	cout << "  prime prod bitsize : " << primesprod.bitsize() << endl;
	cout << "----------------\n\n";
#endif

	Givaro::Modular<Integer> FM(primesprod);

	// 2^16-adic decomposition of a and b over double
	std::vector<double> ax(m * k * pts);
	std::vector<double> bx(k * n * pts);
	std::vector<MatrixP *> cx_i(num_primes);
#ifdef PROFILE
	chrono.stop();
	cout << "  -> Kron:  init=" << chrono.usertime() << endl;
	chrono.start();
#endif
	for (uint64_t i = 0; i < m * k; i++)
		to_kronecker_16(ax.data() + i * pts, a[i], pts);
	for (uint64_t i = 0; i < k * n; i++)
		to_kronecker_16(bx.data() + i * pts, b[i], pts);
#ifdef PROFILE
	chrono.stop();
	cout << "  -> Kron:  tok=" << chrono.usertime() << endl;
	chrono.start();
#endif
	// pointwise multiplication in RNS
	std::vector<ModField> f(num_primes, ModField(2));
	for (uint64_t l = 0; l < num_primes; l++)
		f[l] = ModField(basis[l]);
#ifdef PROFILE
	LinBox::Timer mul;
	double tmul = 0.;
#endif
	MatrixP ai(f[0], m, k, pts);
	MatrixP bi(f[0], k, n, pts);

	for (uint64_t i = 0; i < num_primes; i++)
	{

		LinBox::PolynomialMatrixFFTPrimeMulDomain<ModField> fftdomain(f[i]);
		//MatrixP ai(f[i],m,k,pts);
		//MatrixP bi(f[i],k,n,pts);
		ai.changeField(f[i]);
		bi.changeField(f[i]);
		FFLAS::fassign(f[i], m * k * pts, ax.data(), 1, ai.getWritePointer(), 1);
		FFLAS::fassign(f[i], k * n * pts, bx.data(), 1, bi.getWritePointer(), 1);

		cx_i[i] = new MatrixP(f[i], m, n, pts);
		//cout<<ai<<endl;
		//cout<<bi<<endl;
#ifdef PROFILE
		mul.clear();
		mul.start();
#endif
		fftdomain.mul_fft(lpts, *cx_i[i], ai, bi);
#ifdef PROFILE
		mul.stop();
		tmul += mul.usertime();
#endif
		//cout<<*cx_i[i]<<endl;
	}
#ifdef PROFILE
	chrono.stop();
	cout << "  -> Kron:  mul=" << chrono.usertime() << " (" << tmul << ")" << endl;
	chrono.start();
#endif
	/*
	  std::vector<Integer>  cx   (m*n*dc);
	  // construct contiguous storage for c_i
	  uint64_t ncx= m*n*dc;
	  std::vector<double> cx_rns(ncx*num_primes);
	  for (uint64_t l=0;l<num_primes;l++){
	  for (uint64_t i=0;i<m*n;i++)
	  for (uint64_t j=0;j<dc;j++)
	  cx_rns[l*ncx + (j+i*dc)]= cx_i[l]->get(i,j);
	  delete cx_i[l];
	  }
	  // convert from RNS	       
	  RNS.convert(1,ncx,Integer(0),cx.data(),ncx, cx_rns.data(), ncx);	
	*/

	Givaro::ZRing<uint64_t> Z64;
	Givaro::ZRing<double> Z53;
	// reconstruct the result with MRS
	double alpha;
	!
		// compute cx_i[1]=(cx_i[1]-cx_i[0]) m0^(-1) mod m1
		f[1].init(alpha, basis[0]);
	f[1].invin(alpha);
	FFLAS::fsubin(f[1], m * n * pts, cx_i[0]->getPointer(), 1, cx_i[1]->getWritePointer(), 1);
	FFLAS::fscalin(f[1], m * n * pts, alpha, cx_i[1]->getWritePointer(), 1);
	if (num_primes == 3)
	{
		// compute ((cx_i[2]-cx_i[1]) m0^(-1) - cx_i[1]) m1^(-1) mod m2
		f[2].init(alpha, basis[0]);
		f[2].invin(alpha);
		FFLAS::fsubin(f[2], m * n * pts, cx_i[0]->getPointer(), 1, cx_i[2]->getWritePointer(), 1);
		FFLAS::fscalin(f[2], m * n * pts, alpha, cx_i[2]->getWritePointer(), 1);
		f[2].init(alpha, basis[1]);
		f[2].invin(alpha);
		FFLAS::fsubin(f[2], m * n * pts, cx_i[1]->getPointer(), 1, cx_i[2]->getWritePointer(), 1);
		FFLAS::fscalin(f[2], m * n * pts, alpha, cx_i[2]->getWritePointer(), 1);
	}

	// compute cx= cx_i[0] + m1(cx_i[1];
	std::vector<uint64_t> cx(m * n * pts);
	//FFLAS::faxpby(Z53,m*n*pts,1.0,cx_i[0]->getPointer(),1,basis[1],cx_i[1]->getWritePointer(),1);
	FFLAS::fscalin(Z53, m * n * pts, basis[0], cx_i[1]->getWritePointer(), 1);
	FFLAS::faddin(Z53, m * n * pts, cx_i[0]->getPointer(), 1, cx_i[1]->getWritePointer(), 1);
	FFLAS::finit(Z64, m * n * pts, cx_i[1]->getPointer(), 1, cx.data(), 1);

	if (num_primes == 3)
	{
		// compute cx=cx+ m1m2 cx_i[2]
		std::vector<uint64_t> cx2(m * n * pts);
		uint64_t m1m2 = uint64_t(basis[1]) * uint64_t(basis[2]);
		FFLAS::finit(Z64, m * n * pts, cx_i[2]->getPointer(), 1, cx2.data(), 1);
		FFLAS::faxpy(Z64, m * n * pts, m1m2, cx2.data(), 1, cx.data(), 1);
	}
	for (uint64_t i = 0; i < num_primes; i++)
		delete cx_i[i];
#ifdef PROFILE
	chrono.stop();
	cout << "  -> Kron:  rns=" << chrono.usertime() << endl;
	chrono.start();
#endif
	//for(uint64_t i=0;i<m*n;i++)
	//	from_kronecker_16(c[i],cx.data()+i*pts, dc);
	Integer bound = Integer(n) * Integer(dc) * Integer((uint64_t)1 << 32);

	//from_kronecker_16(m*n, c, cx.data(), dc, pts, bound.bitsize());
	from_kronecker_16(m * n, c, cx.data(), dc, pts, primesprod.bitsize());
#ifdef PROFILE
	chrono.stop();
	cout << "  -> Kron:  frk=" << chrono.usertime() << endl;
	chrono.start();
#endif
	// Givaro::ZRing<Givaro::Integer> Z;
	// write_field(Z,std::cout,a,m,k,1);
	// write_field(Z,std::cout,b,k,n,1);
	// write_field(Z,std::cout,c,m,k,1);
}

//} // end of namespace LinBox

#endif // KRONECKERFFT_H
