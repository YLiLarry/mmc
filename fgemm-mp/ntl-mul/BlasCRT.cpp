#include <math.h>
#include <gmp.h>
#include <NTL/ZZ_p.h>
#include <NTL/vec_ZZ_p.h>
#include <NTL/ZZ_pX.h>
#include <NTL/sp_arith.h>
#include <NTL/SmartPtr.h>
#include <NTL/vec_vec_long.h>

extern "C" {
#include <cblas.h>
}

#include "BlasCRT.h"

using namespace std;
NTL_CLIENT


#define SIZE(p) (((long *) (p))[1])

static inline mp_limb_t * DATA(_ntl_gbigint p) { 
  return ((mp_limb_t *) (((long *) (p)) + 2)); 
}

  
#define BITS_TO_LIMBS(n)  (((n) + 63) / 64)

void mpz_import (ZZ& z, long count, long sizeM, const unsigned long *data){

  z.SetSize(sizeM);
  mp_limb_t * zp = DATA(z.rep.rep);
  mp_limb_t limb, byte;
  unsigned char * dp = (unsigned char *) data;
  int lbits;
    
#define ACCUMULATE(N)                                   \
  do {							\
    limb |= (mp_limb_t) byte << lbits;			\
    lbits += (N);					\
    if (lbits >= 64)					\
      {							\
	*zp++ = limb;					\
	lbits -= 64;					\
	limb = byte >> ((N) - lbits);			\
      }							\
  } while (0)
  
  limb = 0;
  lbits = 0;
  for (long i = 0; i < count; i++){
    for (long j = 0; j < 2; j++){
      byte = *dp++;
      ACCUMULATE (8);
    }
    byte = *dp & 15;
    ACCUMULATE (4);
    dp += 6;
  }
  if (lbits != 0){
    *zp++ = limb;
  }

  zp = DATA(z.rep);
  long sz = sizeM;
  while(sz > 0 && zp[sz-1] == 0)
    sz--;

  SIZE(z.rep.rep) = sz;
}




long BlasCRT::split (double *data, long step, const ZZ& z){

  long zsize = z.size();
  if (zsize == 0)
    return 0;

  long count = ((NumBits(z)-1) / 20) + 1;
  mp_limb_t *zp = DATA(z.rep);

  unsigned long tmp;
  mp_limb_t limb;
  long c16 = count / 16;

  long k = 0, j = 0;
  for (long i = 0; i < c16; i++, j += 16*step, k += 5){
    limb = zp[k];
    data[j] = limb & 1048575;
    limb >>= 20;
    data[j+step] = limb & 1048575;
    limb >>= 20;
    data[j+2*step] = limb & 1048575;
    limb >>= 20;
    tmp = limb & 15;

    limb = zp[k+1];
    data[j+3*step] = ((limb & 65535) << 4) | tmp;
    limb >>= 16;
    data[j+4*step] = limb & 1048575;
    limb >>= 20;
    data[j+5*step] = limb & 1048575;
    limb >>= 20;
    tmp = limb & 255;
    
    limb = zp[k+2];
    data[j+6*step] = ((limb & 4095) << 8) | tmp;
    limb >>= 12;
    data[j+7*step] = limb & 1048575;
    limb >>= 20;
    data[j+8*step] = limb & 1048575;
    limb >>= 20;
    tmp = limb & 4095;
    
    limb = zp[k+3];
    data[j+9*step] = ((limb & 255) << 12) | tmp;
    limb >>= 8;
    data[j+10*step] = limb & 1048575;
    limb >>= 20;
    data[j+11*step] = limb & 1048575;
    limb >>= 20;
    tmp = limb & 65535;

    limb = zp[k+4];
    data[j+12*step] = ((limb & 15) << 16) | tmp;
    limb >>= 4;
    data[j+13*step] = limb & 1048575;
    limb >>= 20;
    data[j+14*step] = limb & 1048575;
    limb >>= 20;
    data[j+15*step] = limb & 1048575;
  }

  long jj = 16*c16;
  limb = (k >= zsize ? 0 : zp[k++]);
  if (jj >= count)
    return count;
  data[j] = limb & 1048575;
  limb >>= 20;
  if (jj+1 >= count)
    return count;
  data[j+1*step] = limb & 1048575;
  limb >>= 20;
  if (jj+2 >= count)
    return count;
  data[j+2*step] = limb & 1048575;
  limb >>= 20;
  tmp = limb & 15;
  
  limb = (k >= zsize ? 0 : zp[k++]);
  if (jj+3 >= count)
    return count;
  data[j+3*step] = ((limb & 65535) << 4) | tmp;
  limb >>= 16;
  if (jj+4 >= count)
    return count;
  data[j+4*step] = limb & 1048575;
  limb >>= 20;
  if (jj+5 >= count)
    return count;
  data[j+5*step] = limb & 1048575;
  limb >>= 20;
  tmp = limb & 255;
  
  limb = (k >= zsize ? 0 : zp[k++]);
  if (jj+6 >= count)
    return count;
  data[j+6*step] = ((limb & 4095) << 8) | tmp;
  limb >>= 12;
  if (jj+7 >= count)
    return count;
  data[j+7*step] = limb & 1048575;
  limb >>= 20;
  if (jj+8 >= count)
    return count;
  data[j+8*step] = limb & 1048575;
  limb >>= 20;
  tmp = limb & 4095;
  
  limb = (k >= zsize ? 0 : zp[k++]);
  if (jj+9 >= count)
    return count;
  data[j+9*step] = ((limb & 255) << 12) | tmp;
  limb >>= 8;
  if (jj+10 >= count)
    return count;
  data[j+10*step] = limb & 1048575;
  limb >>= 20;
  if (jj+11 >= count)
    return count;
  data[j+11*step] = limb & 1048575;
  limb >>= 20;
  tmp = limb & 65535;
    
  limb = (k >= zsize ? 0 : zp[k++]);
  if (jj+12 >= count)
    return count;
  data[j+12*step] = ((limb & 15) << 16) | tmp;
  limb >>= 4;
  if (jj+13 >= count)
    return count;
  data[j+13*step] = limb & 1048575;
  limb >>= 20;
  if (jj+14 >= count)
    return count;
  data[j+14*step] = limb & 1048575;
  if (jj+15 >= count)
    return count;
  limb >>= 20;
  data[j+15*step] = limb;

  return count;
}

long BlasCRT::split (double *data, long step, const ZZ_p& z){
  return split(data, step, z._ZZ_p__rep);
}

// we add the construction of another matrix using powers of 2^60
BlasCRT::BlasCRT(){
  size = NumBits(ZZ_p::modulus()); // in bits
  size = 1 + (size / 20);

  if ((size % 3) != 0) 
    size++;
  if ((size % 3) != 0) 
    size++;
  size60 = size / 3;
  numPrimes = ZZ_p::GetFFTInfo()->NumPrimes;

  // builds the matrix used for reduction in the second method
  reductionMatrix = new double[6*numPrimes*size60];

  long size_of_one_mat = size60 * numPrimes;
  double * reductionMatrix_0 = reductionMatrix;
  double * reductionMatrix_1 = reductionMatrix + 1*size_of_one_mat;
  double * reductionMatrix_2 = reductionMatrix + 2*size_of_one_mat;
  double * reductionMatrix_3 = reductionMatrix + 3*size_of_one_mat;
  double * reductionMatrix_4 = reductionMatrix + 4*size_of_one_mat;
  double * reductionMatrix_5 = reductionMatrix + 5*size_of_one_mat;
  
  // [a0*b0,a1*b1,a2*b2,(a0-a2)*(b0-b2),(a0+a1+a2)*(b0+b1+b2),(a0+a1-a2)*(b0+b1-b2)]
  for (long i = 0; i < numPrimes; i++){
    long p = ZZ_p::GetFFTInfo()->prime[i];
    mulmod_t inv_p = PrepMulMod(p);

    if (p > (1L<<60)){
      cout << "error in building reduction matrix: prime too large\n";
      exit(-1);
    }
    
    long tmp = 1;
    long two_sixty = (1L << 60) % p;
    const long mask = 1048575;

    for (long j = 0; j < size60; j++){
      long b0 = tmp & mask;
      long b1 = (tmp >> 20) & mask;
      long b2 = (tmp >> 40) & mask;

      reductionMatrix_0[i*size60 + j] = b0;
      reductionMatrix_1[i*size60 + j] = b1;
      reductionMatrix_2[i*size60 + j] = b2;
      reductionMatrix_3[i*size60 + j] = b0-b2;
      reductionMatrix_4[i*size60 + j] = b0+b1+b2;
      reductionMatrix_5[i*size60 + j] = b0+b1-b2;

      tmp = MulMod(tmp, two_sixty, p, inv_p); // todo: prepmulmod this
    }
  }

  // builds the matrix used for CRT
  prod_vec.SetLength(numPrimes);
  ZZ tmp = to_ZZ(1);
  for (long i = 0; i < numPrimes; i++)
    tmp *= ZZ_p::GetFFTInfo()->prime[i];
  // using a subproduct tree is usually not really necessary, it seems.

  long nbits = 0;
  for (long i = 0; i < numPrimes; i++){
    ZZ cof = tmp / ZZ_p::GetFFTInfo()->prime[i];
    ZZ_p::GetFFTInfo()->reduce_struct.adjust(cof);
    prod_vec[i] = cof;
    nbits = max(NumBits(cof), nbits);
  }

  sizeM = ceil((double)nbits / 20.0);
  sizeM = 3*ceil((double)sizeM / 3.0); // make it a multiple of 3
  sizeM += 6; // pad with zeros at the end

  double * tmpCoeffs = new double[sizeM];
  crtMatrix = new double[2*sizeM*numPrimes];

  long s3 = sizeM / 3;
  long size_one_mat = s3*numPrimes;
  double *crtMatrix0 = crtMatrix;
  double *crtMatrix1 = crtMatrix + size_one_mat;
  double *crtMatrix2 = crtMatrix + 2*size_one_mat;
  double *crtMatrix3 = crtMatrix + 3*size_one_mat;
  double *crtMatrix4 = crtMatrix + 4*size_one_mat;
  double *crtMatrix5 = crtMatrix + 5*size_one_mat;

  for (long i = 0; i < numPrimes; i++){
    long count = split(tmpCoeffs, 1, prod_vec[i]);
    for (long j = count; j < sizeM; j++)
      tmpCoeffs[j] = 0;

    //diagD=diagonal_matrix([d0-d1+d2-d3, -d1+d2, -d1+d2-d3+d4, d1+d3, -d2+d3, d1-d3])
    for (long c = 0; c < s3; c++){

      double d0, d1, d2, d3, d4;
      if (c > 0){
	d0 = tmpCoeffs[3*c-2];
	d1 = tmpCoeffs[3*c-1];
      }
      else {
	d0 = 0;
	d1 = 0;
      }
      d2 = tmpCoeffs[3*c];
      d3 = tmpCoeffs[3*c+1];
      d4 = tmpCoeffs[3*c+2];

      crtMatrix0[i + c*numPrimes] = d0-d1+d2-d3;
      crtMatrix1[i + c*numPrimes] = -d1+d2;
      crtMatrix2[i + c*numPrimes] = -d1+d2-d3+d4;
      crtMatrix3[i + c*numPrimes] = d1+d3;
      crtMatrix4[i + c*numPrimes] = -d2+d3;
      crtMatrix5[i + c*numPrimes] = d1-d3;
    }

  }

  delete[] tmpCoeffs;
}


BlasCRT::~BlasCRT(){
  delete[] reductionMatrix;
  delete[] crtMatrix;
}

// reduces a vector of ZZ_p modulo the current primes
// second version, using faster pol-mul
void BlasCRT::reduce(Unique2DArray<long>& a, const vec_ZZ_p& coeffs) {

  double tt;

  tt = GetTime();
  long len = coeffs.length();
  double * mat_coeffs = new double[size60 * len * 6]; // was: size * len = size60 * 3 * len
  double * mat_res = new double[numPrimes * len * 6]; // was: 3 * numPrimes * len

  long size_of_one_mat_coeffs = size60 * len;
  double * mat_coeffs_0 = mat_coeffs;
  double * mat_coeffs_1 = mat_coeffs + 1*size_of_one_mat_coeffs;
  double * mat_coeffs_2 = mat_coeffs + 2*size_of_one_mat_coeffs;
  double * mat_coeffs_3 = mat_coeffs + 3*size_of_one_mat_coeffs;
  double * mat_coeffs_4 = mat_coeffs + 4*size_of_one_mat_coeffs;
  double * mat_coeffs_5 = mat_coeffs + 5*size_of_one_mat_coeffs;

  long size_of_one_mat_res = numPrimes * len;
  double * mat_res_0 = mat_res;
  double * mat_res_1 = mat_res + 1*size_of_one_mat_res;
  double * mat_res_2 = mat_res + 2*size_of_one_mat_res;
  double * mat_res_3 = mat_res + 3*size_of_one_mat_res;
  double * mat_res_4 = mat_res + 4*size_of_one_mat_res;
  double * mat_res_5 = mat_res + 5*size_of_one_mat_res;

  long size_of_one_mat_red = numPrimes * size60;
  double * reductionMatrix_0 = reductionMatrix;
  double * reductionMatrix_1 = reductionMatrix + 1*size_of_one_mat_red;
  double * reductionMatrix_2 = reductionMatrix + 2*size_of_one_mat_red;
  double * reductionMatrix_3 = reductionMatrix + 3*size_of_one_mat_red;
  double * reductionMatrix_4 = reductionMatrix + 4*size_of_one_mat_red;
  double * reductionMatrix_5 = reductionMatrix + 5*size_of_one_mat_red;
  cout << " alloc: " << GetTime()-tt << endl;


  tt = GetTime();
  double * buffer = new double[size];
  for (long i = 0; i < len; i++){
    long count = split(buffer, 1, coeffs[i]);
    for (long j = count; j < size; j++)
      buffer[j] = 0;
    for (long j = 0; j < size60; j++){
      long a0 = buffer[3*j];
      long a1 = buffer[3*j+1];
      long a2 = buffer[3*j+2];
      mat_coeffs_0[j*len + i] = a0;
      mat_coeffs_1[j*len + i] = a1;
      mat_coeffs_2[j*len + i] = a2;
    }
  }
  delete[] buffer;
  for (long i = 0; i < size60*len; i++){
    double a0 = mat_coeffs_0[i];
    double a1 = mat_coeffs_1[i];
    double a2 = mat_coeffs_2[i];
    double s = a0+a1;
    mat_coeffs_3[i] = a0-a2;
    mat_coeffs_4[i] = s+a2;
    mat_coeffs_5[i] = s-a2;
  }
  cout << " split: " << GetTime()-tt << endl;

 
  tt = GetTime();
  cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans, 
	      numPrimes, len, size60, 
	      1.0, 
	      reductionMatrix_0, size60, 
	      mat_coeffs_0, len, 
	      0.0, 
	      mat_res_0, len);
  cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans, 
	      numPrimes, len, size60, 
	      1.0, 
	      reductionMatrix_1, size60, 
	      mat_coeffs_1, len, 
	      0.0, 
	      mat_res_1, len);
  cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans, 
	      numPrimes, len, size60, 
	      1.0, 
	      reductionMatrix_2, size60, 
	      mat_coeffs_2, len, 
	      0.0, 
	      mat_res_2, len);
  cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans, 
	      numPrimes, len, size60, 
	      1.0, 
	      reductionMatrix_3, size60, 
	      mat_coeffs_3, len, 
	      0.0, 
	      mat_res_3, len);
  cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans, 
	      numPrimes, len, size60, 
	      1.0, 
	      reductionMatrix_4, size60, 
	      mat_coeffs_4, len, 
	      0.0, 
	      mat_res_4, len);
  cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans, 
	      numPrimes, len, size60, 
	      1.0, 
	      reductionMatrix_5, size60, 
	      mat_coeffs_5, len, 
	      0.0, 
	      mat_res_5, len);
  cout << " matrix product: " << GetTime()-tt << endl;


  tt = GetTime();
  const long two_fourty = 1L << 40;
  for (long j = 0; j < numPrimes; j++){
    long p = ZZ_p::GetFFTInfo()->prime[j];
    mulmod_t inv_p = PrepMulMod(p);
    mulmod_precon_t two_fourty_prec = PrepMulModPrecon(two_fourty, p, inv_p);
    long *aj = &a[j][0];
    
    for (long i = 0; i < len; i++){
      unsigned long c0 = mat_res_0[j*len + i];
      long p1 = mat_res_1[j*len + i];
      long c4 = mat_res_2[j*len + i];
      long p3 = mat_res_3[j*len + i];
      long p4 = mat_res_4[j*len + i];
      long p5 = mat_res_5[j*len + i];

      long t = c0+c4;
      long c2 = t-p3+p1;
      long s = p4-t-c2;

      unsigned long c1 = ((unsigned long) (s+p5-p3-p1)) >> 1;
      long c3 = s-c1;
      const unsigned long mask20 = (1L << 20) - 1;
      c0 = c0 + ((c1 & mask20) << 20);
      c2 = c2 + (c1 >> 20) + ((c3 & mask20) << 20);
      c4 = c4 + (c3 >> 20);
      long res = MulModPrecon(c4, two_fourty, p, two_fourty_prec);
      res = MulModPrecon(res+c2, two_fourty, p, two_fourty_prec);
      aj[i] = AddMod(res, c0, p);
    }
  }
  cout << " unsplit: " << GetTime()-tt << endl;

  tt = GetTime();
  delete[] mat_coeffs;
  delete[] mat_res;
  cout << " free: " << GetTime()-tt << endl;
}






void BlasCRT::CRT(vec_ZZ& coeffs, const vec_vec_long& a) {
  
  double tt;

  tt = GetTime();
  const long len = a[0].length();
  const long s3 = sizeM/3;

  long size_one_mat = s3*numPrimes;
  double *crtMatrix0 = crtMatrix;
  double *crtMatrix1 = crtMatrix + size_one_mat;
  double *crtMatrix2 = crtMatrix + 2*size_one_mat;
  double *crtMatrix3 = crtMatrix + 3*size_one_mat;
  double *crtMatrix4 = crtMatrix + 4*size_one_mat;
  double *crtMatrix5 = crtMatrix + 5*size_one_mat;

  double * mat_residues = new double[numPrimes*6*len];
  long size_of_one_mat1 = numPrimes*len;
  double * mat_residues0 = mat_residues;
  double * mat_residues1 = mat_residues + 1*size_of_one_mat1;
  double * mat_residues2 = mat_residues + 2*size_of_one_mat1;
  double * mat_residues3 = mat_residues + 3*size_of_one_mat1;
  double * mat_residues4 = mat_residues + 4*size_of_one_mat1;
  double * mat_residues5 = mat_residues + 5*size_of_one_mat1;

  double * mat_result = new double[sizeM*2*len];
  long size_of_one_mat2 = s3*len;
  double * mat_result0 = mat_result;
  double * mat_result1 = mat_result + 1*size_of_one_mat2;
  double * mat_result2 = mat_result + 2*size_of_one_mat2;
  double * mat_result3 = mat_result + 3*size_of_one_mat2;
  double * mat_result4 = mat_result + 4*size_of_one_mat2;
  double * mat_result5 = mat_result + 5*size_of_one_mat2;
  cout << "  alloc: " << GetTime()-tt << endl;


  tt = GetTime();
  for (long i = 0; i < numPrimes; i++){
    for (long j = 0; j < len; j++){
      long tmp = a[i][j];
      double a0 = tmp & 1048575;
      double a1 = (tmp >> 20) & 1048575;
      double a2 = (tmp >> 40) & 1048575;
      mat_residues0[i*len + j] = a2;
      mat_residues1[i*len + j] = a1;
      mat_residues2[i*len + j] = a0;
      double a2m0 = a2-a0;
      mat_residues3[i*len + j] = a0+a1+a2;
      mat_residues4[i*len + j] = a2m0;
      mat_residues5[i*len + j] = a2m0+a1;
    }
  }
  cout << "  setup: " << GetTime()-tt << endl;

  tt = GetTime();
  cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans, 
	      s3, len, numPrimes, 
	      1.0, 
	      crtMatrix0, numPrimes, 
	      mat_residues0, len, 
	      0.0, 
	      mat_result0, len);
  cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans, 
	      s3, len, numPrimes, 
	      1.0, 
	      crtMatrix1, numPrimes, 
	      mat_residues1, len, 
	      0.0, 
	      mat_result1, len);
  cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans, 
	      s3, len, numPrimes, 
	      1.0, 
	      crtMatrix2, numPrimes, 
	      mat_residues2, len, 
	      0.0, 
	      mat_result2, len);
  cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans, 
	      s3, len, numPrimes, 
	      1.0, 
	      crtMatrix3, numPrimes, 
	      mat_residues3, len, 
	      0.0, 
	      mat_result3, len);
  cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans, 
	      s3, len, numPrimes, 
	      1.0, 
	      crtMatrix4, numPrimes, 
	      mat_residues4, len, 
	      0.0, 
	      mat_result4, len);
  cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans, 
	      s3, len, numPrimes, 
	      1.0, 
	      crtMatrix5, numPrimes, 
	      mat_residues5, len, 
	      0.0, 
	      mat_result5, len);
  cout << "  matrix product: " << GetTime()-tt << endl;

  tt = GetTime();
  unsigned long * output = new unsigned long[sizeM];
  for (long j = 0; j < len; j++){
    unsigned long res = 0;
    for (long i = 0; i < s3; i++){
      long r0 = (long) mat_result0[i*len+j];
      long r1 = (long) mat_result1[i*len+j];
      long r2 = (long) mat_result2[i*len+j];
      long r3 = (long) mat_result3[i*len+j];
      long r4 = (long) mat_result4[i*len+j];
      long r5 = (long) mat_result5[i*len+j];
      
      long r35 = (r3+r5) >> 1;
      long rm35 = (r3-r5) >> 1;

      unsigned long tmp0 = r0+r4+r35;
      unsigned long tmp1 = r1+r35;
      unsigned long tmp2 = (r2+rm35)-r4;

      tmp0 += res;
      output[3*i] = tmp0 & 1048575;
      tmp1 += tmp0 >> 20;
      output[3*i+1] = tmp1 & 1048575;
      tmp2 += tmp1 >> 20;
      output[3*i+2] = tmp2 & 1048575;
      res = tmp2 >> 20;
    }
    mpz_import(coeffs[j], sizeM, ceil((double)sizeM * 20.0/64.0), output);
  }
  cout << "  combine: " << GetTime()-tt << endl;

  
  tt = GetTime();
  delete[] output;
  delete[] mat_residues;
  delete[] mat_result;
  cout << "  free: " << GetTime()-tt << endl;

}
