
#include <NTL/vec_vec_long.h>
#include <NTL/vec_ZZ.h>
#include <NTL/lzz_p.h>
#include <NTL/ZZ_p.h>
#include <NTL/ZZ_pX.h>
#include <NTL/ZZX.h>
#include <NTL/FFT.h>
#include <gmp.h>

#include "BlasCRT.h"

NTL_CLIENT


static inline mp_limb_t * DATA(_ntl_gbigint p) { 
  return ((mp_limb_t *) (((long *) (p)) + 2)); 
}


/*-----------------------------------------------------------------------*/
/*-----------------------------------------------------------------------*/
/*                 to FFT rep stuff                                      */
/*-----------------------------------------------------------------------*/
/*-----------------------------------------------------------------------*/

// computes an n = 2^k point convolution.
// WAS: if deg(x) >= 2^k, then x is first reduced modulo X^n-1.
// NOW: error
void my_ToFFTRep(FFTRep& y, const ZZ_pX& x, long k, BlasCRT& blas_crt){

   const ZZ_pFFTInfoT *FFTInfo = ZZ_p::GetFFTInfo();
   if (k > FFTInfo->MaxRoot) 
      ResourceError("Polynomial too big for FFT");

   long n, i, j, m;
   long nprimes = FFTInfo->NumPrimes;

   long hi = deg(x);
   y.SetSize(k);
   n = 1L << k;
   m = max(hi + 1, 0);
   if (n < m){
     cerr << "FFT: degree too large" << endl;
     exit(-1);
   }

   double tt;

   tt = GetTime();
   blas_crt.reduce(y.tbl, x.rep);
   cout << "mod: " << GetTime()-tt << endl;

   tt = GetTime();
   for (i = 0; i < nprimes; i++) {
     long *yp = &y.tbl[i][0];
     for (j = m; j < n; j++) {
       yp[j] = 0;
     }
   }
   cout << "pad: " << GetTime()-tt << endl;

   tt = GetTime();
   for (i = 0; i < nprimes; i++) {
      long *yp = &y.tbl[i][0];
      FFTFwd(yp, yp, k, i);
   }
   cout << "FFT: " << GetTime()-tt << endl;

}

/*-----------------------------------------------------------------------*/
/*-----------------------------------------------------------------------*/
/*                   from FFT rep stuff                                  */
/*-----------------------------------------------------------------------*/
/*-----------------------------------------------------------------------*/

// converts from FFT-representation to coefficient representation
// y is destroyed
void my_FromFFTRep(ZZ_pX& x, FFTRep& y, long hi, BlasCRT& blas_crt){

   const ZZ_pFFTInfoT *FFTInfo = ZZ_p::GetFFTInfo();
   long nprimes = FFTInfo->NumPrimes;
   const long *u = FFTInfo->u.elts();
   const long *prime = FFTInfo->prime.elts();
   const mulmod_precon_t  *uqinv = FFTInfo->uqinv.elts();
   const double *prime_recip = FFTInfo->prime_recip.elts();

   if (FFTInfo->crt_struct.special()) {
     cout << "Special mode. I don't know what this is." << endl;
     exit(-1);
   }

   double tt;
   long k, n, i, j, l;
   vec_long ModularRepBuf;

   k = y.k;
   n = (1L << k);

   tt = GetTime();
   for (i = 0; i < nprimes; i++) {
      long *yp = &y.tbl[i][0];
      FFTRev1(yp, yp, k, i);
   }
   cout << "iFFT: " << GetTime()-tt << endl;

   hi = min(hi, n-1);
   l = hi+1;
   l = max(l, 0);
   x.rep.SetLength(l);

   tt = GetTime();
   long * qtab = new long[l];
   vec_vec_long aa;
   aa.SetLength(nprimes);
   for (long i = 0; i < nprimes; i++)
     aa[i].SetLength(l);

   for (j = 0; j < l; j++) {
      double yd = 0.0;
      for (i = 0; i < nprimes; i++) {
	long r = MulModPrecon(y.tbl[i][j], u[i], prime[i], uqinv[i]);
	aa[i][j] = r;
	yd += double(r)*prime_recip[i];
      }
      qtab[j] = long(yd + 0.5);
   }
   cout << "preprocess: " << GetTime()-tt << endl;

   tt = GetTime();
   vec_ZZ vZ;
   vZ.SetLength(l);
   blas_crt.CRT(vZ, aa);
   cout << "CRT: " << GetTime()-tt << endl;

   // there is some Montgomery stuff here
   tt = GetTime();
   for (long j = 0; j < l; j++){
      MulAddTo(vZ[j], FFTInfo->MinusMModP, qtab[j]);
      FFTInfo->reduce_struct.eval(x.rep[j].LoopHole(), vZ[j]);
   }
   cout << "Montgomery: " << GetTime()-tt << endl;
   
   x.normalize();
   delete[] qtab;
}

/*-----------------------------------------------------------------------*/
/*-----------------------------------------------------------------------*/
/*                        FFT mul stuff                                  */
/*-----------------------------------------------------------------------*/
/*-----------------------------------------------------------------------*/

void my_FFTMul(ZZ_pX& x, const ZZ_pX& a, const ZZ_pX& b){
  long k, d;
  
  if (IsZero(a) || IsZero(b)) {
    clear(x);
    return;
  }
  double tt = GetTime();
  BlasCRT blas_crt = BlasCRT();
  cout << "create: " << GetTime()-tt << endl;
  
  d = deg(a) + deg(b);
  k = NextPowerOfTwo(d+1);
  
  FFTRep R1(INIT_SIZE, k), R2(INIT_SIZE, k);
  
  my_ToFFTRep(R1, a, k, blas_crt);
  my_ToFFTRep(R2, b, k, blas_crt);
  mul(R1, R1, R2);
  my_FromFFTRep(x, R1, d, blas_crt);
}

void my_mul(ZZ_pX& c, const ZZ_pX& a, const ZZ_pX& b){
   if (IsZero(a) || IsZero(b)) {
      clear(c);
      return;
   }

   if (&a == &b) {
      sqr(c, a);
      return;
   }

   long k = ZZ_p::ModulusSize();
   long s = min(deg(a), deg(b)) + 1;

   if (s == 1 || (k == 1 && s < 40) || (k == 2 && s < 20) ||
                 (k == 3 && s < 12) || (k <= 5 && s < 8) ||
                 (k <= 12 && s < 4) )  {
      PlainMul(c, a, b);
   }
   else if (s < 16) {
      ZZX A, B, C;
      conv(A, a);
      conv(B, b);
      KarMul(C, A, B);
      conv(c, C);
   }
   else {
      long mbits;
      mbits = NumBits(ZZ_p::modulus());
      double rat = SSRatio(deg(a), mbits, deg(b), mbits);

      if ((k >= 53  && rat < 1.10) || 
	  (k >= 106 && rat < 1.30) || 
	  (k >= 212 && rat < 1.75)) {
	ZZX A, B, C;
	conv(A, a);
	conv(B, b);
	SSMul(C, A, B);
	conv(c, C);
      }
      else {
	my_FFTMul(c, a, b);
      }
   }
}


int main(){

  for (long j = 10000; j < 18000; j += 2000){

    ZZ p;
    p = RandomBits_ZZ(j);
    p = 2*p + 1;
    ZZ_p::init(p);

    cout << "bitsize " << j << endl;

    for (long i = 10000; i < 20000; i += 300){
      cout << "degree " << i << endl;
      ZZ_pX f = random_ZZ_pX(i);
      ZZ_pX g = random_ZZ_pX(i);
      ZZ_pX h;
      
      double t;
      t = GetTime();
      my_mul(h, f, g);
      cout << "mul " << GetTime()-t << endl;
      
      t = GetTime();
      ZZ_pX h2 = f*g;
      cout << "ntl " << GetTime()-t << endl;

      if (h != h2){
	cout << "error " << j << " " << i << endl;
	// cout << "f=" << f << endl;
	// cout << "g=" << g << endl;
	return -1;
      }

    }
  }
  return 0;
}
