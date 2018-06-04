
#ifndef BLASCRT_H_
#define BLASCRT_H_

#include <NTL/ZZ_p.h>
#include <NTL/vector.h>
#include <NTL/vec_ZZ_p.h>
#include <NTL/vec_ZZ.h>
#include <NTL/ZZ_pX.h>
#include <NTL/SmartPtr.h>
#include <NTL/vec_vec_long.h>

NTL_CLIENT

using namespace std;

class BlasCRT{

  vec_ZZ prod_vec;
  double * reductionMatrix; // for words of 60 bits
  double * crtMatrix;
    
public:
  
  // size is in 20 bits, size60 in 60 bits
  long size, size60, sizeM, numPrimes;
  long split (double *data, long step, const ZZ& z);
  long split (double *data, long step, const ZZ_p& z);

  void reduce(Unique2DArray<long>& a, const vec_ZZ_p& coeffs);
  void CRT(vec_ZZ& coeffs, const vec_vec_long& a);

  BlasCRT();
  ~BlasCRT();

};


#endif /* BIGINTMAT_H_ */

