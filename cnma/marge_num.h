#if !defined(H_MARGE_HUM)
#define H_MARGE_HUM

#include <linbox/integer.h>
#include "../nocopy_integer.h"

using namespace LinBox;

// void dc_reduce_minus(LInteger& a, const LInteger& mod);
void dc_reduce_minus(mpz_t a, unsigned long int n);

#endif  // H_MARGE_HUM
