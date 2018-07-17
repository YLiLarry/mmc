#if !defined(H_MARGE_HUM)
#define H_MARGE_HUM

#include <gmp++/gmp++.h>
#include "../nocopy_integer.h"

// void dc_reduce_minus(LInteger& a, const LInteger& mod);
void dc_reduce_minus(mpz_t a, unsigned long int n);

#endif // H_MARGE_HUM
