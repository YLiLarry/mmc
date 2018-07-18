#if !defined(H_PARGE_HUM)
#define H_PARGE_HUM

#include <gmp++/gmp++.h>
#include "../nocopy_integer.h"

// void dc_reduce_plus(LInteger& a, const LInteger& mod);
void dc_reduce_plus(mpz_t a, unsigned long int n);

#endif // H_PARGE_HUM
