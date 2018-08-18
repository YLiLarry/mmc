#include "matrix.h"

namespace CNMA {
void garner_parge_shift(mpz_t a, size_t N, const mpz_t r[], const mpz_t m[], uint64_t coef);

void garner_parge_shift_mixed(mpz_t a,                  // output
                              size_t N,                    // size of r[], bitsize[], m[], Mi[], work[]
                              const mpz_t r[],          // remainders
                              const uint64_t bitsize[], // bitsize of each moduli, ie, n as in 2^n+1 
                              const mpz_t m[],          // moduli: first is 2^n, second is 2^n+3, rest are 2^i+1
                              const mpz_t Mi[],         // precomputed Mi (see paper)
                              uint64_t coef,            // coefficient used to generate these parges
                              mpz_t work[]);            // a work workarray, caller is responsible for initializing and freeing this for efficiency reason
}