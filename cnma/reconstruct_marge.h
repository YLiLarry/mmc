#include "matrix.h"

void garner(mpz_t a, int N, const mpz_t r[], const mpz_t m[], const mpz_t Mi[]);
void prod(mpz_t M, mpz_t Mi, mpz_t* m, int N);
void precompute_Mi(mpz_t Mi[], const mpz_t m[], const size_t N);
