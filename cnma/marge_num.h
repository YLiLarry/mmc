#if !defined(H_MARGE_HUM)
#define H_MARGE_HUM

namespace CNMA {
void dc_reduce_minus(mpz_t a, unsigned long int n);
void minadd(mpz_t res, mpz_t a, mpz_t b, int n);
void minsub(mpz_t res, mpz_t a, mpz_t b, int n);
void minmul(mpz_t res, mpz_t a, mpz_t b, int n);
void get_mod(mpz_t m, int n);
void bits(mpz_t r, mpz_t a, unsigned long int n, mp_bitcnt_t b);
}
#endif // H_MARGE_HUM
