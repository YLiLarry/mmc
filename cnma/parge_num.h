#if !defined(H_PARGE_HUM)
#define H_PARGE_HUM

void dc_reduce_minus(mpz_t a, unsigned long int n);
void dc_reduce_plus(mpz_t a, long unsigned int n);
void plusadd(mpz_t res, mpz_t a, mpz_t b, int n);
void plussub(mpz_t res, mpz_t a, mpz_t b, int n);
void plusmul(mpz_t res, mpz_t a, mpz_t b, int n);
void get_mod_plus(mpz_t m, int n);
void bits(mpz_t r, mpz_t a, unsigned long int n, mp_bitcnt_t b);

#endif // H_PARGE_HUM
