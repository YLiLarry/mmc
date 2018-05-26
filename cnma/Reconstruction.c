#include "matrix.h"

void garner(mpz_t a, int N, mpz_t *r, mpz_t *m, mpz_t M, mpz_t Mi);
void prod(mpz_t M, mpz_t Mi, mpz_t *m, int N);

void garner(mpz_t a, int N, mpz_t *r, mpz_t *m, mpz_t M, mpz_t Mi){
	mpz_t arr[N];
	for(int i = 0; i < N; i++) 
		mpz_init(arr[i]);
	mpz_set(arr[0], r[0]);
	mpz_t t;
	mpz_init(t);
	for(int i = 1; i <= N-1; i++){
		mpz_set(t, arr[i-1]);
		for(int j = i-2; j >= 0; j--){
			mpz_mul(t, t, m[j]);
			mpz_add(t, t, arr[j]);
		}
		mpz_sub(t, r[i], t);
		mpz_t temp; 
		mpz_init(temp);
		mpz_mod(temp, Mi, m[i]);
		mpz_mul(arr[i], t, temp);
	}
	mpz_set(a, arr[N-1]);
	for(int i = N-1; i >= 0; i--){
		mpz_mul(a, a, m[i]);
		mpz_add(a, a, arr[i]);
	}		
}

void prod(mpz_t M, mpz_t Mi, mpz_t *m, int N){
	mpz_set_ui(M, 1);
	mpz_t plchldr;
	mpz_init(plchldr);
	for(int i = 1; i <= N-1; i++){
		mpz_mul(M, M, m[i-1]);
		mpz_gcdext(plchldr, Mi, plchldr, M, m[i]);
	}
}
