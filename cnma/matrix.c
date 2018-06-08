/*
 * matrix.c
 *
 *  Created on: Sep 15, 2017
 *      Author: lidiia
 */

#include "matrix.h"

void init_matrix(int n, mpz_t (*m)[n]);
void build_matrix(int n, mpz_t (*m)[n]);
void mul(int n, mpz_t (*m)[n], mpz_t (*x)[n], mpz_t (*c)[n]);
void timing(int iter, struct time_mul *s_mul, struct time_convert *s_conv);
int compare(int n, mpz_t (*a)[n], mpz_t (*b)[n]);
void bits(mpz_t r, mpz_t a, unsigned long int n, mp_bitcnt_t b);

void init_matrix(int n, mpz_t (*m)[n]){
	time_t t;
	int bit_n = 1024*8;

	for (int i = 0; i < n; i++) {
		for(int j = 0; j < n; j++){
		    mpz_init2(m[i][j], bit_n);
		}
	}
}

void build_matrix(int n, mpz_t (*m)[n]){
	mpz_t r_num;
	unsigned long int seed;
	gmp_randstate_t state;
	time_t t;

	srand(time(&t));
	for(int i = 0; i < n; i++){
		for(int j = 0; j < n; j++){
			seed = rand() % 50000;

			gmp_randinit_default (state);
			gmp_randseed_ui(state, seed);

			mpz_init(r_num);

			mpz_urandomb(r_num,state,50);
			mpz_set_ui(m[i][j], mpz_get_ui(r_num));
		}
	}
	gmp_randclear(state);
	mpz_clear(r_num);

}

void mul(int n, mpz_t (*m)[n], mpz_t (*x)[n], mpz_t (*c)[n]){
	 int i,j,k;

	 for (i = 0; i<n; i++){
		 for(j = 0; j < n; j++){
			 for (k = 0; k < n; k++){
				mpz_addmul(c[i][j], m[i][k], x[k][j]);
			 }
		 }
	 }
}

void timing(int iter, struct time_mul *s_mul, struct time_convert *s_conv){
	printf("Multiplication Times:\n");	
	printf("%30s%15s\n", "Method", "Time (s)");
	for (int i = 1; i <= iter; i++){
		printf("%30s%15f\n", (*(s_mul+i-1)).method, (*(s_mul+i-1)).seconds);
	}
	printf("\nConversion Times\n");
	for (int i = 1; i <= iter; i++){
		printf("%30s%15f\n", (*(s_conv+i-1)).method, (*(s_conv+i-1)).seconds);
	}
}


int compare(int n, mpz_t (*a)[n], mpz_t (*b)[n]){
	int check;
	for(int i = 0; i < n; i++){
		for(int j = 0; j < n; j++){
            check =  mpz_cmp(a[i][j], b[i][j]);
			if(check != 0) return 0;
		}
	}
	return 1;
}


void bits(mpz_t r, mpz_t a, unsigned long int n, mp_bitcnt_t b){
	mpz_tdiv_r_2exp(r, a, b);
	mpz_tdiv_q_2exp(r, r, n);
}


/*
long int generate_element(){
	mpz_t r_num;
	unsigned long int seed;
	gmp_randstate_t state;

	int devr = open("/dev/random", O_RDONLY);
	read(devr, &seed, 8);
	close(devr);

	gmp_randinit_default (state);
	gmp_randseed_ui(state, seed);

	mpz_init(r_num);

	mpz_urandomb(r_num,state,50);
//	gmp_printf("%Zd", r_num);
	long int n = mpz_get_ui (r_num);

	gmp_randclear(state);
	mpz_clear(r_num);

	return n;
}
*/

