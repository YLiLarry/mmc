/*
 * SargeNum.c
 * Uses moduli of forms 2- and 2+ and mpz_mod to reduce
 *
 *  Created on: Dec 20, 2017
 *      Author: Lidiia Zhitova
 */
#include "matrix.h"

void convert_mpz(int n, mpz_t (*)[n], mpz_t m); // converts entries of a matrix into modular representation
void modadd(mpz_t res, mpz_t u, mpz_t v, mpz_t m);
void modsub(mpz_t res, mpz_t u, mpz_t v, mpz_t m);
void modmul(mpz_t res, mpz_t u, mpz_t v, mpz_t m); // multiplies 2 numbers in modular representation
void mx_mod_mul(int n, mpz_t (*A)[n], mpz_t (*B)[n], mpz_t (*C)[n], mpz_t mod); // multiplies 2 matrices using modmul


void convert_mpz(int n, mpz_t (*matrix)[n], mpz_t m){
	for(int i = 0; i < n; i++){
		for(int j = 0; j < n; j++){
			mpz_mod(*(*(matrix + i) + j), *(*(matrix + i) + j), m);
		}
	}
}

void modadd(mpz_t res, mpz_t u, mpz_t v, mpz_t m){
	mpz_add(res, u, v); // (u mod m + v mod m)
	mpz_mod(res, res, m); // (u mod m + v mod m) mod m = u+v mod m
}

void modsub(mpz_t res, mpz_t u, mpz_t v, mpz_t m){
	mpz_sub(res, u, v); // (u mod m - v mod m)
	mpz_mod(res, res, m); // (u mod m - v mod m) mod m = u-v mod m
}

void modmul(mpz_t res, mpz_t u, mpz_t v, mpz_t m){
	mpz_mul(res, u, v); // (u mod m * v mod m)
	mpz_mod(res, res, m); // (u mod m * v mod m) mod m = u*v mod m
}

void mx_mod_mul(int n, mpz_t (*A)[n], mpz_t (*B)[n], mpz_t (*C)[n], mpz_t mod){
	 for (int i = 0; i<n; i++){
		 for(int j = 0; j < n; j++){
			 mpz_t sum;
			 mpz_init(sum);
			 mpz_set_ui(sum,0);
			 for (int k = 0; k < n; k++){
				 modmul(sum, A[i][k], B[k][j], mod);
				 modadd(C[i][j], sum, C[i][j], mod);
			 }
			 mpz_clear(sum);
		 }
	 }
}
