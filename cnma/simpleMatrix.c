/*
 ============================================================================
 Name        : simpleMatrix.c
 Author      : Lidiia Zhitova
 Version     :
 Copyright   : Your copyright notice
 Description :
 ============================================================================
 */

#include "matrix.h"

//long int generate_element();

int n;
long int d;
int output;
int a[0][0];

int main(int argc, char *argv[]) {
	char ans;
	printf("Enter testing mode? [y/n] ");
	ans = getchar();
	if (ans == 'y'){
		test();	
	}
	else{
		/*int N = 3;
		mpz_t P, Pi, m[N];
		mpz_init(P);
		mpz_init(Pi);
		for (int i = 0; i < N; i++){
			mpz_init(m[i]);
		}
		mpz_set_ui(m[0], 7);
		mpz_set_ui(m[1], 13);
		mpz_set_ui(m[2], 19);	
		prod(P, Pi, &m, N);
		gmp_printf("%Zd, %Zd\n", P, Pi);*/

		mpz_t r_num;
		unsigned long int seed;
		gmp_randstate_t state;
		long unsigned int n;
		
		printf("Enter n: ");
		scanf("%d", &n);	

		seed = rand() % 500000000;
		gmp_randinit_default (state);
		gmp_randseed_ui(state, seed);
		mpz_init(r_num);
		mpz_urandomb(r_num,state,200);
		gmp_printf("\nGenerated number: %Zd\n", r_num);
		
		mpz_t res, m;
		mpz_init(res);
		mpz_init(m);
		get_mod(m, n);

		clock_t start = clock();
		mpz_mod(res, r_num, m);
		clock_t end = clock();
		float sec = (float)(end - start) / CLOCKS_PER_SEC;
		gmp_printf("GMP reduction result (2n^-1): %Zd\n", res);
		printf("Time: %f\n", sec);
		printf("\n");
		
		mpz_t cp_num;
		mpz_init(cp_num);
		mpz_set(cp_num, r_num);
		start = clock();
		dc_reduce_minus(r_num, n);
		end = clock();
		sec = (float)(end - start) / CLOCKS_PER_SEC;
		gmp_printf("DC_reduce_minus result: %Zd\n", r_num); 
		printf("Time: %f\n", sec);
		printf("\n");

		get_mod_plus(m, n);
		start = clock();
		mpz_mod(res, cp_num, m);
		end = clock();
		sec = (float)(end - start) / CLOCKS_PER_SEC;
		gmp_printf("GMP reduction result(2^n+1): %Zd\n", res);
		printf("Time: %f\n", sec);
		printf("\n");

		start = clock();
		dc_reduce_plus(cp_num, n);
		end = clock();
		sec = (float)(end - start) / CLOCKS_PER_SEC;
		gmp_printf("DC_reduce_plus result: %Zd\n", cp_num); 
		printf("Time: %f\n", sec);
		printf("\n"); 
	
		/*if (argc != 2){
			fprintf(stderr, "USAGE: %s matrix_size\n", argv[0]);
	    		exit(EXIT_FAILURE);*/
		}


		//argv[1] = matrix dimension
		/*n = atoi(argv[1]);
		mpz_t A[n][n];
		mpz_t B[n][n];
		mpz_t C[n][n];
		mpz_t mod, b;

		mpz_init(mod);
		mpz_set_ui(mod, 2047);

		init_matrix(n, A);
		init_matrix(n, B);
		init_matrix(n, C);
		build_matrix(n, A);
		build_matrix(n, B);


		mpz_t reduced_a;
		mpz_init(reduced_a);
		mpz_set_ui(reduced_a, A[0][0]);
		printf("\n");
		printf("dc_reduce_minus function testing\n");
		gmp_printf("Number before reduction: %Zd\n", reduced_a);
		dc_reduce_minus(reduced_a,11);
		gmp_printf("After reduction: %Zd\n", reduced_a);
		printf("\n	***\n\n");

		int num_of_methods = 2;
		struct time_mul mul_sec[num_of_methods];
		struct time_convert conv_sec[num_of_methods];

		//AxB = C, then mod C
		clock_t start = clock();
		mul(n, A, B, C);
		clock_t end = clock();
		float sec = (float)(end - start) / CLOCKS_PER_SEC;
		mul_sec[0].method = "AxB = C, then mod C";
		mul_sec[0].seconds = sec;

		start = clock();
		convert_mpz(n, C, mod);
		end = clock();
		sec = (float)(end - start) / CLOCKS_PER_SEC;
		conv_sec[0].method = "AxB = C, then mod C";
		conv_sec[0].seconds = sec;

		//mod A, mod B, then AxB = D
		mpz_t D[n][n];
		init_matrix(n, D);
		start = clock();
		convert_mpz(n, A, mod);
		convert_mpz(n, B, mod);
		end = clock();
		sec = (float)(end - start) / CLOCKS_PER_SEC;
		conv_sec[1].method = "mod A, mod B, then AxB = D";
		conv_sec[1].seconds = sec;

		start = clock();
		mx_mod_mul(n, A, B, D, mod);
		end = clock();
		sec = (float)(end - start) / CLOCKS_PER_SEC;
		mul_sec[1].method = "mod A, mod B, then AxB = D";
		mul_sec[1].seconds = sec;

		timing(num_of_methods, mul_sec, conv_sec);
		int val = compare(n, C, D);
		if(val)printf("\nmatrices C and D are identical\n");
		else printf("\nmatrices C and D are not identical\n");*/


		/*printf("\nMatrix C - method 1: (AxB = C, then mod C)\n");
		for (int i = 0; i < n; i++) {
			for(int j = 0; j < n; j++){
				gmp_printf("%Zd  ", C[i][j]);
			}
			printf("\n");
		}

		printf("\nMatrix D - method 2: (mod A, mod B, then AxB = D)\n");
		for (int i = 0; i < n; i++) {
			for(int j = 0; j < n; j++){
				gmp_printf("%Zd  ", D[i][j]);
			}
			printf("\n");
		}
	}*/
}

