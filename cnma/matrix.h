/*
 * matrix.h
 *
 *  Created on: Sep 15, 2017
 *      Author: lidiia
 */

#ifndef MATRIX_H_
#define MATRIX_H_

#include <gmp.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>

#include <fcntl.h>
#include <inttypes.h>
#include <sys/stat.h>
#include <sys/types.h>
#include <time.h>

void init_matrix(int n, mpz_t (*m)[n]);
void build_matrix(int n, mpz_t (*m)[n]);
void mul(int n, mpz_t (*m)[n], mpz_t (*x)[n], mpz_t (*c)[n]);
void timing(int iter, struct time_mul *s_mul, struct time_convert *s_conv);
int compare(int n, mpz_t (*a)[n], mpz_t (*b)[n]);
void bits(mpz_t r, mpz_t a, unsigned long int n, mp_bitcnt_t b);

struct time_mul {
    char *method;
    float seconds;
};

struct time_convert {
    char *method;
    float seconds;
};

#endif /* MATRIX_H_ */
