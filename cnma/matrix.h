/*
 * matrix.h
 *
 *  Created on: Sep 15, 2017
 *      Author: lidiia
 */

#ifndef MATRIX_H_
#define MATRIX_H_

#include <stdio.h>
#include <stdlib.h>
#include <gmp.h>
#include <math.h>

#include <sys/types.h>
#include <sys/stat.h>
#include <fcntl.h>
#include <inttypes.h>
#include <time.h>

struct time_mul{
	char *method;
	float seconds;
};

struct time_convert{
	char *method;
	float seconds;
};

#endif /* MATRIX_H_ */
