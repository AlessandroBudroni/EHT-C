/*  This file is part of EHT-C.
 *
 *   name file: matrix.h
 *   Author: Alessandro Budroni
 *
 *  EHT-C is free software: you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation, either version 3 of the License, or
 *  (at your option) any later version.
 *
 *  EHT-C is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public License
 *  along with Nome-Programma.  If not, see <http://www.gnu.org/licenses/>
 */


#ifndef MATRIX_H
#define MATRIX_H

#include "config.h"
#include "utils.h"

/* sample random square NxN matrix */
void random_matrix(FP A[][N]);

/* test utilities */
int vector_equal(u16 n, FP A[n], FP B[n]);
// int matrix_equal(u16 n_rows, u16 n_cols, FP A[][n_cols], FP B[][n_cols]);

/* used in keygeneration and decryption */
int invert_matrix(FP inverseC[][N], FP C[][N]);
void fast_matrix_mul_times_vector(FP s[N], FP B[][N], FP b[N]);

/* stdout utilities */
// void print_matrix(u16 n_rows, u16 n_cols, FP Matrix[][n_cols]);
// void print_vector(u16 lenght, FP Vector[lenght]);

#endif