/*  This file is part of EHT-C.
 *
 *   name file: ciphersuite.h
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


#ifndef CIPHERSUITE_H
#define CIPHERSUITE_H

#include "config.h"
#include "utils.h"

#ifdef FULL_STACK

typedef struct
{
    // FP RepC[M][LAMBDA2][2];  // CA = TB mod Q
    u16 P1[M];
    u16 P2[M];
    FP T[M];
    FP inverseB[N][N];
    double pce[Q];  // distribution associated to the row-vectors of C
    FP BCode[KDIM][N];
} privateKey;

typedef struct
{
    FP A[N][M]; // in transposed form
} publicKey;

typedef struct
{
    FP y[M];  // y = As + e mod q
} cipherText;

#else

typedef struct
{
    u16 *P1;  // CA = TB mod Q
    u16 *P2;
    FP *T;
    matrix *inverseB;
    double *pce;  // distribution associated to the row-vectors of C
    matrix *BCode;
} privateKey;

typedef struct
{
    matrix *A;
} publicKey;

typedef struct
{
    FP *y;  // y = As + e mod q
} cipherText;

#endif /* FULL_STACK */

#endif