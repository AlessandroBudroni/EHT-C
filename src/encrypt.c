/*  This file is part of EHT-C.
 *
 *   name file: encrypt.c
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


#include <math.h>

#include "encrypt.h"

static FP CDF_TABLE[(int)CDF_TABLE_LEN];

/* Precompute CDF table */
void precompute_cdf_table()
{

    int len_x = 16;
    double sum;

    CDF_TABLE[0] = (1<<(len_x-1))*normal_PDF(0, 0, SIGMA)-1;
    for (int z = 1; z < CDF_TABLE_LEN; z++)
    {
        sum = 0;
        for (int i = 1; i <= z; i++)
            sum += normal_PDF(i, 0, SIGMA);
        CDF_TABLE[z] = CDF_TABLE[0]+(1<<len_x)*sum;
    }

}

/* Generate random secret - NOT SECURE
	generate N-KDIM entries at random, then encode.
*/
void generate_secret(FP Secret[N], FP G[][N-KDIM])
{

    // get seed for RNG
    int seed = get_seed();

    // initialize RNG
    srand((unsigned) seed);

    for (u16 row = 0; row < N-KDIM; row++)
        Secret[row] = (FP)(rand()%Q);

    encode(Secret, G);
}

/* Adapted from https://github.com/Microsoft/PQCrypto-LWEKE, under licence MIT. NOT CONSTANT TIME. */
void sample_error(FP *s, const size_t n)
{
    // Fills vector s with n samples from the noise distribution which requires 16 bits to sample.
    // The distribution is specified by its CDF.
    // Input: The input is overwritten by the output.

    // get seed for RNG
    int seed = get_seed();

    // initialize RNG
    srand((unsigned) seed);

    for (u16 row = 0; row < n; row++)
        s[row] = (FP)(rand());

    unsigned int i, j;

    for (i = 0; i < n; ++i)
    {
        FP sample = 0;
        FP prnd = s[i] >> 1;    // Drop the least significant bit
        FP sign = s[i] & 0x1;    // Pick the least significant bit

        // No need to compare with the last value.
        for (j = 0; j < (unsigned int)(CDF_TABLE_LEN - 1); j++)
        {
            // Constant time comparison: 1 if CDF_TABLE[j] < s, 0 otherwise. Uses the fact that CDF_TABLE[j] and s fit in 15 bits.
            sample += (FP)(CDF_TABLE[j] - prnd) >> 15;
        }
        // Assuming that sign is either 0 or 1, flips sample iff sign = 1. WARNING: NOT CONSTANT TIME
        s[i] = sign ? sample : Q-sample;
    }
}

/* Faster encryption y = Ax+e. Generate error on the fly using method adapted from: http://c-faq.com/lib/gaussian.html. */
void mul_times_secret_plus_error(FP y[M], FP A[][N], FP x[N], FP e[M])
{

    FP tmp_entry[4];

    for (u16 row = 0; row < M; row+=4)
    {
        tmp_entry[0] = 0;
        tmp_entry[1] = 0;
        tmp_entry[2] = 0;
        tmp_entry[3] = 0;
        for (u16 k = 0; k < N; k+=8)
        {
            tmp_entry[0] = (tmp_entry[0] + A[row  ][k]*x[k  ] + A[row  ][k+1]*x[k+1] + A[row  ][k+2]*x[k+2] + A[row  ][k+3]*x[k+3] + A[row  ][k+4]*x[k+4] + A[row  ][k+5]*x[k+5] + A[row  ][k+6]*x[k+6] + A[row  ][k+7]*x[k+7]) % Q;
            tmp_entry[1] = (tmp_entry[1] + A[row+1][k]*x[k  ] + A[row+1][k+1]*x[k+1] + A[row+1][k+2]*x[k+2] + A[row+1][k+3]*x[k+3] + A[row+1][k+4]*x[k+4] + A[row+1][k+5]*x[k+5] + A[row+1][k+6]*x[k+6] + A[row+1][k+7]*x[k+7]) % Q;
            tmp_entry[2] = (tmp_entry[2] + A[row+2][k]*x[k  ] + A[row+2][k+1]*x[k+1] + A[row+2][k+2]*x[k+2] + A[row+2][k+3]*x[k+3] + A[row+2][k+4]*x[k+4] + A[row+2][k+5]*x[k+5] + A[row+2][k+6]*x[k+6] + A[row+2][k+7]*x[k+7]) % Q;
            tmp_entry[3] = (tmp_entry[3] + A[row+3][k]*x[k  ] + A[row+3][k+1]*x[k+1] + A[row+3][k+2]*x[k+2] + A[row+3][k+3]*x[k+3] + A[row+3][k+4]*x[k+4] + A[row+3][k+5]*x[k+5] + A[row+3][k+6]*x[k+6] + A[row+3][k+7]*x[k+7]) % Q;
        }
        y[row    ] = (e[row    ] + tmp_entry[0])% Q;
        y[row + 1] = (e[row + 1] + tmp_entry[1])% Q;
        y[row + 2] = (e[row + 2] + tmp_entry[2])% Q;
        y[row + 3] = (e[row + 3] + tmp_entry[3])% Q;
    }

}

/* given the MasterSecret, sample error and perform encryption: y=Ax+e */
void EHT_encrypt(FP MasterSecret[N], cipherText *CipherText, publicKey *PublicKey)
{

    FP e[M];
    sample_error(e, M);
    mul_times_secret_plus_error(CipherText->y, PublicKey->A, MasterSecret, e);

}