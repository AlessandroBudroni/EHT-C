/*  This file is part of EHT-C.
 *
 *   name file: keygen.c
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


#include "keygen.h"
#include "ciphersuite.h"
#include "matrix.h"

#include <immintrin.h>

#include <string.h>
#include <stdint.h>
#include <stdlib.h>
#include <math.h>

/* convolution formula for computing distirbution table
   compute distribution of e1+e2. store distribution in pce1
*/
void convolution_sum(double *pce1, double *pce2)
{

    double sum[Q];
    memset(sum, 0, Q*sizeof(double));
    for (int w = 0; w < Q; w++)
    {
        for (int a = 0; a < Q; a++)
            sum[w] += pce1[a]*pce2[(int)make_FP(w-a)];
    }

    for (int w = 0; w < Q; w++)
        pce1[w] = sum[w];
}

/* recursive formula to compute distribution fast
   compute distribution e1+e2+...+en, where ei~pe */
void recursive(u16 length, double *pce, double *pe)
{

    if(length == 1)
    {
        for (int w = 0; w < Q; w++)
            pce[w] = pe[w];
        return;
    }

    recursive(length/2, pce, pe);
    convolution_sum(pce, pce);

    if (length%2 == 0)
        return;
    else
    {
        convolution_sum(pce, pe);
        return;
    }

}

/* precompute distribution of e1+e2+...+elambda2, needed in decryption */
void precompute_distribuion(privateKey *PrivateKey)
{

    double sigma = SIGMA;

#ifdef PRINT_STAMP
    char stamp[200];
    sprintf(stamp, "Compute distribution table");
    time_stamp(stamp);
#endif

    // else, compute the distribution
    // error distribution table
    double pe[Q];
    memset(pe, 0, Q*sizeof(double));
    pe[0] = normal_PDF(0, 0, sigma);
    for(int i = 1; i<PRECISION; i++)
    {
        pe[i] = normal_PDF(i, 0, sigma);
        pe[Q-i] = pe[i];
    }

    memset(PrivateKey->pce, 0, Q*sizeof(double));
    recursive(LAMBDA2, PrivateKey->pce, pe);

    for (int i = 0; i < Q; ++i)
        PrivateKey->pce[i] = PrivateKey->pce[i]<0.00000000001 ? -9999999 : log(Q*PrivateKey->pce[i]);

#ifdef PRINT_STAMP
    sprintf(stamp, "Done");
    time_stamp(stamp);
#endif
}

// Key generation using Hadamard matrices

/* instead of permutation matrix, we store a list that represent the permutation matrix. This is used to speed up key generation and avoid large matrix multiplication */
void get_permutations(u16 P1[M], u16 P2[M])
{

    // seed for random generation
    int seed = get_seed();
    // initialize RNG
    srand((unsigned) seed);

    // Method: set non-zero entries using sampling without replacement
    int list1[M], list2[M];
    for (int i = 0; i < M; ++i)
    {
        list1[i] = i;
        list2[i] = i;
    }

    int rand_pos1, rand_value1, tmp_value1, rand_pos2, rand_value2, tmp_value2;

    for (int i = 0; i < M; ++i)
    {
        // get random from 0 to n-1 without replacement
        rand_pos1 = rand()%(M-i);
        rand_value1 = list1[rand_pos1];

        // swap values to remove from list
        tmp_value1 = list1[M-i-1];
        list1[M-i-1] = rand_value1;
        list1[rand_pos1] = tmp_value1;

        P1[i] = rand_value1;

        // get random from 0 to n-1 without replacement
        rand_pos2 = rand()%(M-i);
        rand_value2 = list2[rand_pos2];

        // swap values to remove from list
        tmp_value2 = list2[M-i-1];
        list2[M-i-1] = rand_value2;
        list2[rand_pos2] = tmp_value2;

        P2[i] = rand_value2;
    }
}

/* inplace row-permutation of matrix in its transposed form (so it is actually column permutation) */
void inplace_row_permutation_in_transposed_form(FP A[N][M], u16 P[M])
{

    FP tempv;
    uint8_t past_pos[M];
    memset(past_pos, 0, M*sizeof(uint8_t));
    int pos = 0, start = 0, dest = 0, processed = 0, point = 0;
    dest = P[pos];

    while(processed < M)
    {
        if(past_pos[dest] != 1)
        {
            for (int i = 0; i < N; ++i)
            {
                tempv = A[i][pos];
                A[i][pos] = A[i][dest];
                A[i][dest] = tempv;
            }
        }
        processed++;
        past_pos[pos] = 1;

        if (start == dest && processed < M)
        {
            while(past_pos[point] == 1 && point<M)
                point++;

            pos = point;
            start = pos;
        }
        else
            pos = dest;
        dest = P[pos];
    }
}

/* inplace row-permutation of matrix in its transposed form (so it is actually column permutation) and multiply times lambda2^-1*/
void inplace_row_permutation_in_transposed_form_and_mul_ilambda2(FP A[N][M], u16 P[M])
{
	FP invlambda2 = inverse_mod((FP)LAMBDA2);
    FP tempv;
    uint8_t past_pos[M];
    memset(past_pos, 0, M*sizeof(uint8_t));
    int pos = 0, start = 0, dest = 0, processed = 0, point = 0;
    dest = P[pos];

    while(processed < M)
    {
        if(past_pos[dest] != 1)
        {
            for (int i = 0; i < N; ++i)
            {
                tempv = A[i][pos];
                A[i][pos] = A[i][dest];
                A[i][dest] = tempv;
                A[i][pos] = (A[i][pos]*invlambda2)%Q;
            }
        }
        else
        {
        	for (int i = 0; i < N; ++i)
        		 A[i][pos] = (A[i][pos]*invlambda2)%Q;
        }
        processed++;
        past_pos[pos] = 1;

        // for (int i = 0; i < N; ++i)
    	   //  A[i][pos] = (A[i][pos]*invlambda2)%Q;

        if (start == dest && processed < M)
        {
            while(past_pos[point] == 1 && point<M)
                point++;

            pos = point;
            start = pos;
        }
        else
            pos = dest;
        dest = P[pos];
    }
}

/* row permutation of matrix in trasposed form (so it's a column permutation) */
void row_permutation_in_transposed_form(FP A[N][M], FP tempA[N][M], FP P[M])
{

    FP invlambda2 = inverse_mod((FP)LAMBDA2);

    for (int row = 0; row < M; row++)
        for (int col = 0; col < N; col++)
            A[col][row] = (tempA[col][P[row]]*invlambda2)%Q;
}

/* tempA is given in transpose form. Apply P1, multiply times lambda2^-1, return in the right form (not transpose) */
void row_permutation_and_transpose(FP A[M][N], FP tempA[N][M], FP P[M])
{

    FP invlambda2 = inverse_mod((FP)LAMBDA2);

    for (int row = 0; row < M; row++)
        for (int col = 0; col < N; col++)
            A[row][col] = (tempA[col][P[row]]*invlambda2)%Q;
}


/* generate public/private keypair with the Hadamard matrix-based approach */
void EHT_keygen(privateKey *PrivateKey, publicKey *PublicKey, FP H[][N])
{

    // generate private B and its inverse
    FP B[N][N];
    do
    {
        random_matrix(B);
    }
    while(!invert_matrix(PrivateKey->inverseB, B) );

    // generate private T TODO make this random properly
    int seed = get_seed();
    srand((unsigned) seed);
    for (u16 i = 0; i < M; i++)
        PrivateKey->T[i] = (i%K == 0) ? 1 : make_FP((FP)rand()%Q); // set 1 to the first entry of each chunk

    // generate random permutations P1, P2
    get_permutations(PrivateKey->P1, PrivateKey->P2);

    // get inverse of permutations, needed when computing A
    u16 iP2[M], iP1[M];
    for (int i = 0; i < M; ++i)
    {
        iP1[PrivateKey->P1[i]] = i;
        iP2[PrivateKey->P2[i]] = i;
    }

    // multiply C^-1 times TB to get A
    // multiply T times B and transpose automatically. This is needed to use the FWHT later.
    // FP tempA[N][M];
    for (int row = 0; row < M; row++)
        for (int col = 0; col < N; col++)
            PublicKey->A[col][row] = (PrivateKey->T[row]*B[row/K][col]) % Q;

    // permute first TB with the inverse of P2
    inplace_row_permutation_in_transposed_form(PublicKey->A, iP1);

    // apply FHWT to the transpose of TB already permuted by iP1
    for (int row = 0; row < N; row++)
        for (int col = 0; col < M/LAMBDA2; col++)
            FWHT(PublicKey->A[row]+LAMBDA2*col, LAMBDA2);

    // row permutation of tempA in its transpose form, then multiply by lambda2^-1
    inplace_row_permutation_in_transposed_form_and_mul_ilambda2(PublicKey->A, iP2);

    // precompute parity check for decryption, given paritycheck H, compute H*B^-1
    u32 t_entry;
    for (u16 row = 0; row < KDIM; row++)
    {
        for (u16 col = 0; col < N; col++)
        {
            t_entry = 0;
            for (u16 k = 0; k < N; k+=8)
                t_entry = (t_entry + H[row][k  ]*PrivateKey->inverseB[k  ][col]
                           + H[row][k+1]*PrivateKey->inverseB[k+1][col]
                           + H[row][k+2]*PrivateKey->inverseB[k+2][col]
                           + H[row][k+3]*PrivateKey->inverseB[k+3][col]
                           + H[row][k+4]*PrivateKey->inverseB[k+4][col]
                           + H[row][k+5]*PrivateKey->inverseB[k+5][col]
                           + H[row][k+6]*PrivateKey->inverseB[k+6][col]
                           + H[row][k+7]*PrivateKey->inverseB[k+7][col]);
            PrivateKey->BCode[row][col] = t_entry % Q;
        }
    }

}




