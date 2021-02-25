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

}


// Key generation using Hadamard matrices

/* instead of permutation matrix, we store a list that represent the permutation matrix. This is used to speed up key generation and avoid large matrix multiplication */
#ifdef FULL_STACK
void get_permutations(u16 P1[M], u16 P2[M])
#else
void get_permutations(u16 *P1, u16 *P2)
#endif
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
#ifdef FULL_STACK
void inplace_row_permutation_in_transposed_form(FP A[N][M], u16 P[M])
#else
void inplace_row_permutation_in_transposed_form(matrix *A, u16 *P)
#endif
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
#ifdef FULL_STACK
                tempv = A[i][pos];
                A[i][pos] = A[i][dest];
                A[i][dest] = tempv;
#else
                tempv = get_matrix_entry(A, i, pos);
                set_matrix_entry(A, i, pos, get_matrix_entry(A, i, dest));
                set_matrix_entry(A, i, dest, tempv);
#endif /* FULL_STACK */
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
#ifdef FULL_STACK
void inplace_row_permutation_in_transposed_form_and_mul_ilambda2(FP A[N][M], u16 P[M])
#else
void inplace_row_permutation_in_transposed_form_and_mul_ilambda2(matrix *A, u16 *P)
#endif /* FULL_STACK */
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
#ifdef FULL_STACK
                tempv = A[i][pos];
                A[i][pos] = A[i][dest];
                A[i][dest] = tempv;
                A[i][pos] = (A[i][pos]*invlambda2)%Q;
#else
                tempv = get_matrix_entry(A, i, pos);
                set_matrix_entry(A, i, pos, get_matrix_entry(A, i, dest));
                set_matrix_entry(A, i, dest, tempv);
                set_matrix_entry(A, i, pos, (get_matrix_entry(A, i, pos)*invlambda2)%Q );
#endif /* FULL_STACK */
            }
        }
        else
        {
            for (int i = 0; i < N; ++i)
#ifdef FULL_STACK
                A[i][pos] = (A[i][pos]*invlambda2)%Q;
#else
                set_matrix_entry(A, i, pos, (get_matrix_entry(A, i, pos)*invlambda2)%Q );
#endif /* FULL_STACK */
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


/* generate public/private keypair with the Hadamard matrix-based approach */
#ifdef FULL_STACK
void EHT_keygen(privateKey *PrivateKey, publicKey *PublicKey, FP H[][N])
#else
void EHT_keygen(privateKey *PrivateKey, publicKey *PublicKey, matrix *H)
#endif
{
    // generate private B and its inverse
#ifdef FULL_STACK
    FP B[N][N];
    do
    {
        random_matrix(B);
    }
    while(!invert_matrix(PrivateKey->inverseB, B) );
#else
    matrix B;
    calloc_matrix(&B, N, N);
    do
    {
        random_matrix(&B);
    }
    while(!invert_matrix(PrivateKey->inverseB, &B) );
#endif

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
    for (int row = 0; row < M; row++)
        for (int col = 0; col < N; col++)
#ifdef FULL_STACK
            PublicKey->A[col][row] = (PrivateKey->T[row]*B[row/K][col]) % Q;
#else
            set_matrix_entry(PublicKey->A, col, row, (PrivateKey->T[row]*get_matrix_entry(&B, row/K, col)) % Q);
#endif /* FULL_STACK */

    // permute first TB with the inverse of P2
    inplace_row_permutation_in_transposed_form(PublicKey->A, iP1);

    // apply FHWT to the transpose of TB already permuted by iP1
    for (int row = 0; row < N; row++)
        for (int col = 0; col < M/LAMBDA2; col++)
#ifdef FULL_STACK
            FWHT(PublicKey->A[row]+LAMBDA2*col, LAMBDA2);
#else
            FWHT(PublicKey->A->buff+row*M+LAMBDA2*col, LAMBDA2);
#endif /* FULL_STACK */

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
#ifdef FULL_STACK
                t_entry += H[row][k]*PrivateKey->inverseB[k][col]
                           + H[row][k+1]*PrivateKey->inverseB[k+1][col]
                           + H[row][k+2]*PrivateKey->inverseB[k+2][col]
                           + H[row][k+3]*PrivateKey->inverseB[k+3][col]
                           + H[row][k+4]*PrivateKey->inverseB[k+4][col]
                           + H[row][k+5]*PrivateKey->inverseB[k+5][col]
                           + H[row][k+6]*PrivateKey->inverseB[k+6][col]
                           + H[row][k+7]*PrivateKey->inverseB[k+7][col];
            PrivateKey->BCode[row][col] = t_entry % Q;

#else
                t_entry += get_matrix_entry(H, row, k)*get_matrix_entry(PrivateKey->inverseB, k, col)
                           + get_matrix_entry(H, row, k+1)*get_matrix_entry(PrivateKey->inverseB, k+1, col)
                           + get_matrix_entry(H, row, k+2)*get_matrix_entry(PrivateKey->inverseB, k+2, col)
                           + get_matrix_entry(H, row, k+3)*get_matrix_entry(PrivateKey->inverseB, k+3, col)
                           + get_matrix_entry(H, row, k+4)*get_matrix_entry(PrivateKey->inverseB, k+4, col)
                           + get_matrix_entry(H, row, k+5)*get_matrix_entry(PrivateKey->inverseB, k+5, col)
                           + get_matrix_entry(H, row, k+6)*get_matrix_entry(PrivateKey->inverseB, k+6, col)
                           + get_matrix_entry(H, row, k+7)*get_matrix_entry(PrivateKey->inverseB, k+7, col);
            set_matrix_entry(PrivateKey->BCode, row, col, t_entry % Q);
#endif /* FULL_STACK */
        }
    }
#ifndef FULL_STACK
    free_matrix(&B);
#endif
}


#ifndef FULL_STACK

void calloc_privateKey(privateKey *PrivateKey)
{

    PrivateKey->P1 = (FP*)calloc(M, sizeof(FP));
    PrivateKey->P2 = (FP*)calloc(M, sizeof(FP));
    PrivateKey->T = (FP*)calloc(M, sizeof(FP));
    PrivateKey->inverseB = (matrix*)calloc(1, sizeof(matrix));
    PrivateKey->BCode = (matrix*)calloc(1, sizeof(matrix));

    calloc_matrix(PrivateKey->inverseB, N, N);
    calloc_matrix(PrivateKey->BCode, KDIM, N);

    ASSERT(PrivateKey->P1 != NULL && PrivateKey->P2 != NULL &&
           PrivateKey->T != NULL && PrivateKey->inverseB != NULL, "failed allocating privateKey");
}

void calloc_publicKey(publicKey *PublicKey)
{

    PublicKey->A = (matrix*)calloc(1, sizeof(matrix));

    calloc_matrix(PublicKey->A, N, M); // transposed

    ASSERT(PublicKey->A != NULL, "failed allocating publicKey");
}

void calloc_cipherText(cipherText *CipherText)
{

    CipherText->y = (FP*)calloc(M, sizeof(FP));

    ASSERT(CipherText->y != NULL, "failed allocating cipherText");
}

void calloc_distribution(privateKey *PrivateKey)
{

    PrivateKey->pce = (double*)calloc(Q, sizeof(double));

    ASSERT(PrivateKey->pce != NULL, "failed allocating distribution for privateKey");
}

void free_privateKey(privateKey *PrivateKey)
{

    free(PrivateKey->P1);
    free(PrivateKey->P2);
    free(PrivateKey->T);
    free_matrix(PrivateKey->inverseB);
    free(PrivateKey->inverseB);
    free_matrix(PrivateKey->BCode);
    free(PrivateKey->BCode);
}

void free_distribution(privateKey *PrivateKey)
{

    free(PrivateKey->pce);
}

void free_publicKey(publicKey *PublicKey)
{

    free_matrix(PublicKey->A);
    free(PublicKey->A);
}

void free_cipherText(cipherText *CipherText)
{

    free(CipherText->y);
}

#endif /* FULL_STACK */

