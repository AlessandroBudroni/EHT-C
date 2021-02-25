/*  This file is part of EHT-C.
 *
 *   name file: matrix.c
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

#include "matrix.h"
#include "string.h"


#ifndef FULL_STACK
/* (c)allocate memory for matrix*/
void calloc_matrix(matrix *Matrix, u16 N_rows, u16 N_cols)
{

    Matrix->buff = (FP*)calloc((u32)N_cols*(u32)N_rows, sizeof(FP));
    ASSERT(Matrix->buff != NULL, "failed allocating memory for matrix");

    Matrix->n_cols = N_cols;
    Matrix->n_rows = N_rows;
}

/* free memory for a matrix struct */
void free_matrix(matrix *Matrix)
{

    if(!Matrix->buff)
        return;
    if(Matrix->buff)
        free(Matrix->buff);

    Matrix->n_cols = 0;
    Matrix->n_rows = 0;
}

/* Set matrix entry, count from 0 to n_rows-1, from 0 to n_cols-1 */
void set_matrix_entry(matrix *Matrix, u16 row, u16 col, FP value)
{

    Matrix->buff[Matrix->n_cols*row + col] = value;
}

/* Get matrix entry, count from 0 to n_rows-1, from 0 to n_cols-1 */
FP get_matrix_entry(matrix *Matrix, u16 row, u16 col)
{

    return Matrix->buff[Matrix->n_cols*row + col];
}

#endif /* FULL_STACK */

/* Generate random matrix N_rows x n_cols, allocation must be done externally */
#ifdef FULL_STACK
void random_matrix(FP A[][N])
#else
void random_matrix(matrix *A)
#endif
{
    // get seed for RNG
    int seed = get_seed();

    // initialize RNG
    srand((unsigned) seed);

    for (u16 row = 0; row < N; row++)
        for (u16 col = 0; col < N; col++)
#ifdef FULL_STACK
            A[row][col] = (FP)(rand()%Q);
#else
            set_matrix_entry(A, row, col, (FP)(rand()%Q));
#endif
}

/* Return s = B*b */
#ifdef FULL_STACK
void fast_matrix_mul_times_vector(FP s[N], FP B[][N], FP b[N])
#else
void fast_matrix_mul_times_vector(FP *s, matrix *B, FP *b)
#endif
{
    FP tmp_entry[4];

    for (u16 row = 0; row < N; row+=4)
    {
        tmp_entry[0] = 0;
        tmp_entry[1] = 0;
        tmp_entry[2] = 0;
        tmp_entry[3] = 0;
        for (u16 k = 0; k < N; k+=8)
        {
#ifdef FULL_STACK
            tmp_entry[0] = (tmp_entry[0] + B[row  ][k]*b[k  ] + B[row  ][k+1]*b[k+1] + B[row  ][k+2]*b[k+2] + B[row  ][k+3]*b[k+3] + B[row  ][k+4]*b[k+4] + B[row  ][k+5]*b[k+5] + B[row  ][k+6]*b[k+6] + B[row  ][k+7]*b[k+7]) % Q;
            tmp_entry[1] = (tmp_entry[1] + B[row+1][k]*b[k  ] + B[row+1][k+1]*b[k+1] + B[row+1][k+2]*b[k+2] + B[row+1][k+3]*b[k+3] + B[row+1][k+4]*b[k+4] + B[row+1][k+5]*b[k+5] + B[row+1][k+6]*b[k+6] + B[row+1][k+7]*b[k+7]) % Q;
            tmp_entry[2] = (tmp_entry[2] + B[row+2][k]*b[k  ] + B[row+2][k+1]*b[k+1] + B[row+2][k+2]*b[k+2] + B[row+2][k+3]*b[k+3] + B[row+2][k+4]*b[k+4] + B[row+2][k+5]*b[k+5] + B[row+2][k+6]*b[k+6] + B[row+2][k+7]*b[k+7]) % Q;
            tmp_entry[3] = (tmp_entry[3] + B[row+3][k]*b[k  ] + B[row+3][k+1]*b[k+1] + B[row+3][k+2]*b[k+2] + B[row+3][k+3]*b[k+3] + B[row+3][k+4]*b[k+4] + B[row+3][k+5]*b[k+5] + B[row+3][k+6]*b[k+6] + B[row+3][k+7]*b[k+7]) % Q;
#else
            tmp_entry[0] = (tmp_entry[0] + get_matrix_entry(B, row, k)*b[k  ] + get_matrix_entry(B, row, k+1)*b[k+1] + get_matrix_entry(B, row, k+2)*b[k+2] + get_matrix_entry(B, row, k+3)*b[k+3] + get_matrix_entry(B, row, k+4)*b[k+4] + get_matrix_entry(B, row, k+5)*b[k+5] + get_matrix_entry(B, row, k+6)*b[k+6] + get_matrix_entry(B, row, k+7)*b[k+7]) % Q;
            tmp_entry[1] = (tmp_entry[1] + get_matrix_entry(B, row +1, k)*b[k  ] + get_matrix_entry(B, row +1, k+1)*b[k+1] + get_matrix_entry(B, row +1, k+2)*b[k+2] + get_matrix_entry(B, row +1, k+3)*b[k+3] + get_matrix_entry(B, row +1, k+4)*b[k+4] + get_matrix_entry(B, row +1, k+5)*b[k+5] + get_matrix_entry(B, row +1, k+6)*b[k+6] + get_matrix_entry(B, row +1, k+7)*b[k+7]) % Q;
            tmp_entry[2] = (tmp_entry[2] + get_matrix_entry(B, row +2, k)*b[k  ] + get_matrix_entry(B, row +2, k+1)*b[k+1] + get_matrix_entry(B, row +2, k+2)*b[k+2] + get_matrix_entry(B, row +2, k+3)*b[k+3] + get_matrix_entry(B, row +2, k+4)*b[k+4] + get_matrix_entry(B, row +2, k+5)*b[k+5] + get_matrix_entry(B, row +2, k+6)*b[k+6] + get_matrix_entry(B, row +2, k+7)*b[k+7]) % Q;
            tmp_entry[3] = (tmp_entry[3] + get_matrix_entry(B, row +3, k)*b[k  ] + get_matrix_entry(B, row +3, k+1)*b[k+1] + get_matrix_entry(B, row +3, k+2)*b[k+2] + get_matrix_entry(B, row +3, k+3)*b[k+3] + get_matrix_entry(B, row +3, k+4)*b[k+4] + get_matrix_entry(B, row +3, k+5)*b[k+5] + get_matrix_entry(B, row +3, k+6)*b[k+6] + get_matrix_entry(B, row +3, k+7)*b[k+7]) % Q;
#endif
        }
        s[row    ] = tmp_entry[0];
        s[row + 1] = tmp_entry[1];
        s[row + 2] = tmp_entry[2];
        s[row + 3] = tmp_entry[3];
    }

}

/* return 1 if A==B, otherwise 0 */
#ifdef FULL_STACK
int vector_equal(u16 n, FP A[n], FP B[n])
#else
int vector_equal(u16 n, FP *A, FP *B)
#endif
{
    for (int i = 0; i < n; i++)
    {
        if (A[i] != B[i])
        {
            // printf("Matrix not equal (%d %d): %d %d\n", i, j, A[i][j], B[i][j]);
            return 0;
        }
    }
    return 1;
}


/* return 1 if A==B, otherwise 0
int matrix_equal(u16 n_rows, u16 n_cols, FP A[][n_cols], FP B[][n_cols]){

    for (int i = 0; i < n_cols; i++)
    {
        for (int j = 0; j < n_rows; j++)
        {
            if (A[i][j] != B[i][j])
            {
                // printf("Matrix not equal (%d %d): %d %d\n", i, j, A[i][j], B[i][j]);
                return 0;
            }
        }
    }
    return 1;
}
*/

/* swap i-th and j-th rows of A */
#ifdef FULL_STACK
void swap_rows(u16 i, u16 j, u16 n_cols, FP A[][n_cols])
{

    FP tmpRow[n_cols];
    memcpy(tmpRow, A[i], n_cols*sizeof(FP));
    memcpy(A[i], A[j], n_cols*sizeof(FP));
    memcpy(A[j], tmpRow, n_cols*sizeof(FP));
}

/* Use Gaussian Elimination on the concatenation matrix C|I to get I|C^-1. Save C^-1 in inverseC
    return 1 if success, return 0 if matrix is not invertible.
*/
int invert_matrix(FP inverseC[][N], FP C[][N])
{

    FP CI[N][2*N];

    // copy matrix C in CI, then concatenate with I
    for (u16 i = 0; i < N; ++i)
    {
        memcpy(CI[i], C[i], N*sizeof(FP));
        memset(CI[i]+N, 0, N*sizeof(FP));
        CI[i][N+i] = 1;
    }

    //loop for the generation of diagonal matrix I|C^-1
    int pt;
    FP c;
    for(u16 j = 0; j < N; j++)
    {
        for(u16 i = 0; i < N; i++)
        {
            if (i > j) // upper triangolar
            {
                c = make_FP(CI[i][j] * inverse_mod(CI[j][j]));
                for(u16 k = 0; k < 2*N; k++)
                    CI[i][k] = make_FP(CI[i][k] -(int32_t)(c * CI[j][k]) %Q);
            }
            else if( i == j) // put 1 in diagonal
            {
                if (CI[i][j] == 0)
                {
                    pt = i;
                    while(1)
                    {
                        if (pt == N)
                        {
                            return 0; // matrix non invertible
                        }

                        if(CI[pt][j] == 0)
                            pt++;
                        else
                            break;
                    }
                    swap_rows(i, pt, 2*N, CI);
                }
                c = inverse_mod(CI[i][j]);
                for(u16 k = 0; k < 2*N; k++)
                {
                    CI[i][k] = make_FP((c * CI[i][k]) %Q);
                }
            }
        }
        for (int l = j-1; l >= 0; l--) // make diagonal
        {
            c = CI[l][j];
            if (c != 0)
            {
                for (int k = 0; k < 2*N; k++)
                {
                    CI[l][k] = make_FP(CI[l][k] -(int32_t)(c * CI[j][k]) %Q );
                }
            }
        }
    }



    // check that the matrix was actually invertible, or if for some reason it failed in putting the diagonal I
    for (int i = 0; i < N; i++)
    {
        if (CI[i][i] != 1)
        {
            return 0;
        }
    }

    // copy inverse of C
    for (int i = 0; i < N; i++)
        for (int j = 0; j < N; j++)
            inverseC[i][j] = CI[i][N+j];

    return 1;
}
#else
void swap_rows(matrix *A, u16 i, u16 j)
{

    u16 n_cols = A->n_cols;
    FP tmpRow[n_cols];
    memcpy(tmpRow, A->buff + n_cols*i, n_cols*sizeof(FP));
    memcpy(A->buff + n_cols*i, A->buff + n_cols*j, n_cols*sizeof(FP));
    memcpy(A->buff + n_cols*j, tmpRow, n_cols*sizeof(FP));
}

/* Use Gaussian Elimination on the concatenation matrix C|I to get I|C^-1. Save C^-1 in inverseC
    return 1 if success, return 0 if matrix is not invertible.
*/
int invert_matrix(matrix *inverseC, matrix *C)
{

    ASSERT(C->n_cols == C->n_rows && C->n_cols != 0, "trying to invert a non-square or an empty matrix");

    u16 n = C->n_cols; // order of the matrix

    matrix CI;
    calloc_matrix(&CI, n, 2*n);

    // copy matrix C in CI, then concatenate with I

    for (u16 i = 0; i < C->n_rows; i++)
        for (u16 j = 0; j < C->n_cols; j++)
            set_matrix_entry(&CI, i, j, get_matrix_entry(C, i, j));


    for (u16 i = 0; i < n; ++i)
        set_matrix_entry(&CI, i, n+i, 1);

    //loop for the generation of diagonal matrix I|C^-1
    int pt;
    FP c;
    for(u16 j = 0; j < n; j++)
    {
        for(u16 i = 0; i < n; i++)
        {
            if (i > j) // upper triangolar
            {
                c = make_FP(get_matrix_entry(&CI, i, j) * inverse_mod(get_matrix_entry(&CI, j, j)));
                for(u16 k = 0; k < 2*n; k++)
                    set_matrix_entry(&CI, i, k, make_FP(get_matrix_entry(&CI, i, k) -(c * get_matrix_entry(&CI, j, k)) %Q) );
            }
            else if( i == j) // put 1 in diagonal
            {
                if (get_matrix_entry(&CI, i, j) == 0)
                {
                    pt = i;
                    while(1)
                    {
                        if (pt == n)
                        {
                            free_matrix(&CI);
                            return 0; // matrix non invertible
                        }

                        if(get_matrix_entry(&CI, pt, j) == 0)
                            pt++;
                        else
                            break;
                    }
                    swap_rows(&CI, i, pt);
                }
                c = inverse_mod(get_matrix_entry(&CI, i, j));
                for(u16 k = 0; k < 2*n; k++)
                {
                    set_matrix_entry(&CI, i, k, make_FP((c * get_matrix_entry(&CI, i, k)) %Q) );
                }
            }
        }
        for (int l = j-1; l >= 0; l--) // make diagonal
        {
            c = get_matrix_entry(&CI, l, j);
            if (c != 0)
            {
                for (int k = 0; k < 2*n; k++)
                {
                    set_matrix_entry(&CI, l, k, make_FP(get_matrix_entry(&CI, l, k) -(c * get_matrix_entry(&CI, j, k)) %Q) );
                }
            }
        }
    }

    // check that the matrix was actually invertible, or if for some reason it failed in putting the diagonal I
    for (int i = 0; i < n; i++)
    {
        if (get_matrix_entry(&CI, i, i) != 1)
        {
            free_matrix(&CI);
            return 0;
        }
    }

    // copy inverse of C
    for (int i = 0; i < n; i++)
        for (int j = 0; j < n; j++)
            set_matrix_entry(inverseC, i, j, get_matrix_entry(&CI, i, n+j) );

    free_matrix(&CI);
    return 1;
}
#endif


// /* print matrix to stdout */
// void print_matrix(u16 n_rows, u16 n_cols, FP Matrix[][n_cols]){

// 	if(Matrix == NULL){
// 		printf("Empty matrix\n");
// 		return;
// 	}

// 	for (u16 i = 0; i < n_rows; i++)
// 	{
// 		for (u16 j = 0; j < n_cols; j++)
// 			printf("%lu ", (unsigned long)Matrix[i][j]);
// 		printf("\n");
// 	}
// 	printf("\n");
// }

// /* optimeized to print column vectors */
// void print_vector(u16 lenght, FP Vector[lenght]){

//     for (int i = 0; i < lenght; i++)
//         printf("%lu ", (unsigned long)Vector[i]);
//     printf("\n");
// }

