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


/* Generate random matrix N_rows x n_cols, allocation must be done externally */
void random_matrix(FP A[][N])
{

    // get seed for RNG
    int seed = get_seed();

    // initialize RNG
    srand((unsigned) seed);

    for (u16 row = 0; row < N; row++)
        for (u16 col = 0; col < N; col++)
            A[row][col] = (FP)(rand()%Q);
}

/* Return s = B*b */
void fast_matrix_mul_times_vector(FP s[N], FP B[][N], FP b[N])
{

    // ASSERT( N % 8 == 0, "m is not a multiple of 8.");

    FP tmp_entry[4];

    for (u16 row = 0; row < N; row+=4)
    {
        tmp_entry[0] = 0;
        tmp_entry[1] = 0;
        tmp_entry[2] = 0;
        tmp_entry[3] = 0;
        for (u16 k = 0; k < N; k+=8)
        {
            tmp_entry[0] = (tmp_entry[0] + B[row  ][k]*b[k  ] + B[row  ][k+1]*b[k+1] + B[row  ][k+2]*b[k+2] + B[row  ][k+3]*b[k+3] + B[row  ][k+4]*b[k+4] + B[row  ][k+5]*b[k+5] + B[row  ][k+6]*b[k+6] + B[row  ][k+7]*b[k+7]) % Q;
            tmp_entry[1] = (tmp_entry[1] + B[row+1][k]*b[k  ] + B[row+1][k+1]*b[k+1] + B[row+1][k+2]*b[k+2] + B[row+1][k+3]*b[k+3] + B[row+1][k+4]*b[k+4] + B[row+1][k+5]*b[k+5] + B[row+1][k+6]*b[k+6] + B[row+1][k+7]*b[k+7]) % Q;
            tmp_entry[2] = (tmp_entry[2] + B[row+2][k]*b[k  ] + B[row+2][k+1]*b[k+1] + B[row+2][k+2]*b[k+2] + B[row+2][k+3]*b[k+3] + B[row+2][k+4]*b[k+4] + B[row+2][k+5]*b[k+5] + B[row+2][k+6]*b[k+6] + B[row+2][k+7]*b[k+7]) % Q;
            tmp_entry[3] = (tmp_entry[3] + B[row+3][k]*b[k  ] + B[row+3][k+1]*b[k+1] + B[row+3][k+2]*b[k+2] + B[row+3][k+3]*b[k+3] + B[row+3][k+4]*b[k+4] + B[row+3][k+5]*b[k+5] + B[row+3][k+6]*b[k+6] + B[row+3][k+7]*b[k+7]) % Q;
        }
        s[row    ] = tmp_entry[0];
        s[row + 1] = tmp_entry[1];
        s[row + 2] = tmp_entry[2];
        s[row + 3] = tmp_entry[3];
    }

}

/* return 1 if A==B, otherwise 0 */
int vector_equal(u16 n, FP A[n], FP B[n])
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

