/*  This file is part of EHT-C.
 *
 *   name file: encoding.c
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


#include <string.h>
#include <math.h>
#include "encoding.h"

/* Precompute generator and parity-check matrices */
#ifdef FULL_STACK
void precompute_G_H(FP G[][N-KDIM], FP H[][N])
#else
void precompute_G_H(matrix *G, matrix *H)
#endif
{

    // Generator
#ifdef FULL_STACK
    memset(G, 0, N*(N-KDIM)*sizeof(FP));
#endif

    for (int i = 0; i < N-KDIM; ++i)
#ifdef FULL_STACK
        G[i][i] = 1;
#else
        set_matrix_entry(G, i, i, 1);
#endif


    int seed = 0;
    srand((unsigned) seed);

    for (int i = N-KDIM; i < N; ++i)
    {
        for (int j = 0; j < N-KDIM; ++j)
#ifdef FULL_STACK
            G[i][j] = (FP)rand()%Q;
#else
            set_matrix_entry(G, i, j, (FP)rand()%Q);
#endif
    }

    // Parity Check
#ifdef FULL_STACK
    memset(H, 0, N*KDIM*sizeof(FP));
#endif

    for (int i = 0; i < KDIM; ++i)
#ifdef FULL_STACK
        H[i][N-KDIM+i] = 1;
#else
        set_matrix_entry(H, i, N-KDIM+i, 1);
#endif

    for (int i = 0; i < KDIM; ++i)
    {
        for (int j = 0; j < N-KDIM; ++j)
        {
#ifdef FULL_STACK
            H[i][j] = make_FP(-(int32_t)G[N-KDIM+i][j]);
#else
            set_matrix_entry(H, i, j, make_FP(-(int32_t) get_matrix_entry(G, N-KDIM+i, j)));
#endif
        }
    }
}

/* Encode message - compute only redoundancy part */
#ifdef FULL_STACK
void encode(FP Secret[N], FP G[][N-KDIM])
#else
void encode(FP *Secret, matrix *G)
#endif
{
    u32 sum = 0;
    for (u16 i = 0; i < KDIM; ++i)
    {
        sum = 0;
        Secret[N-KDIM+i] = 0;
        for (u16 j = 0; j < (N-8); j+=8)
        {
#ifdef FULL_STACK
            sum = (sum + G[N-KDIM+i][j  ]*Secret[j  ]
                   + G[N-KDIM+i][j+1]*Secret[j+1]
                   + G[N-KDIM+i][j+2]*Secret[j+2]
                   + G[N-KDIM+i][j+3]*Secret[j+3]
                   + G[N-KDIM+i][j+4]*Secret[j+4]
                   + G[N-KDIM+i][j+5]*Secret[j+5]
                   + G[N-KDIM+i][j+6]*Secret[j+6]
                   + G[N-KDIM+i][j+7]*Secret[j+7]) %Q;
#else
            sum = (sum + get_matrix_entry(G, N-KDIM+i, j)*Secret[j  ]
                   + get_matrix_entry(G, N-KDIM+i, j+1)*Secret[j+1]
                   + get_matrix_entry(G, N-KDIM+i, j+2)*Secret[j+2]
                   + get_matrix_entry(G, N-KDIM+i, j+3)*Secret[j+3]
                   + get_matrix_entry(G, N-KDIM+i, j+4)*Secret[j+4]
                   + get_matrix_entry(G, N-KDIM+i, j+5)*Secret[j+5]
                   + get_matrix_entry(G, N-KDIM+i, j+6)*Secret[j+6]
                   + get_matrix_entry(G, N-KDIM+i, j+7)*Secret[j+7]) %Q;
#endif
        }
        for (u16 j = (N-8); j < N-KDIM; ++j)
#ifdef FULL_STACK
            sum = (sum + G[N-KDIM+i][j]*Secret[j]) %Q;
#else
            sum = (sum + get_matrix_entry(G, N-KDIM+i, j)*Secret[j]) %Q;
#endif
        Secret[N-KDIM+i] = sum;
    }
}

/* decode message end detect if there have been errors */
#ifdef FULL_STACK
int decode(FP Codeword[N], FP H[][N])
#else
int decode(FP *Codeword, matrix *H)
#endif
{
    u32 sum = 0;
    for (u16 i = 0; i < KDIM; ++i)
    {
        sum = 0;
        for (u16 j = 0; j < N; j+=8)
#ifdef FULL_STACK
            sum = (sum + (uint32_t)H[i][j]*Codeword[j]
                   + H[i][j+1]*Codeword[j+1]
                   + H[i][j+2]*Codeword[j+2]
                   + H[i][j+3]*Codeword[j+3]
                   + H[i][j+4]*Codeword[j+4]
                   + H[i][j+5]*Codeword[j+5]
                   + H[i][j+6]*Codeword[j+6]
                   + H[i][j+7]*Codeword[j+7]) % Q;
#else
            sum = (sum + (uint32_t)get_matrix_entry(H, i, j)*Codeword[j]
                   + get_matrix_entry(H, i, j+1)*Codeword[j+1]
                   + get_matrix_entry(H, i, j+2)*Codeword[j+2]
                   + get_matrix_entry(H, i, j+3)*Codeword[j+3]
                   + get_matrix_entry(H, i, j+4)*Codeword[j+4]
                   + get_matrix_entry(H, i, j+5)*Codeword[j+5]
                   + get_matrix_entry(H, i, j+6)*Codeword[j+6]
                   + get_matrix_entry(H, i, j+7)*Codeword[j+7]) % Q;
#endif
        if (sum)
            return 0;
    }
    return 1;
}

