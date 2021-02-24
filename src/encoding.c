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

void precompute_G_H(FP G[][N-KDIM], FP H[][N])
{

    // Generator

    memset(G, 0, N*(N-KDIM)*sizeof(FP));

    for (int i = 0; i < N-KDIM; ++i)
        G[i][i] = 1;

    int seed = 0;
    srand((unsigned) seed);

    for (int i = N-KDIM; i < N; ++i)
    {
        for (int j = 0; j < N-KDIM; ++j)
            G[i][j] = (FP)rand()%Q;
    }

    // Parity Check

    memset(H, 0, N*KDIM*sizeof(FP));
    for (int i = 0; i < KDIM; ++i)
        H[i][N-KDIM+i] = 1;

    for (int i = 0; i < KDIM; ++i)
    {
        for (int j = 0; j < N-KDIM; ++j)
        {
            H[i][j] = make_FP(-(int)G[N-KDIM+i][j]);
        }
    }
}


void encode(FP Secret[N], FP G[][N-KDIM])
{

    u32 sum = 0;
    for (u16 i = 0; i < KDIM; ++i)
    {
        sum = 0;
        Secret[N-KDIM+i] = 0;
        for (u16 j = 0; j < (N-8); j+=8)
        {
            sum = (sum + G[N-KDIM+i][j  ]*Secret[j  ]
                   + G[N-KDIM+i][j+1]*Secret[j+1]
                   + G[N-KDIM+i][j+2]*Secret[j+2]
                   + G[N-KDIM+i][j+3]*Secret[j+3]
                   + G[N-KDIM+i][j+4]*Secret[j+4]
                   + G[N-KDIM+i][j+5]*Secret[j+5]
                   + G[N-KDIM+i][j+6]*Secret[j+6]
                   + G[N-KDIM+i][j+7]*Secret[j+7]) %Q;
        }
        for (u16 j = (N-8); j < N-KDIM; ++j)
        {
            sum = (sum + G[N-KDIM+i][j]*Secret[j]) %Q;
        }
        Secret[N-KDIM+i] = sum;
    }
}


int decode(FP Codeword[N], FP H[][N])
{

    u32 sum = 0;
    for (u16 i = 0; i < KDIM; ++i)
    {
        sum = 0;
        for (u16 j = 0; j < N; j+=8)
            sum = (sum + (uint32_t)H[i][j]*Codeword[j]
                   + (uint32_t)H[i][j+1]*Codeword[j+1]
                   + (uint32_t)H[i][j+2]*Codeword[j+2]
                   + (uint32_t)H[i][j+3]*Codeword[j+3]
                   + (uint32_t)H[i][j+4]*Codeword[j+4]
                   + (uint32_t)H[i][j+5]*Codeword[j+5]
                   + (uint32_t)H[i][j+6]*Codeword[j+6]
                   + (uint32_t)H[i][j+7]*Codeword[j+7]) % Q;
        if (sum)
            return 0;
    }
    return 1;
}

