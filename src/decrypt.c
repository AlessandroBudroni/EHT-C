/*  This file is part of EHT-C.
 *
 *   name file: decrypt.c
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
#include <string.h>
#include <time.h>

#include "decrypt.h"
#include "encoding.h"

#ifdef FULL_STACK
void vector_permutation(FP result[M], FP vector[M], u16 P[M])
#else
void vector_permutation(FP *result, FP *vector, u16 *P)
#endif /* FULL_STACK */
{
    for (int i = 0; i < M; ++i)
        result[i] = vector[P[i]];
}

/* inplace row-permutation */
#ifdef FULL_STACK
void inplace_vector_permutation(FP vector[M], u16 P[M])
#else
void inplace_vector_permutation(FP *vector, u16 *P)
#endif /* FULL_STACK */
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
            tempv = vector[pos];
            vector[pos] = vector[dest];
            vector[dest] = tempv;
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


static FP b[N];
static u16 tot_extra = 0;
static FP extra_residue[MAX_EXTRA];
static FP extra_position[MAX_EXTRA];
#ifdef FULL_STACK
static FP (*ptrBCode)[N];
#else
static matrix *ptrBCode;
#endif

int recursive_find_right_b(FP index)
{

    u16 pos = extra_position[index];
    u16 index_init = index;
    int ret;

    while(index < tot_extra && extra_position[index] == pos)
    {
        index++;
    }

    if (index == tot_extra)
    {
        for (int i = index_init; i < index; ++i)
        {
            b[extra_position[i]] = extra_residue[i];
            ret = decode(b, ptrBCode);
            if (ret)
                return 1;
        }
    }
    else
    {
        for (int i = index_init; i < index; ++i)
        {
            b[extra_position[i]] = extra_residue[i];
            ret = recursive_find_right_b(index);
            if (ret)
                return 1;
        }
    }

    return 0;
}


/* decrypt and decode the message - return 1 on success, return 0 on failure
	works only when T has 1 as the first element of each chunk
 */
int EHT_decrypt(FP MasterSecret[N], cipherText *CipherText, privateKey *PrivateKey)
{

    ptrBCode = PrivateKey->BCode;

    // memset b[N];
    tot_extra = 0;
    memset(extra_residue, 0, MAX_EXTRA*sizeof(FP));
    memset(extra_position, 0, MAX_EXTRA*sizeof(FP));


    FP z[M];
    // fast z = C*y
    vector_permutation(z, CipherText->y, PrivateKey->P2);
    for (int i = 0; i < M/LAMBDA2; ++i)
        FWHT(z+LAMBDA2*i, LAMBDA2);
    inplace_vector_permutation(z, PrivateKey->P1);

    double S, pi;
    int skip, ik, last_stop = 0;

    int accepted, ret;
    FP entry_pce, r_min, r_max;

    int c = ceil(Q/4);

    FP saved_entry[K];

    for (int i = 0; i < N; ++i)
    {
        accepted = 0;
        ik = i*K;

        r_min = make_FP(z[ik]-c);
        r_max = make_FP(z[ik]+c);

        for (int residue = r_min; residue !=r_max; residue=(residue+1)%Q)
        {
            S = 0;
            skip = 0;

            for (u16 j = 0; j < K; j++)
            {
                if(residue == r_min)
                {
                    entry_pce = (PrivateKey->T[ik +j]*residue+Q-z[ik +j])%Q;
                    saved_entry[j] = entry_pce;
                }
                else if(j<=last_stop)
                {
                    entry_pce = (saved_entry[j] +PrivateKey->T[ik +j])%Q;
                    saved_entry[j] = entry_pce;

                }
                else
                {
                    entry_pce = (PrivateKey->T[ik +j]*residue+Q-z[ik +j])%Q;
                    saved_entry[j] = entry_pce;
                }

                pi = PrivateKey->pce[entry_pce];

                if (pi < -9999998)  // be careful to this condition since it is a double
                {
                    skip = 1;
                    last_stop = j;
                    break;
                }
                S += pi;
            }

            if((!skip) && (S > 0))  // it is good, extend the solution
            {
                last_stop = K-1;
                if (accepted >= 1)
                {
                    if (tot_extra == MAX_EXTRA || tot_extra == MAX_EXTRA -1 && accepted == 1)
                    {
                        printf("too many false accepted %d\n", tot_extra);
                        exit(-1);
                    }
                    if (accepted == 1)
                    {
                        extra_residue[tot_extra] = b[i];
                        extra_position[tot_extra] = i;
                        tot_extra++;
                    }
                    extra_residue[tot_extra] = residue;
                    extra_position[tot_extra] = i;
                    tot_extra++;
                    accepted++;
                }
                else
                {
                    b[i] = residue;
                    accepted++;
                }
            }
        }
        if(accepted == 0) // right solution rejected
            return -1;
    }

    if (!tot_extra)
        ret = decode(b, PrivateKey->BCode);
    else
        ret = recursive_find_right_b(0);

    if (!ret)
        return 0;

    fast_matrix_mul_times_vector(MasterSecret, PrivateKey->inverseB, b);

    return 1;
}
