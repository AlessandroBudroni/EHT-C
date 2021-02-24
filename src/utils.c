/*  This file is part of EHT-C.
 *
 *   name file: utils.c
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

#include "utils.h"

#include <stdlib.h>
#include <math.h>
#include <stdarg.h>
#include <unistd.h>
#include <fcntl.h>

/* if condition is false, then print error message and exit */
void ASSERT(int condition, const char* string, ...)
{

    if (!condition)
    {
        va_list args;
        va_start(args, string);
        char err_msg[200];
        vsprintf(err_msg, string, args);
        printf("ERROR: %s\n", err_msg);
        va_end(args);
        exit(-1);
    }
}

/*
	time_stamp: make a time stamp
	INPUT
	- string: string to print
*/
void time_stamp (const char* string, ...)
{
    time_t my_time = time(NULL);
    va_list args;
    va_start(args, string);
    char msg[256];
    vsprintf(msg, string, args);
    printf("%.24s - %s\n", ctime(&my_time), msg);
    va_end(args);
}


FP make_FP (int32_t input)
{

    while (input < 0)
    {
        input += Q;
    }

    if (input >= Q)
        return (input % Q);

    return input;
}

/* get seed for RNG from /dev/urandom */
int get_seed()
{

    // get seed for RNG
    int randomData = open("/dev/urandom", O_RDONLY);
    ASSERT(randomData  > 0, "could not open /dev/urandom in get_seed");
    int seed, ret;
    ret = read(randomData, &seed, sizeof(int));
    ASSERT(ret > 0, "could not read from /dev/urandom in get_seed");
    close(randomData);

    return seed;
}

/* Access system counter for benchmarking */
int64_t cpucycles(void)
{
    unsigned int hi, lo;
    __asm__ volatile ("rdtsc\n\t" : "=a" (lo), "=d"(hi));
    return ((int64_t)lo) | (((int64_t)hi) << 32);
}

const double pi = 3.141592653589793;
const double e = 2.718281828459045;

/* Gaussian Distribution PDF */
double normal_PDF(int x, double mean, double sigma)
{
    return (1/(sigma*sqrt(2*pi)))*exp(-0.5*(((x-mean)/sigma)*((x-mean)/sigma)));
}

// Recursive - return x^y modulo Q
FP power_mod(FP x, FP y)
{
    if (y == 0)
        return 1;
    FP p = power_mod(x, y/2) % Q;
    p = (p * p) % Q;

    return (y%2 == 0)? p : (x * p) % Q;
}

/* Return the inverse of a modulo Q
   Assumption: m is prime, it uses Fermat's Little Theorem
*/
FP inverse_mod(int a)
{
    return power_mod(a, Q-2);
}

// Fast Walsh Hadamard Transform
void FWHT (FP* data, int size)
{
    u8 n = log2(size);
    FP tmp;
    uint16_t a, b;
    for (u8 i = 0; i < n; ++i)
    {
        for (u16 j = 0; j < (1 << n); j += 1 << (i+1))
        {
            for (u16 k = 0; k < (1<<i); ++k)
            {
                a = j + k;
                b = j + k + (1<<i);

                tmp = data[a];
                data[a] = (data[a] + data[b]) %Q;
                data[b] = (tmp +Q -data[b]) %Q;
            }
        }
    }
}