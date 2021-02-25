/*  This file is part of EHT-C.
 *
 *   name file: utils.h
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

#ifndef UTILS_H
#define UTILS_H

#include "config.h"

#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <math.h>

/* check and print */
void ASSERT(int condition, const char* string, ...);
void time_stamp (const char* string, ...);
int64_t cpucycles(void);


/* FP utils */
FP make_FP (int32_t input);
FP inverse_mod(int a);

/* Random utils */
int get_seed();

/* others */
double normal_PDF(int x, double mean, double sigma);
void FWHT(FP* data, int size);

#endif