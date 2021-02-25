/*  This file is part of EHT-C.
 *
 *   name file: encrypt.h
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

#ifndef ENCRYPT_H
#define ENCRYPT_H

#include <stdint.h>
#include <stdlib.h>
#include <math.h>

#include "config.h"
#include "utils.h"
#include "matrix.h"
#include "ciphersuite.h"
#include "encoding.h"

void precompute_cdf_table();

#ifdef FULL_STACK
void generate_secret(FP Secret[N], FP G[][N-KDIM]);
#else
void generate_secret(FP Secret[N], matrix *G);
#endif
void EHT_encrypt(FP MasterSecret[N], cipherText *CipherText, publicKey *PublicKey);

#endif