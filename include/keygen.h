/*  This file is part of EHT-C.
 *
 *   name file: keygen.h
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


#ifndef KEYGEN_H
#define KEYGEN_H

#include <stdio.h>
#include <stdlib.h>

#include "config.h"
#include "utils.h"
#include "matrix.h"
#include "ciphersuite.h"

// Precomputation
void precompute_distribuion(privateKey *PrivateKey);

// Hadamard key generation
#ifdef FULL_STACK
void EHT_keygen(privateKey *PrivateKey, publicKey *PublicKey, FP H[][N]);
#else
void EHT_keygen(privateKey *PrivateKey, publicKey *PublicKey, matrix *H);

void calloc_privateKey(privateKey *PrivateKey);
void calloc_publicKey(publicKey *PublicKey);
void calloc_cipherText(cipherText *CipherText);
void calloc_distribution(privateKey *PrivateKey);

void free_privateKey(privateKey *PrivateKey);
void free_distribution(privateKey *PrivateKey);
void free_publicKey(publicKey *PublicKey);
void free_cipherText(cipherText *CipherText);

#endif

#endif