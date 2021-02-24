/*  This file is part of EHT-C.
 *
 *   name file: decrypt.h
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


#ifndef DECRYPT_H
#define DECRYPT_H

#include <stdio.h>
#include <stdlib.h>

#include "config.h"
#include "utils.h"
#include "matrix.h"
#include "ciphersuite.h"
#include "encoding.h"

int EHT_decrypt(FP MasterSecret[N], cipherText *CipherText, privateKey *PrivateKey);

#endif