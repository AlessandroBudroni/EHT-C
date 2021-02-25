/*  This file is part of EHT-C.
 *
 *   name file: test_correctness.c
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
#include "matrix.h"
#include "encrypt.h"
#include "decrypt.h"
#include "ciphersuite.h"
#include "keygen.h"
#include "encoding.h"

#include <string.h>
#include <sys/stat.h>

/* test correctness */
int main(int argc, char const *argv[])
{
    time_stamp("Parameters: n %d, q %d, k %d, m %d, lambda^2 %d, sigma %f", (int)N, (int)Q, (int)K, (int)M, (int)LAMBDA2, (float)SIGMA);

    privateKey PrivateKey;
    publicKey PublicKey;
    cipherText CipherText;

#ifndef FULL_STACK
    calloc_privateKey(&PrivateKey);
    calloc_distribution(&PrivateKey);
    calloc_publicKey(&PublicKey);
    calloc_cipherText(&CipherText);
#endif

    time_stamp("Precompute distribution");
    precompute_distribuion(&PrivateKey);

    time_stamp("Precompute CDF table");
    precompute_cdf_table();

    time_stamp("Precompute generator G and parity-check H");
#ifdef FULL_STACK
    FP G[N][N-KDIM], H[KDIM][N];
    precompute_G_H(G, H);
#else
    matrix G, H;
    calloc_matrix(&G, N, N-KDIM);
    calloc_matrix(&H, KDIM, N);
    precompute_G_H(&G, &H);
#endif
    time_stamp("Done\n");

    // generate client's master secret
    FP ClientMasterSecret[N], ServerMasterSecret[N];

    int ret, failed = 0, detected = 0, success = 0, keygens = 10, iterations = 10;
    time_stamp("Start %d encryption/decryptions", iterations*keygens);

    for (int j = 0; j < keygens; j++)
    {
#ifdef FULL_STACK
        EHT_keygen(&PrivateKey, &PublicKey, H);
#else
        EHT_keygen(&PrivateKey, &PublicKey, &H);
#endif
        for (int i = 0; i < iterations; ++i)
        {
#ifdef FULL_STACK
            generate_secret(ClientMasterSecret, G);
#else
            generate_secret(ClientMasterSecret, &G);
#endif
            EHT_encrypt(ClientMasterSecret, &CipherText, &PublicKey);

            ret = EHT_decrypt(ServerMasterSecret, &CipherText, &PrivateKey);
            if (!ret)
                detected++;
            else if(!vector_equal(N, ServerMasterSecret, ClientMasterSecret))
                failed++;
            else
                success++;
        }
    }
    time_stamp("Failed %d, detected %d, success %d out of %d", failed, detected, success, iterations*keygens);

    if( keygens*iterations == success)
        time_stamp("Test passed");
    else
        time_stamp("Test failed");

#ifndef FULL_STACK
    free_distribution(&PrivateKey);
    free_privateKey(&PrivateKey);
    free_publicKey(&PublicKey);
    free_cipherText(&CipherText);
    free_matrix(&G);
    free_matrix(&H);
#endif

    return 0;
}