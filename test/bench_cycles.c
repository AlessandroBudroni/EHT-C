/*  This file is part of EHT-C.
 *
 *   name file: bench_cycles.c
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

#include <sys/stat.h>
#include <sys/time.h>

#define BENCH_LOOPS 10

/* run benchmark test - count CPU loops */
static void benchmark()
{

    time_stamp("Parameters: n %d, q %d, k %d, m %d, lambda^2 %d, sigma %f", (int)N, (int)Q, (int)K, (int)M, (int)LAMBDA2, (float)SIGMA);

    unsigned long long cycles_keygen = 0, cycles_encaps = 0, cycles_decaps = 0, cycles1, cycles2;

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

    // generate client's master secret
    FP ClientMasterSecret[N], ServerMasterSecret[N];

    for (int i = 0; i < BENCH_LOOPS; ++i)
    {
        cycles1 = cpucycles();
#ifdef FULL_STACK
        EHT_keygen(&PrivateKey, &PublicKey, H);
#else
        EHT_keygen(&PrivateKey, &PublicKey, &H);
#endif
        cycles2 = cpucycles();
        cycles_keygen = cycles_keygen+(cycles2-cycles1);

        // Benchmarking encapsulation
        cycles1 = cpucycles();
#ifdef FULL_STACK
        generate_secret(ClientMasterSecret, G);
#else
        generate_secret(ClientMasterSecret, &G);
#endif
        EHT_encrypt(ClientMasterSecret, &CipherText, &PublicKey);
        cycles2 = cpucycles();
        cycles_encaps = cycles_encaps+(cycles2-cycles1);

        // Benchmarking decapsulation
        cycles1 = cpucycles();
        EHT_decrypt(ServerMasterSecret, &CipherText, &PrivateKey);
        cycles2 = cpucycles();
        cycles_decaps = cycles_decaps+(cycles2-cycles1);
    }

    unsigned long long keygen_cycles = ((cycles_keygen/BENCH_LOOPS)/1000);
    unsigned long long encryption_cycles = ((cycles_encaps/BENCH_LOOPS)/1000);
    unsigned long long decryption_cycles = ((cycles_decaps/BENCH_LOOPS)/1000);
    unsigned long long tot_cycles = encryption_cycles + decryption_cycles;


    printf("  Key generation runs in ....................................... %10lld ", keygen_cycles  );
    printf("cycles");
    printf("\n");
    printf("  Encryption runs in ........................................... %10lld ", encryption_cycles  );
    printf("cycles");
    printf("\n");
    printf("  Decryption runs in ........................................... %10lld ", decryption_cycles  );
    printf("cycles");
    printf("\n");
    printf("  Encryption + Decryption runs in .............................. %10lld ", tot_cycles  );
    printf("cycles");
    printf("\n");

#ifndef FULL_STACK
    free_distribution(&PrivateKey);
    free_privateKey(&PrivateKey);
    free_publicKey(&PublicKey);
    free_cipherText(&CipherText);
    free_matrix(&G);
    free_matrix(&H);
#endif

}

int main()
{

    benchmark();
    return 0;
}
