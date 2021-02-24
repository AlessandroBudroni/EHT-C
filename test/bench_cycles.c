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

static void benchmark()
{

    privateKey PrivateKey;
    publicKey PublicKey;
    cipherText CipherText;
    FP ClientMasterSecret[N], ServerMasterSecret[N];

    // precomputation
    FP G[KDIM][N-KDIM], H[KDIM][N];
    precompute_cdf_table();
    precompute_distribuion(&PrivateKey);
    precompute_G_H(G, H);

    unsigned long long cycles_keygen = 0, cycles_encaps = 0, cycles_decaps = 0, cycles1, cycles2;

    time_stamp("Parameters: n %d, q %d, k %d, m %d, lambda^2 %d, sigma %f", (int)N, (int)Q, (int)K, (int)M, (int)LAMBDA2, (float)SIGMA);

    for (int i = 0; i < BENCH_LOOPS; ++i)
    {
        cycles1 = cpucycles();
        EHT_keygen(&PrivateKey, &PublicKey, H);
        cycles2 = cpucycles();
        cycles_keygen = cycles_keygen+(cycles2-cycles1);

        // Benchmarking encapsulation
        cycles1 = cpucycles();
        generate_secret(ClientMasterSecret, G);
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

}

int main()
{

    benchmark();
    return 0;
}
