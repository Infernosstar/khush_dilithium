#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "randombytes.h"
#include "fips202.h"
#include "externalfunctions.h"
#include "params.h"
#include "key_gen.h"

void keygen(uint8_t *pk, uint8_t *sk) {
   // printf("reached rthis function");
    // Key generation variables
    uint8_t xi[SEEDBYTES];
    uint8_t rho[SEEDBYTES], rhoprime[CRHBYTES], key[SEEDBYTES];
    uint8_t tr[CRHBYTES];
    
    // Polynomial arrays
    int32_t A[ML_K][ML_L][256];
    int32_t s1[ML_L][256], s2[ML_K][256];
    int32_t t[ML_K][256], t1[ML_K][256], t0[ML_K][256];
    
    // 1. Generate random seed
    randombytes(xi, SEEDBYTES);
    
    // 2. Expand seed (H(xi||k||l, 128))
    uint8_t input[SEEDBYTES+2];
    memcpy(input, xi, SEEDBYTES);
    input[SEEDBYTES] = ML_K;
    input[SEEDBYTES+1] = ML_L;
    
    uint8_t output[128];
    shake128(output, 128, input, SEEDBYTES+2);
    memcpy(rho, output, SEEDBYTES);
    memcpy(rhoprime, output+SEEDBYTES, CRHBYTES);
    memcpy(key, output+SEEDBYTES+CRHBYTES, SEEDBYTES);
    
    // 3. Generate matrix A
    expandA(rho, A);
    
    // 4. Generate secret vectors
    uint8_t s_seed[(SEEDBYTES<<1)+2];
    memcpy(s_seed, rhoprime, SEEDBYTES<<1);
    expandS(s_seed, (int32_t *)s1, (int32_t *)s2);
    
    // 5. Compute t = As1 + s2
    for (int i = 0; i < ML_K; i++) {
        for (int j = 0; j < 256; j++) {
            t[i][j] = 0;
            for (int k = 0; k < ML_L; k++) {
                // Multiply A[i][k] with s1[k] using NTT
                int32_t temp[256];
                NTT(s1[k]);
                for (int l = 0; l < 256; l++) {
                    temp[l] = (A[i][k][l] * s1[k][l]) % Q;
                }
                NTT_Inverse(temp);
                t[i][j] = (t[i][j] + temp[j]) % Q;
            }
            t[i][j] = (t[i][j] + s2[i][j]) % Q;
        }
    }
    
    // 6. Power2Round decomposition
    for (int i = 0; i < ML_K; i++) {
        for (int j = 0; j < 256; j++) {
            Power2Round(t[i][j], &t1[i][j], &t0[i][j]);
        }
    }
    
    // 7. Encode keys
    pkEncode(pk, rho, t1);
    shake128(tr, CRHBYTES, pk, PKBYTES);
    skEncode(sk, rho, key, tr, s1, s2, t0);
    
    printf("Public key generated (%d bytes)\n", PKBYTES);
    printf("Secret key generated (%d bytes)\n", SKBYTES);
}