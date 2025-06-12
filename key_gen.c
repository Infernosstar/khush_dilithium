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
    
    // 2. Expand seed using SHAKE256 with proper domain separation
    keccak_state state;
    shake256_init(&state);
    
    // Absorb all seed material with proper padding
    uint8_t domain_sep[2] = {ML_K, ML_L};
    shake256_absorb(&state, xi, SEEDBYTES);
    shake256_absorb(&state, domain_sep, 2);
    shake256_finalize(&state);
    
    // Squeeze all required material at once
    uint8_t squeezed[(SEEDBYTES + CRHBYTES + SEEDBYTES)];
    shake256_squeeze(squeezed, sizeof(squeezed), &state);
    
    // Split the squeezed output without reinitializing
    memcpy(rho, squeezed, SEEDBYTES);
    memcpy(rhoprime, squeezed + SEEDBYTES, CRHBYTES);
    memcpy(key, squeezed + SEEDBYTES + CRHBYTES, SEEDBYTES);
    
    // 3. Generate matrix A
    expandA(rho, A);
    
    // 4. Generate secret vectors s1, s2
    keccak_state s_state;
    shake256_init(&s_state);
    shake256_absorb(&s_state, rhoprime, CRHBYTES);
    shake256_finalize(&s_state);
    
    uint8_t s_seed_expanded[(SEEDBYTES<<1)+2];
    shake256_squeeze(s_seed_expanded, sizeof(s_seed_expanded), &s_state);
    expandS(s_seed_expanded, (int32_t *)s1, (int32_t *)s2);
    
    // 5. Compute t = A*s1 + s2
    for (int i = 0; i < ML_K; i++) {
        for (int j = 0; j < 256; j++) {
            t[i][j] = 0;
            for (int k = 0; k < ML_L; k++) {
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
    
    // Compute tr = SHAKE256(pk) with proper finalization
    shake256_init(&state);
    shake256_absorb(&state, pk, PKBYTES);
    shake256_finalize(&state);
    shake256_squeeze(tr, CRHBYTES, &state);
    
    skEncode(sk, rho, key, tr, s1, s2, t0);
    
    // Debug output to verify no leading zeros
   /* printf("Public key (first 32 bytes): ");
    for (int i = 0; i < 32; i++) printf("%02x", pk[i]);
    printf("\nSecret key (first 32 bytes): ");
    for (int i = 0; i < 32; i++) printf("%02x", sk[i]);
    printf("\n");*/
}