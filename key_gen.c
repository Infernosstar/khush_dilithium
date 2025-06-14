#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "randombytes.h"
#include "fips202.h"
#include "externalfunctions.h"
#include "params.h"
#include "key_gen.h"
#include "key_sign.h"
// Helper function to compute bit length (copied from key_sign.c)
int compute_bitlen(int value) {
    int bitlen = 0;
    while (value > 0) {
        bitlen++;
        value >>= 1;
    }
    return bitlen;
}

// // skEncode function: reverse of skDecode
// void skEncode(const uint8_t rho[SEEDBYTES],
//               const uint8_t K[SEEDBYTES],
//               const uint8_t tr[CRHBYTES],
//               const int32_t s1[ML_L][256],
//               const int32_t s2[ML_K][256],
//               const int32_t t0[ML_K][256],
//               uint8_t *sk) {

//     int offset = 0;
//     memcpy(sk + offset, rho, SEEDBYTES);
//     offset += SEEDBYTES;
//     memcpy(sk + offset, K, SEEDBYTES);
//     offset += SEEDBYTES;
//     memcpy(sk + offset, tr, CRHBYTES);
//     offset += CRHBYTES;

//     int bitlen_2eta = compute_bitlen(2 * ETA);
//     int bytes_per_poly_2eta = (256 * bitlen_2eta + 7) / 8;

//     for (int i = 0; i < ML_L; i++) {
//         BitPack(sk + offset, s1[i], 2 * ETA, ETA);
//         offset += bytes_per_poly_2eta;
//     }

//     for (int i = 0; i < ML_K; i++) {
//         BitPack(sk + offset, s2[i], 2 * ETA, ETA);
//         offset += bytes_per_poly_2eta;
//     }

//     int bitlen_t0 = D;
//     int bytes_per_poly_t0 = (256 * bitlen_t0 + 7) / 8;
//     for (int i = 0; i < ML_K; i++) {
//         BitPack(sk + offset, t0[i], (1 << (D - 1)) - 1, (1 << (D - 1)) - 1);
//         offset += bytes_per_poly_t0;
//     }
// }

// Function to compare two int32_t arrays
int compare_poly(const int32_t a[256], const int32_t b[256]) {
    for (int i = 0; i < 256; i++) {
        if (a[i] != b[i]) return 0;
    }
    return 1;
}

// Function to compare 2D arrays
int compare_poly_2d(int32_t a[][256], int32_t b[][256], int rows) {
    for (int i = 0; i < rows; i++) {
        if (!compare_poly(a[i], b[i])) return 0;
    }
    return 1;
}
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
        // Decode
    uint8_t rho_dec[SEEDBYTES];
    uint8_t K_dec[SEEDBYTES];
    uint8_t tr_dec[CRHBYTES];
    int32_t s1_dec[ML_L][256];
    int32_t s2_dec[ML_K][256];
    int32_t t0_dec[ML_K][256];
int32_t t1_dec[ML_K][256];
    skDecode(sk, rho_dec, K_dec, tr_dec, s1_dec, s2_dec, t0_dec);

    // Compare results
    int success = 1;
    if (memcmp(rho, rho_dec, SEEDBYTES) != 0) {
        printf("rho mismatch\n");
        success = 0;
    }
    if (memcmp(key, K_dec, SEEDBYTES) != 0) {
        printf("K mismatch\n");
        success = 0;
    }
    if (memcmp(tr, tr_dec, CRHBYTES) != 0) {
        printf("tr mismatch\n");
        success = 0;
    }
    if (!compare_poly_2d(s1, s1_dec, ML_L)) {
        printf("s1 mismatch\n");
        success = 0;
    }
    if (!compare_poly_2d(s2, s2_dec, ML_K)) {
        printf("s2 mismatch\n");
        success = 0;
    }
    if (!compare_poly_2d(t0, t0_dec, ML_K)) {
        printf("t0 mismatch\n");
        success = 0;
    }

    if (success) {
        printf("Test passed: skDecode(skEncode(x)) == x\n");
    } else {
        printf("Test failed: mismatch found\n");
    }

    free(sk);
    
    // Debug output to verify no leading zeros
   /* printf("Public key (first 32 bytes): ");
    for (int i = 0; i < 32; i++) printf("%02x", pk[i]);
    printf("\nSecret key (first 32 bytes): ");
    for (int i = 0; i < 32; i++) printf("%02x", sk[i]);
    printf("\n");*/
}