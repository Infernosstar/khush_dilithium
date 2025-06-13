#include <stdio.h>
#include <stdint.h>
#include <string.h>
#include <stdlib.h>  // Added for malloc and free
#include "key_sign.h"
#include "params.h"
#include "key_gen.h"
#include "externalfunctions.h"
// Forward declaration of BitPack from key_sign.c
// void BitPack(uint8_t *output, const int32_t w[256], int32_t a, int32_t b);

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

int main() {
    // Sample test data initialization
    uint8_t rho[SEEDBYTES] = {0};
    uint8_t K[SEEDBYTES] = {0};
    uint8_t tr[CRHBYTES] = {0};
    int32_t s1[ML_L][256] = {{0}};
    int32_t s2[ML_K][256] = {{0}};
    int32_t t0[ML_K][256] = {{0}};

    // Initialize sample data with some values
    for (int i = 0; i < SEEDBYTES; i++) {
        rho[i] = (uint8_t)(i + 1);
        K[i] = (uint8_t)(i + 2);
    }
    for (int i = 0; i < CRHBYTES; i++) {
        tr[i] = (uint8_t)(i + 3);
    }
    for (int i = 0; i < ML_L; i++) {
        for (int j = 0; j < 256; j++) {
            s1[i][j] = (i + j) % (2 * ETA + 1) - ETA;
        }
    }
    for (int i = 0; i < ML_K; i++) {
        for (int j = 0; j < 256; j++) {
            s2[i][j] = (i + j) % (2 * ETA + 1) - ETA;
            t0[i][j] = (i + j) % ((1 << (D - 1)) - 1);
        }
    }

    // Buffer for encoded secret key
    size_t sk_len = SEEDBYTES * 2 + CRHBYTES +
                    ML_L * ((256 * (compute_bitlen(2 * ETA)) + 7) / 8) * 2 +
                    ML_K * ((256 * D + 7) / 8);
    uint8_t *sk = (uint8_t *)malloc(sk_len);
    if (!sk) {
        printf("Memory allocation failed\n");
        return 1;
    }

    // Encode
    skEncode(sk,rho, K, tr, s1,s2,t0);

    // Decode
    uint8_t rho_dec[SEEDBYTES];
    uint8_t K_dec[SEEDBYTES];
    uint8_t tr_dec[CRHBYTES];
    int32_t s1_dec[ML_L][256];
    int32_t s2_dec[ML_K][256];
    int32_t t0_dec[ML_K][256];

    skDecode(sk, rho_dec, K_dec, tr_dec, s1_dec, s2_dec, t0_dec);

    // Compare results
    int success = 1;
    if (memcmp(rho, rho_dec, SEEDBYTES) != 0) {
        printf("rho mismatch\n");
        success = 0;
    }
    if (memcmp(K, K_dec, SEEDBYTES) != 0) {
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
    return success ? 0 : 1;
}
