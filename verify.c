#include <stdint.h>
#include <string.h>
#include "params.h"
#include "fips202.h"
#include "externalfunctions.h"
#include "Conversions.h"
#include <stdio.h>

// Function to verify a signature
int ml_dsa_verify_internal(const uint8_t *pk, const uint8_t *M_prime, size_t M_prime_len, const uint8_t *sigma) {
   int ml_dsa_verify_internal(const uint8_t *pk, const uint8_t *M_prime, size_t M_prime_len, const uint8_t *sigma) {
    // Step 1: Decode public key
    uint8_t rho[SEEDBYTES];
    int32_t t1[ML_K][256];
    pkDecode(pk, rho, t1);

    // Step 2: Decode signature
    uint8_t c_tilde[32];
    int32_t z[ML_L][256];
    int32_t h[ML_K][256];
    if (sigDecode(sigma, c_tilde, z, h) != 0) {
        return 0;
    }

    // Step 5: Expand matrix A in NTT domain
    int32_t A_hat[ML_K][ML_L][256];
    expandA(rho, A_hat);

    // Step 6: Compute tr = H(pk, 64)
    uint8_t tr[64];
    {
        keccak_state state;
        shake256_init(&state);
        shake256_absorb(&state, pk, PKBYTES);
        shake256_squeeze(tr, 64, &state);
    }

    // Step 7: Compute μ = H(BytesToBits(tr) || M', 64)
    uint8_t mu[64];
    {
        uint8_t tr_bits[CRHBYTES * 8];
        BytesToBits(tr, CRHBYTES, tr_bits);
        uint8_t input[CRHBYTES * 8 + M_prime_len];
        memcpy(input, tr_bits, CRHBYTES * 8);
        memcpy(input + CRHBYTES * 8, M_prime, M_prime_len);
        
        keccak_state state;
        shake256_init(&state);
        shake256_absorb(&state, input, CRHBYTES * 8 + M_prime_len);
        shake256_squeeze(mu, 64, &state);
    }

    // Step 8: Sample c from c_tilde
    int32_t c[256];
    {
        keccak_state state;
        shake256_init(&state);
        shake256_absorb(&state, c_tilde, 32);
        shake256_squeeze((uint8_t *)c, 256 * sizeof(int32_t), &state);
        // Additional processing may be needed here depending on SampleInBall implementation
    }


    // Step 9: Compute w_approx' = NTT^-1(A_hat ◦ NTT(z) - NTT(c) ◦ NTT(t1 * 2^d))
    int32_t w_approx_prime[ML_K][256];
    int32_t z_hat[ML_L][256];
    
    // Convert z to NTT domain
    for (int i = 0; i < ML_L; i++) {
        NTT(z_hat[i], z[i]);
    }

    // Compute A_hat ◦ NTT(z)
    for (int i = 0; i < ML_K; i++) {
        memset(w_approx_prime[i], 0, sizeof(w_approx_prime[i]));
        for (int j = 0; j < ML_L; j++) {
            int32_t temp[256];
            pointwise_mul(temp, A_hat[i][j], z_hat[j]);
            NTT_Inverse(temp);
            for (int k = 0; k < 256; k++) {
                w_approx_prime[i][k] = (w_approx_prime[i][k] + temp[k]) % Q;
                if (w_approx_prime[i][k] < 0) w_approx_prime[i][k] += Q;
            }
        }
    }

    // Subtract c ◦ t1 * 2^d
    for (int i = 0; i < ML_K; i++) {
        int32_t temp[256];
        pointwise_mul(temp, c_hat, t1_hat[i]);
        NTT_Inverse(temp);
        for (int j = 0; j < 256; j++) {
            w_approx_prime[i][j] = (w_approx_prime[i][j] - temp[j] * (1 << D)) % Q;
            if (w_approx_prime[i][j] < 0) w_approx_prime[i][j] += Q;
        }
    }

    // Step 10: Reconstruct w1' using hints
    int32_t w1_prime[ML_K][256];
    UseHint(h, w_approx_prime, w1_prime);

    // Step 12: Compute c_tilde' = H(μ || w1Encode(w1'), λ/4)
    uint8_t w1_encoded[ML_K * ((256 * (BITLEN_GAMMA1 * 2) + 7) / 8];
    w1Encode(w1_encoded, w1_prime);
    
    uint8_t c_tilde_prime[32];
    {
        uint8_t input[64 + sizeof(w1_encoded)];
        memcpy(input, mu, 64);
        memcpy(input + 64, w1_encoded, sizeof(w1_encoded));
        shake256(c_tilde_prime, 32, input, 64 + sizeof(w1_encoded));
    }

    // Step 13: Check conditions
    // 1. Check that c_tilde matches c_tilde_prime
    if (memcmp(c_tilde, c_tilde_prime, 32) != 0) {
        return 0;
    }

    // 2. Check that ||z||_∞ < γ1 - β
    for (int i = 0; i < ML_L; i++) {
        for (int j = 0; j < 256; j++) {
            if (z[i][j] >= GAMMA1 - BETA || z[i][j] <= -GAMMA1 + BETA) {
                return 0;
            }
        }
    }

    return 1; // Signature is valid
}

// Wrapper function that formats the message before verification
int ml_dsa_verify(const uint8_t *pk, const uint8_t *M, size_t M_len,
                 const uint8_t *tx, size_t tx_len, const uint8_t *sigma) {
    // Check context length
    if (tx_len > 255) {
        return 0;
    }

    // Format M′ = (0 || tx_len || tx || M)
    size_t M_prime_len = 2 + tx_len + M_len;
    uint8_t *M_prime = (uint8_t *)malloc(M_prime_len);
    if (!M_prime) return 0;

    M_prime[0] = 0;
    M_prime[1] = (uint8_t)tx_len;
    memcpy(M_prime + 2, tx, tx_len);
    memcpy(M_prime + 2 + tx_len, M, M_len);

    int result = ml_dsa_verify_internal(pk, M_prime, M_prime_len, sigma);
    free(M_prime);
    return result;
}