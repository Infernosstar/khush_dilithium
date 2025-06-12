#include <stdint.h>
#include <string.h>
#include "params.h"
#include "fips202.h"
#include "externalfunctions.h"
#include "Conversions.h"
#include <stdio.h>
#include <stdlib.h>
#include "key_sign.h"
#include "key_gen.h"
static void SimpleBitUnpack(int32_t *poly, const uint8_t *bytes, int bitlen) {
    int offset = 0;
    for (int i = 0; i < 256; i++) {
        poly[i] = 0;
        for (int j = 0; j < bitlen; j++) {
            if (bytes[offset / 8] & (1 << (offset % 8))) {
                poly[i] |= (1 << j);
            }
            offset++;
        }
    }
}

// Algorithm 23: pkDecode
void pkDecode(const uint8_t *pk, uint8_t *rho, int32_t t1[ML_K][256]) {
    // Step 1: Extract ρ and z_i components from public key
    memcpy(rho, pk, SEEDBYTES); // ρ is first 32 bytes
    
    // Each t1 polynomial is packed in 32*(bitlen(q-1)-d) bits
    int bitlen = BITLEN_Q_MINUS_1 - D;
    int bytes_per_poly = (256 * bitlen + 7) / 8;
    
    // Step 2-4: Unpack each t1 polynomial
    for (int i = 0; i < ML_K; i++) {
        const uint8_t *z_i = pk + SEEDBYTES + i * bytes_per_poly;
        SimpleBitUnpack(t1[i], z_i, bitlen);
    }
}

// Helper function for signature decoding
static void BitUnpack(int32_t *poly, const uint8_t *bytes, int a, int b) {
    int range = a + b;
    int bitlen = 0;
    while ((1 << bitlen) <= range) bitlen++;
    
    int offset = 0;
    for (int i = 0; i < 256; i++) {
        poly[i] = 0;
        for (int j = 0; j < bitlen; j++) {
            if (bytes[offset / 8] & (1 << (offset % 8))) {
                poly[i] |= (1 << j);
            }
            offset++;
        }
        // Map from [0, a+b] to [-a, b]
        poly[i] -= a;
    }
}

// Helper function to unpack hints
static void HintBitUnpack(const uint8_t *y, int32_t h[ML_K][256]) {
    // Initialize all hints to 0
    for (int i = 0; i < ML_K; i++) {
        for (int j = 0; j < 256; j++) {
            h[i][j] = 0;
        }
    }
    
    int index = 0;
    for (int i = 0; i < ML_K; i++) {
        // Read omega_i (number of non-zero hints for this polynomial)
        int omega_i = y[OMEGA + i];
        
        // Mark the specified positions as 1
        for (int j = index; j < omega_i; j++) {
            int pos = y[j];
            if (pos < 256) {  // Safety check
                h[i][pos] = 1;
            }
        }
        index = omega_i;
    }
}

// Algorithm 27: sigDecode
int sigDecode(const uint8_t *sigma, uint8_t *c_tilde, int32_t z[ML_L][256], int32_t h[ML_K][256]) {
    // Step 1: Extract components from signature
    // c_tilde is first 32 bytes (λ/4 bits where λ=256)
    memcpy(c_tilde, sigma, 32);
    
    // Each z polynomial is packed in 32*(1+bitlen(γ1-1)) bytes
    int bitlen_z = 1 + BITLEN_GAMMA1_MINUS_1;
    int bytes_per_z_poly = (256 * bitlen_z + 7) / 8;
    
    // Step 2-4: Unpack each z polynomial
    for (int i = 0; i < ML_L; i++) {
        const uint8_t *x_i = sigma + 32 + i * bytes_per_z_poly;
        BitUnpack(z[i], x_i, GAMMA1 - 1, GAMMA1);
    }
    
    // Step 5: Unpack hints
    const uint8_t *y = sigma + 32 + ML_L * bytes_per_z_poly;
    HintBitUnpack(y, h);
    
    // Check if hints were properly encoded
    int hint_count = 0;
    for (int i = 0; i < ML_K; i++) {
        for (int j = 0; j < 256; j++) {
            hint_count += h[i][j];
        }
    }
    
    // If hint count exceeds ω, return error
    if (hint_count > OMEGA) {
        return -1;
    }
    
    return 0;
}
// Function to verify a signature
   int ml_dsa_verify_internal(const uint8_t *pk, const uint8_t *M_prime, size_t M_prime_len, const uint8_t *sigma) {
    // Step 1: Decode public key
  int ml_dsa_verify_internal(const uint8_t *pk, const uint8_t *M_prime,
                          size_t M_prime_len, keccak_state *state,
                          const uint8_t *sigma) {
    // Step 1: Decode public key
    uint8_t rho[SEEDBYTES];
    int32_t t1[ML_K][256];
    pkDecode(pk, rho, t1, state);

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
    uint8_t w1_encoded[ML_K * ((256 * (BITLEN_GAMMA1 * 2) + 7) / 8)];
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

int ml_dsa_verify(const uint8_t *pk, const uint8_t *M, size_t M_len,
                 keccak_state *state, const uint8_t *sigma) {
    // ML-DSA doesn't use tx/tx_len - format M′ as (0 || M)
    size_t M_prime_len = 1 + M_len;
    uint8_t *M_prime = (uint8_t *)malloc(M_prime_len);
    if (!M_prime) return 0;

    M_prime[0] = 0;  // Single byte prefix
    memcpy(M_prime + 1, M, M_len);

    int result = ml_dsa_verify_internal(pk, M_prime, M_prime_len, state, sigma);
    free(M_prime);
    return result;
}}