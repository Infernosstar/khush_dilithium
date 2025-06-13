#include <stdint.h>
#include <string.h>
#include "params.h"
#include "fips202.h"
#include "randombytes.h"
#include "externalfunctions.h"
#include "Conversions.h"
#include <stdlib.h>
#include "key_gen.h"
#include <stdio.h>

void print_Hex(const uint8_t *data, size_t size) {
    for (size_t i = 0; i < size; i++) {
        printf("%02X", data[i]);
    }
    printf("\n");
}

static int compute_bitlen(int value) {
    int bitlen = 0;
    while (value > 0) {
        bitlen++;
        value >>= 1;
    }
    return bitlen;
}

void BitUnpack(int32_t *poly, const uint8_t *bytes, int b, int eta) {
    int bitlen = compute_bitlen(b);
    int offset = 0;

    for (int i = 0; i < 256; i++) {
        poly[i] = 0;
        for (int j = 0; j < bitlen; j++) {
            if (bytes[offset / 8] & (1 << (offset % 8))) {
                poly[i] |= (1 << j);
            }
            offset++;
        }

        if (poly[i] > eta) {
            poly[i] -= 2 * eta + 1;
        }
    }
}

static void NTT_Multiply(const uint8_t c[256], const int32_t poly_hat[256], int32_t *result) {
    for (int i = 0; i < 256; i++) {
        result[i] = (poly_hat[i] * c[i % 32]) % Q;
    }
}

void skDecode(const uint8_t *sk,
              uint8_t rho[SEEDBYTES],
              uint8_t K[SEEDBYTES],
              uint8_t tr[CRHBYTES],
              int32_t s1[ML_L][256],
              int32_t s2[ML_K][256],
              int32_t t0[ML_K][256]) {

    memcpy(rho, sk, SEEDBYTES);
    memcpy(K, sk + SEEDBYTES, SEEDBYTES);
    memcpy(tr, sk + 2 * SEEDBYTES, CRHBYTES);

    int offset = 2 * SEEDBYTES + CRHBYTES;
    int bitlen_2eta = compute_bitlen(2 * ETA);
    int bytes_per_poly_2eta = (256 * bitlen_2eta + 7) / 8;

    for (int i = 0; i < ML_L; i++) {
        BitUnpack(s1[i], sk + offset, 2 * ETA, ETA);
        offset += bytes_per_poly_2eta;
    }

    for (int i = 0; i < ML_K; i++) {
        BitUnpack(s2[i], sk + offset, 2 * ETA, ETA);
        offset += bytes_per_poly_2eta;
    }

    int bitlen_t0 = D;
    int bytes_per_poly_t0 = (256 * bitlen_t0 + 7) / 8;
    for (int i = 0; i < ML_K; i++) {
        BitUnpack(t0[i], sk + offset, (1 << (D - 1)) - 1, (1 << (D - 1)) - 1);
        offset += bytes_per_poly_t0;
    }
}

void Decompose(uint32_t r, uint32_t *q1, uint32_t *q2) {
    uint32_t r_plus = r % Q;
    *q1 = r_plus % (2 * GAMMA2);
    if (r_plus - *q1 == Q - 1) {
        *q2 = 0;
        *q1 = *q1 - 1;
    } else {
        *q2 = (r_plus - *q1) / (2 * GAMMA2);
    }
}

void HighBits(uint32_t r, uint32_t *r0) {
    uint32_t r1;
    Decompose(r, &r1, r0);
}

void w1Encode(uint8_t *w1_encoded, const int32_t w1[ML_K][256]) {
    int32_t max_coeff = (Q - 1) / (2 * GAMMA2) - 1;
    int bitlen = compute_bitlen(max_coeff);
    size_t bytes_per_poly = (256 * bitlen + 7) / 8;

    for (int i = 0; i < ML_K; i++) {
        SimpleBitPack(w1_encoded + i * bytes_per_poly, (int32_t *)w1[i], max_coeff);
    }
}

void expandMask(uint8_t rho[SEEDBYTES + 2], int32_t u, int32_t y[ML_L][256]) {
    uint8_t c = 1;
    int32_t temp = GAMMA1 - 1;
    while (temp >>= 1) c++;

    keccak_state state;
    uint8_t v[32 * c];
    
    for (int r = 0; r < ML_L; r++) {
        IntegerToBytes(u + r, 2, rho + SEEDBYTES);
        
        // Proper SHAKE256 sponge usage
        shake256_init(&state);
        shake256_absorb(&state, rho, SEEDBYTES + 2);
        shake256_finalize(&state);
        shake256_squeeze(v, 32 * c, &state);
        
        BitUnpack(y[r], v, GAMMA1 - 1, GAMMA1);
    }
}

void sampleinBall(uint8_t *c, uint8_t rho[SEEDBYTES]) {
    keccak_state state;
    uint8_t s[8];
    uint8_t j;

    // Proper SHAKE256 sponge usage
    shake256_init(&state);
    shake256_absorb(&state, rho, SEEDBYTES);
    shake256_finalize(&state);
    shake256_squeeze(s, 8, &state);

    memset(c, 0, 256);
    for (int i = 256 - TU; i < 256; i++) {
        do {
            shake256_squeeze(&j, 1, &state);
        } while (j > i);
        c[i] = c[j];
        c[j] = 1;
    }
}

void LowBits(uint32_t r, uint32_t *r1) {
    uint32_t r0;
    Decompose(r, r1, &r0);
}

int makeHint(uint32_t z, uint32_t r) {
    uint32_t r1, v1;
    HighBits(r, &r1);
    HighBits(z + r, &v1);
    return (v1 != r1) ? 1 : 0;
}

void BitPack(uint8_t *output, const int32_t w[256], int32_t a, int32_t b) {
    int32_t range = a + b;
    int bitlen = compute_bitlen(range);
    size_t total_bits = 256 * bitlen;
    size_t total_bytes = (total_bits + 7) / 8;

    memset(output, 0, total_bytes);

    for (int i = 0; i < 256; i++) {
        int32_t val = b - w[i];
        val = (val < 0) ? 0 : ((val > range) ? range : val);

        for (int j = 0; j < bitlen; j++) {
            int bit_pos = i * bitlen + j;
            if (val & (1 << j)) {
                output[bit_pos / 8] |= (1 << (bit_pos % 8));
            }
        }
    }
}

void HintBitPack(int32_t h[ML_K][256], uint8_t *y) {
    int Index = 0;
    for (int i = 0; i < ML_K; i++) {
        for (int j = 0; j < 256; j++) {
            if (h[i][j] != 0) {
                y[Index] = j;
                Index++;
            }
        }
        y[OMEGA + i] = Index;
    }
}

void sigEncode(uint8_t *C_hash, int32_t z[ML_L][256], int32_t h[ML_K][256], uint8_t *sigma) {
    size_t size = 32;
    memcpy(sigma, C_hash, size);
    uint8_t *y = sigma + size;

    for (int i = 0; i < ML_L; i++) {
        BitPack(y, z[i], GAMMA1 - 1, GAMMA1);
        y += (256 * (BITLEN_GAMMA1 * 2) + 7) / 8;
    }

    HintBitPack(h, y);
}

int32_t mul_mod(int32_t a, int32_t b, int32_t q) {
    return (int32_t)a * b % q;
}

void pointwise_mul(int32_t *res, const int32_t *a, const int32_t *b) {
    for (int i = 0; i < N; i++) {
        res[i] = mul_mod(a[i], b[i], Q);
    }
}

int sign_internal(const uint8_t *sk, const uint8_t *M_dash, size_t M_len,
                  const uint8_t rnd[32], keccak_state *ctx, uint8_t *sigma) {
    // Decode private key
    uint8_t rho[SEEDBYTES], K[SEEDBYTES], tr[CRHBYTES];
    int32_t s1[ML_L][256], s2[ML_K][256], t0[ML_K][256];
    skDecode(sk, rho, K, tr, s1, s2, t0);

    // Convert to NTT domain
    for (int i = 0; i < ML_L; i++) NTT(s1[i]);
    for (int i = 0; i < ML_K; i++) {
        NTT(s2[i]);
        NTT(t0[i]);
    }

    // Expand matrix A in NTT domain
    int32_t A_hat[ML_K][ML_L][256];
    expandA(rho, A_hat);

    // Compute μ = H(tr || M′) using proper sponge
    uint8_t mu[64];
    {
        // uint8_t input[CRHBYTES + M_len + 1]; // +1 for safe padding
        // memcpy(input, tr, CRHBYTES);
        // memcpy(input + CRHBYTES, M_dash, M_len);
        uint8_t *temp;
        BytesToBits(tr,64,temp);
        shake256_init(ctx);
        shake256_absorb(ctx,temp, CRHBYTES);
        shake256_absorb(ctx, M_dash, M_len);
        shake256_finalize(ctx);
        shake256_squeeze(mu, 64, ctx);
    }

    // Compute ρ′′ = H(K || rnd || μ) with proper sponge
    uint8_t rho_prime_prime[64];
    {
        // uint8_t input[SEEDBYTES + 32 + 64];
        // memcpy(input, K, SEEDBYTES);
        // memcpy(input + SEEDBYTES, rnd, 32);
        // memcpy(input + SEEDBYTES + 32, mu, 64);
        
        shake256_init(ctx);
        shake256_absorb(ctx, K,SEEDBYTES);
        shake256_absorb(ctx, rnd, SEEDBYTES);
        shake256_absorb(ctx,mu,CRHBYTES);
        shake256_finalize(ctx);
        shake256_squeeze(rho_prime_prime, CRHBYTES, ctx);
    }

    uint8_t kappa = 0;
    int32_t z[ML_L][256], h[ML_K][256];
    uint8_t c[32];
    int32_t y[ML_L][256], w[ML_K][256], w1[ML_K][256];
    uint8_t w1_encoded[ML_K * ((256 * (BITLEN_GAMMA1 * 2) + 7) / 8)];
    int attempts = 0;
    const long long int max_attempts = 1024;

    while (attempts < max_attempts) {
        attempts++;
        
        expandMask(rho_prime_prime, kappa, y);

        memset(w, 0, sizeof(w));
        for (int i = 0; i < ML_K; i++) {
            for (int j = 0; j < ML_L; j++) {
                int32_t temp[256];
                NTT(y[j]);
                pointwise_mul(temp, A_hat[i][j], y[j]);
                NTT_Inverse(temp);
                
                for (int k = 0; k < 256; k++) {
                    w[i][k] = (w[i][k] + temp[k]) % Q;
                    if (w[i][k] < 0) w[i][k] += Q;
                }
            }
        }

        for (int i = 0; i < ML_K; i++) {
            for (int j = 0; j < 256; j++) {
                HighBits(w[i][j],(uint32_t *)&w1[i][j]);
            }
        }

        w1Encode(w1_encoded, w1);
        {
            uint8_t input[64 + sizeof(w1_encoded)];
            /*
            memcpy(input, mu, 64);
            memcpy(input + 64, w1_encoded, sizeof(w1_encoded));
            */
            shake256_init(ctx);
            shake256_absorb(ctx, mu,64);
            shake256_absorb(ctx, w1_encoded, sizeof(w1_encoded));
            shake256_finalize(ctx);
            shake256_squeeze(c, 32, ctx);
        }

        sampleinBall(c, rho_prime_prime);

        int valid = 1;
        int hint_count = 0;
        
        for (int i = 0; i < ML_L; i++) {
            for (int j = 0; j < 256; j++) {
                z[i][j] = y[i][j] + s1[i][j];
                if (z[i][j] < -(GAMMA1 - BETA) || z[i][j] > GAMMA1 - BETA) {
                    valid = 0;
                    break;
                }
            }
            if (!valid) break;
        }

        if (!valid) {
            kappa += ML_L;
            continue;
        }

        for (int i = 0; i < ML_K; i++) {
            for (int j = 0; j < 256; j++) {
                int32_t temp = (w[i][j] - s2[i][j]) % Q;
                if (temp < 0) temp += Q;
                uint32_t r0;
                LowBits(temp, &r0);
                if (r0 >= GAMMA2 - BETA) {
                    valid = 0;
                    break;
                }
            }
            if (!valid) break;
        }

        if (!valid) {
            kappa += ML_L;
            continue;
        }

        for (int i = 0; i < ML_K; i++) {
            for (int j = 0; j < 256; j++) {
                h[i][j] = makeHint(t0[i][j], (w[i][j] - s2[i][j] + t0[i][j]) % Q);
                hint_count += h[i][j];
            }
        }

        if (hint_count > OMEGA) {
            kappa += ML_L;
            continue;
        }

        break;
    }

    if (attempts >= max_attempts) {
        return -1;
    }

    sigEncode(c, z, h, sigma);
    return 0;
}

int ml_dsa_sign(const uint8_t *sk, const uint8_t *M, size_t M_len, 
                keccak_state *ctx, uint8_t *sigma) {
    // printf("\nsigning function started");
    uint8_t rnd[32];
    randombytes(rnd, 32);
    
    size_t M_prime_len = 1 + M_len;
    uint8_t *M_prime ;
    // if (!M_prime)return -1;

    M_prime[0] = 0;
    memcpy(M_prime + 1, M, M_len);

    int ret = sign_internal(sk, M_prime, M_prime_len, rnd, ctx, sigma);

    free(M_prime);
    // printf("key signing returned successfully\n");
    // fflush(stdout);
    return ret;
}