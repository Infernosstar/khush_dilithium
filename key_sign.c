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
void print_Hex(const uint8_t *data, size_t size)
{
    for (size_t i = 0; i < size; i++)
    {
        printf("%02X", data[i]);
    }
    printf("\n");
}
static int compute_bitlen(int value)
{
    int bitlen = 0;
    while (value > 0)
    {
        bitlen++;
        value >>= 1;
    }
    return bitlen;
}

void BitUnpack(int32_t *poly, const uint8_t *bytes, int b, int eta)
{
    int bitlen = compute_bitlen(b);
    int offset = 0;

    for (int i = 0; i < 256; i++)
    {
        poly[i] = 0;
        for (int j = 0; j < bitlen; j++)
        {
            if (bytes[offset / 8] & (1 << (offset % 8)))
            {
                poly[i] |= (1 << j);
            }
            offset++;
        }

        if (poly[i] > eta)
        {
            poly[i] -= 2 * eta + 1;
        }
    }
}

static void NTT_Multiply(const uint8_t c[32], const int32_t poly_hat[256], int32_t *result)
{
    for (int i = 0; i < 256; i++)
    {
        result[i] = (poly_hat[i] * c[i % 32]) % Q;
    }
}

void skDecode(const uint8_t *sk,
              uint8_t rho[SEEDBYTES],
              uint8_t K[SEEDBYTES],
              uint8_t tr[CRHBYTES],
              int32_t s1[ML_L][256],
              int32_t s2[ML_K][256],
              int32_t t0[ML_K][256])
{

    memcpy(rho, sk, SEEDBYTES);
    memcpy(K, sk + SEEDBYTES, SEEDBYTES);
    memcpy(tr, sk + 2 * SEEDBYTES, CRHBYTES);

    int offset = 2 * SEEDBYTES + CRHBYTES;
    int bitlen_2eta = compute_bitlen(2 * ETA);
    int bytes_per_poly_2eta = (256 * bitlen_2eta + 7) / 8;

    // Unpack s1 (L polynomials)
    for (int i = 0; i < ML_L; i++)
    {
        BitUnpack(s1[i], sk + offset, 2 * ETA, ETA);
        offset += bytes_per_poly_2eta;
    }

    // Unpack s2 (K polynomials)
    for (int i = 0; i < ML_K; i++)
    {
        BitUnpack(s2[i], sk + offset, 2 * ETA, ETA);
        offset += bytes_per_poly_2eta;
    }

    // Unpack t0 (K polynomials)
    int bitlen_t0 = D;
    int bytes_per_poly_t0 = (256 * bitlen_t0 + 7) / 8;
    for (int i = 0; i < ML_K; i++)
    {
        BitUnpack(t0[i], sk + offset, (1 << (D - 1)) - 1, (1 << (D - 1)) - 1);
        offset += bytes_per_poly_t0;
    }
}

void Decompose(uint32_t r, uint32_t *q1, uint32_t *q2)
{
    uint32_t r_plus = r % Q;
    *q1 = r_plus % (2 * GAMMA2);
    if (r_plus - *q1 == Q - 1)
    {
        *q2 = 0;
        *q1 = *q1 - 1;
    }
    else
    {
        *q2 = (r_plus - *q1) / (2 * GAMMA2);
    }
}

void HighBits(uint32_t r, uint32_t *r0)
{
    uint32_t r1;
    Decompose(r, &r1, r0);
}

void w1Encode(uint8_t *w1_encoded, const int32_t w1[ML_K][256])
{
    int32_t max_coeff = (Q - 1) / (2 * GAMMA2) - 1;
    int bitlen = compute_bitlen(max_coeff);
    size_t bytes_per_poly = (256 * bitlen + 7) / 8;

    for (int i = 0; i < ML_K; i++)
    {
        SimpleBitPack(w1_encoded + i * bytes_per_poly, w1[i], max_coeff);
    }
}

void expandMask(uint8_t rho[SEEDBYTES + 2], int32_t u, int32_t y[ML_L][256])
{
    uint8_t c = 1;
    int32_t temp = GAMMA1 - 1;
    while (temp >>= 1)
        c++;

    uint8_t v[32 * c];
    for (int r = 0; r < ML_L; r++)
    {
        IntegerToBytes(u + r, 2, rho + SEEDBYTES);
        shake256(v, 32 * c, rho, SEEDBYTES + 2);
        BitUnpack(y[r], v, GAMMA1 - 1, GAMMA1);
    }
}

void sampleinBall(uint8_t *c, uint8_t rho[SEEDBYTES])
{
    keccak_state ctx;
    uint8_t s[8];
    uint8_t j;

    shake256_init(&ctx);
    shake256_absorb(&ctx, rho, SEEDBYTES);
    shake256_squeeze(s, 8, &ctx);

    for (int i = 256 - TU; i < 256; i++)
    {
        do
        {
            shake256_squeeze(&j, 1, &ctx);
        } while (j > i);
        c[i] = c[j];
        c[j] = 1;
    }
}

void LowBits(uint32_t r, uint32_t *r1)
{
    uint32_t r0;
    Decompose(r, r1, &r0);
}

int makeHint(uint32_t z, uint32_t r)
{
    uint32_t r1, v1;
    HighBits(r, &r1);
    HighBits(z + r, &v1);
    return (v1 != r1) ? 1 : 0;
}

void BitPack(uint8_t *output, const int32_t w[256], int32_t a, int32_t b)
{
    int32_t range = a + b;
    int bitlen = compute_bitlen(range);
    size_t total_bits = 256 * bitlen;
    size_t total_bytes = (total_bits + 7) / 8;

    memset(output, 0, total_bytes);

    for (int i = 0; i < 256; i++)
    {
        int32_t val = b - w[i];
        val = (val < 0) ? 0 : ((val > range) ? range : val);

        for (int j = 0; j < bitlen; j++)
        {
            int bit_pos = i * bitlen + j;
            if (val & (1 << j))
            {
                output[bit_pos / 8] |= (1 << (bit_pos % 8));
            }
        }
    }
}

void HintBitPack(int32_t h[ML_K][256], uint8_t *y)
{
    int Index = 0;
    for (int i = 0; i < ML_K; i++)
    {
        for (int j = 0; j < 256; j++)
        {
            if (h[i][j] != 0)
            {
                y[Index] = j;
                Index++;
            }
        }
        y[OMEGA + i] = Index;
    }
}

void sigEncode(uint8_t *C_hash, int32_t z[ML_L][256], int32_t h[ML_K][256], uint8_t *sigma)
{
    size_t size = 32;
    memcpy(sigma, C_hash, size);
    uint8_t *y = sigma + size;

    for (int i = 0; i < ML_L; i++)
    {
        BitPack(y, z[i], GAMMA1 - 1, GAMMA1);
        y += (256 * (BITLEN_GAMMA1 * 2) + 7) / 8;
    }

    HintBitPack(h, y);
}
int32_t mul_mod(int32_t a, int32_t b, int32_t q)
{
    return (int32_t)a * b % q;
}
void pointwise_mul(int32_t *res, const int32_t *a, const int32_t *b)
{
    for (int i = 0; i < N; i++)
    {
        res[i] = mul_mod(a[i], b[i], Q);
    }
}
int sign_internal(const uint8_t *sk, const uint8_t *M_dash, size_t M_len,
                  const uint8_t rnd[32], uint8_t *sigma)
{
    // Step 1: Decode private key
    uint8_t rho[SEEDBYTES], K[SEEDBYTES], tr[CRHBYTES];
    int32_t s1[ML_L][256], s2[ML_K][256], t0[ML_K][256];
    skDecode(sk, rho, K, tr, s1, s2, t0);

    // Step 2: Convert s1, s2, t0 to NTT domain
    for (int i = 0; i < ML_L; i++)
        NTT(s1[i]);
    for (int i = 0; i < ML_K; i++)
    {
        NTT(s2[i]);
        NTT(t0[i]);
    }

    // Step 3: Expand matrix A in NTT domain
    int32_t A_hat[ML_K][ML_L][256];
    expandA(rho, A_hat);

    // Step 4: Compute μ = H(tr || M′)
    uint8_t mu[64];
    {
        uint8_t tr_bits[CRHBYTES * 8];
        BytesToBits(tr, CRHBYTES, tr_bits);
        uint8_t input[CRHBYTES * 8 + M_len];
        memcpy(input, tr_bits, CRHBYTES * 8);
        memcpy(input + CRHBYTES * 8, M_dash, M_len);
        shake256(mu, 64, input, CRHBYTES * 8 + M_len);
    }

    // Step 5: Compute ρ′′ = H(K || rnd || μ)
    uint8_t rho_prime_prime[64];
    {
        uint8_t input[SEEDBYTES + 32 + 64];
        memcpy(input, K, SEEDBYTES);
        memcpy(input + SEEDBYTES, rnd, 32);
        memcpy(input + SEEDBYTES + 32, mu, 64);
        shake256(rho_prime_prime, 64, input, SEEDBYTES + 32 + 64);
    }

    // Step 6: Rejection sampling loop
    uint8_t kappa = 0;
    int32_t z[ML_L][256], h[ML_K][256];
    uint8_t c[32];
    int32_t y[ML_L][256], w[ML_K][256], w1[ML_K][256];
    uint8_t w1_encoded[ML_K * ((256 * (BITLEN_GAMMA1 * 2) + 7) / 8)];
int o=0;
    while (1 && o<1)
    {
        o++;
        // Step 7: Expand mask to get y
        expandMask(rho_prime_prime, kappa, y);

        // Step 8: Compute w = NTT^-1(A_hat ◦ NTT(y))
        for (int i = 0; i < ML_K; i++)
        {
            for (int j = 0; j < ML_L; j++)
            {
                int32_t temp[256];
                NTT(y[j]); // Convert y to NTT domain
                pointwise_mul(temp, A_hat[i][j], y[j]);
                NTT_Inverse(temp); // Convert back to normal domain
                if (j == 0)
                {
                    memcpy(w[i], temp, 256 * sizeof(int32_t));
                }
                else
                {
                    printf("printing and generating w :\n");
                    for (int k = 0; k < 256; k++)
                    {
                        w[i][k] = (w[i][k] + temp[k]) % Q;
                        printf("%d ",w[i][k]);
                    }
                    printf("\n");
                }
            }
        }

        // Step 9: Compute w1 = HighBits(w)
        for (int i = 0; i < ML_K; i++)
        {
            printf("printing the highbits values:\n");
            for (int j = 0; j < 256; j++)
            {
                HighBits(w[i][j], &w1[i][j]);
                printf("%d ",w[i][j]);
            
            }
            printf("\n");
        }

        // Step 10: Compute c = H(μ || w1_encoded)
        w1Encode(w1_encoded, w1);
        //this part is just for fun 
        for(int i =0;i<5120;i++){
            printf("%d ",w1_encoded[i]);
        }
        {
            uint8_t input[64 + sizeof(w1_encoded)];
            memcpy(input, mu, 64);
            memcpy(input + 64, w1_encoded, sizeof(w1_encoded));
            shake256(c, 32, input, 64 + sizeof(w1_encoded));
        }

        // Step 11: Convert c to NTT domain
        int32_t c_ntt[256];
        sampleinBall(c, rho_prime_prime); // Ensure c is in proper format
        NTT(c_ntt);

        // Step 12: Compute ⟨⟨cs1⟩⟩ and ⟨⟨cs2⟩⟩
        int32_t cs1[ML_L][256], cs2[ML_K][256];
        for (int i = 0; i < ML_L; i++)
        {
            pointwise_mul(cs1[i], c_ntt, s1[i]);
            NTT_Inverse(cs1[i]);
        }
        for (int i = 0; i < ML_K; i++)
        {
            pointwise_mul(cs2[i], c_ntt, s2[i]);
            NTT_Inverse(cs2[i]);
        }

        // Step 13: Compute z = y + cs1
        for (int i = 0; i < ML_L; i++)
        {
            for (int j = 0; j < 256; j++)
            {
                z[i][j] = y[i][j] + cs1[i][j];
                // Apply mod± reduction
                if (z[i][j] > GAMMA1 - BETA)
                    z[i][j] -= Q;
                else if (z[i][j] < -(GAMMA1 - BETA))
                    z[i][j] += Q;
            }
        }
       uint8_t pk[PKBYTES], sk[SKBYTES];
/*
        // Generate keys
        keygen(pk, sk);

        // Print keys
        printf("Public Key:\n");
        print_Hex(pk, PKBYTES);

        printf("\nSecret Key:\n");
        print_Hex(sk, SKBYTES);
*/
       // return 0;
        // Step 14: Compute r0 = LowBits(w - cs2)
        int32_t r0[ML_K][256];
        int valid = 1;
        for (int i = 0; i < ML_K; i++)
        {
            for (int j = 0; j < 256; j++)
            {
                int32_t temp = (w[i][j] - cs2[i][j]) % Q;
                LowBits(temp, &r0[i][j]);
                // Check validity condition
                if (abs(r0[i][j]) >= GAMMA2 - BETA)
                {
                    valid = 0;
              
                    break;
                }
            }
            if (!valid)
                break;
        }

        if (!valid)
        {
            kappa += ML_L;
            continue;
        }

        // Step 15: Compute ⟨⟨ct0⟩⟩
        int32_t ct0[ML_K][256];
        for (int i = 0; i < ML_K; i++)
        {
            pointwise_mul(ct0[i], c_ntt, t0[i]);
            NTT_Inverse(ct0[i]);
        }

        // Step 16: Make hints
        uint32_t hint[ML_K][256];
        int hint_count = 0;
        for (int i = 0; i < ML_K; i++)
        {
            for (int j = 0; j < 256; j++)
            {
                // Check ct0 validity condition
                if (abs(ct0[i][j]) >= GAMMA2)
                {
                    valid = 0;
                    break;
                }
                hint[i][j] = makeHint(-ct0[i][j],
                                      (w[i][j] - cs2[i][j] + ct0[i][j]) % Q);
                hint_count += hint[i][j];
            }
            if (!valid)
                break;
        }
        printf("valid statement:%d\n",valid);
        printf("hint_count=%d\n",hint_count);
        if (!valid || hint_count > OMEGA)
        {
            kappa += ML_L;
            continue;
        }

        // If we get here, we have a valid signature
        break;
    }

    // Step 17: Encode signature
    sigEncode(c, z, h, sigma);
    return 0;
}
int ml_dsa_sign(const uint8_t *sk, const uint8_t *M, size_t M_len, 
               const uint8_t *ctx, size_t ctx_len, uint8_t *sigma) {
    // Step 1: Check context length
    if (ctx_len > 255) {
        return -1; // Error: context too long
    }

    // Step 2: Generate random bytes (for randomized variant)
    uint8_t rnd[32];
    randombytes(rnd, 32);
    
    // For deterministic variant, use all zeros instead:
    // memset(rnd, 0, 32);

    // Step 3: Format M′
    size_t M_prime_len = 2 + ctx_len + M_len;
    uint8_t *M_prime = (uint8_t *)malloc(M_prime_len);
    if (!M_prime) return -1;

    // First byte: 0 (as integer)
    M_prime[0] = 0;
    // Second byte: ctx length
    M_prime[1] = (uint8_t)ctx_len;
    // Context bytes
    memcpy(M_prime + 2, ctx, ctx_len);
    // Message bytes
    memcpy(M_prime + 2 + ctx_len, M, M_len);

    // Step 4: Call internal signing function
    int ret = sign_internal(sk, M_prime, M_prime_len, rnd, sigma);

    free(M_prime);
    return ret;
}
