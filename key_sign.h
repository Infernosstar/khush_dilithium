#include<stdint.h>
#include "params.h"
#include "fips202.h"
#include "randombytes.h"

int key_sign_internal(const uint8_t *sk,
                          const uint8_t *M_prime,
                          size_t M_prime_len,
                          const uint8_t rnd[32],
                          uint8_t *sigma);
void BitUnpack(int32_t *poly, const uint8_t *bytes, int b, int eta);
void NTT_Multiply(const uint8_t c[32], const int32_t poly_hat[256], int32_t *result);
// void skDecode(uint8_t *sk,uint8_t *rho,uint8_t *K,uint8_t *tr);
void skDecode(const uint8_t *sk,
              uint8_t rho[SEEDBYTES],
              uint8_t K[SEEDBYTES],
              uint8_t tr[CRHBYTES],
              int32_t s1[ML_L][256],
              int32_t s2[ML_K][256],
              int32_t t0[ML_K][256]);
void BitUnpack(int32_t *poly, const uint8_t *bytes, int b, int eta);
void Decompose(uint32_t r, uint32_t *q1, uint32_t *q2);
void HighBits(uint32_t r, uint32_t *r0);
void w1Encode(uint8_t *w1_encoded, const int32_t w1[ML_K][256]);
void expandMask(uint8_t rho[66], int32_t u, int32_t y[ML_L][256]);
void sampleinBall(uint8_t *c, uint8_t rho[64]);
void LowBits(uint32_t r, uint32_t *r1);

int makeHint(uint32_t z, uint32_t r);
void BitPack(uint8_t *output, const int32_t w[256], int32_t a, int32_t b);
void HintBitPack(uint32_t h[ML_K][256], uint8_t *y);
void sigEncode(uint8_t *C_hash, int32_t z[ML_L][256], uint32_t h[ML_K][256], uint8_t *sigma);
int key_sign_internal(const uint8_t *sk,
                      const uint8_t *M_prime,
                      size_t M_prime_len,
                      const uint8_t rnd[32],
                      uint8_t *sigma);
int ml_dsa_sign(const uint8_t *sk, const uint8_t *M, size_t M_len,
                const uint8_t *ctx, size_t ctx_len, uint8_t *sigma);
void print_Hex(const uint8_t *data, size_t size);