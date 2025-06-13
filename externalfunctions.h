#ifndef EXTERNALFUNCTIONS_H
#define EXTERNALFUNCTIONS_H

#include <stdint.h>
#include "params.h"

#define INVALID_COEFFICIENT -1

// Polynomial and NTT operations
void NTT(int32_t w[256]);
void NTT_Inverse(int32_t w_cap[256]);
void Power2Round(int32_t r, int32_t *r1, int32_t *r0);

// Expansion functions
void expandA(uint8_t rho[SEEDBYTES], int32_t (*A)[ML_L][256]);
void expandS(uint8_t rho[(SEEDBYTES<<1)+2], int32_t *s1, int32_t *s2);

// Encoding functions
void pkEncode(uint8_t *pk, uint8_t rho[SEEDBYTES], int32_t t1[ML_K][256]);
void skEncode(uint8_t *sk, uint8_t rho[SEEDBYTES], uint8_t key[SEEDBYTES], 
              uint8_t tr[CRHBYTES], int32_t s1[ML_L][256], 
              int32_t s2[ML_K][256], int32_t t0[ML_K][256]);

// Helper functions
int32_t CoeffFromThreeBytes(uint8_t b0, uint8_t b1, uint8_t b2);
void RejNttPoly(uint8_t rho[SEEDBYTES+2], int32_t *a);
void RejBoundedPoly(uint8_t rho[SEEDBYTES+2], int32_t *s0);
void IntegerToBytes(int x, int alpha, uint8_t *y);
int SimpleBitPack(uint8_t *output, int32_t w[256], int b);

#endif