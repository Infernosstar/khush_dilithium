#include <stdint.h>
#include "externalfunctions.h"
#include "params.h"
#include "Conversions.h"
#include "fips202.h"
#include <string.h>
#include <stdio.h>

#include <stdint.h>
#include "externalfunctions.h"
#include "params.h"
#include "Conversions.h"
#include "fips202.h"
#include <string.h>
#include <math.h>


int CoeffFromHalfByte(int b) {
    if (b < 0 || b > 15) {
        return INVALID_COEFFICIENT;
    }

    if (ETA == 2) {
        if (b < 15) {
            return 2 - (b % 5);  
        }
    }
    else if (ETA == 4) {
        if (b < 9) {
            return 4 - b; 
        }
    }
    
    return INVALID_COEFFICIENT;
}
void expandA(uint8_t rho[SEEDBYTES], int32_t (*A)[ML_L][256]) {
    uint8_t rho_dash[SEEDBYTES + 2];
    memcpy(rho_dash, rho, SEEDBYTES);
    
    for(uint8_t r = 0; r < ML_K; r++) {
        for(uint8_t s = 0; s < ML_L; s++) {
             IntegerToBytes(s,1,rho_dash + SEEDBYTES);
            IntegerToBytes(r,1,rho_dash + SEEDBYTES+1);
            RejNttPoly(rho_dash, A[r][s]);
        }
    }
}

int32_t CoeffFromThreeBytes(uint8_t b0, uint8_t b1, uint8_t b2) {
    uint32_t b2_prime = b2;
    if (b2 > 127) {
        b2_prime = b2 - 128;
    }
    uint32_t z = (b2_prime << 16) | (b1 << 8) | b0;
    return (z < Q) ? (int32_t)z : -1; // Explicit cast to int32_t
}

void RejNttPoly(uint8_t rho[SEEDBYTES+2], int32_t *a) {
    keccak_state ctx;
    shake128_init(&ctx);
    shake128_absorb(&ctx, rho, SEEDBYTES+2);
    uint8_t s[3];
    int j=0;
    while (j<256){
        shake128_squeeze(s, 3, &ctx);
        int32_t coeff = CoeffFromThreeBytes(s[0], s[1], s[2]);
        if (coeff != -1 && coeff != 0) {
            a[j] = coeff;
            j++;  // Accept this coefficient and move to next
        }
    }
}

void expandS(uint8_t rho[(SEEDBYTES<<1)+2], int32_t *s1, int32_t *s2) {
    uint8_t extended_seed[(SEEDBYTES<<1) + 2]; // Main seed + 2 bytes for counter
    
    // Copy main seed
    memcpy(extended_seed, rho, SEEDBYTES<<1);
    
    // Generate s1 (ℓ polynomials)
    for(int r = 0; r < ML_L; r++) {
        // Append 2-byte counter (little-endian)
        extended_seed[SEEDBYTES<<1] = r & 0xFF;
        extended_seed[(SEEDBYTES<<1)+1] = (r >> 8) & 0xFF;
        
        RejBoundedPoly(extended_seed, &s1[r*256]);
    }
    
    // Generate s2 (k polynomials)
    for(int r = 0; r < ML_K; r++) {
        // Append 2-byte counter (r + ℓ, little-endian)
        uint16_t counter = r + ML_L;
        extended_seed[SEEDBYTES<<1] = counter & 0xFF;
        extended_seed[(SEEDBYTES<<1)+1] = (counter >> 8) & 0xFF;
        
        RejBoundedPoly(extended_seed, &s2[r*256]);
    }
}

void RejBoundedPoly(uint8_t rho[SEEDBYTES+2], int32_t *poly) {
    keccak_state ctx;
    
    // Initialize hash context with extended seed (ρ||r)
    shake256_init(&ctx);
    shake256_absorb(&ctx, rho, SEEDBYTES+2);
    shake256_finalize(&ctx);
    
    int j = 0;
    while(j < 256) {
        uint8_t z;
        shake256_squeeze(&z, 1, &ctx);
        
        // Process both half-bytes
        int z0 = z & 0x0F;        // Lower 4 bits
        int z1 = (z >> 4) & 0x0F;  // Upper 4 bits
        
        int coeff;
        
        // Process first half-byte
        coeff = CoeffFromHalfByte(z0);
        if(coeff != INVALID_COEFFICIENT) {
            poly[j++] = coeff;
            if(j == 256) break;
        }
        
        // Process second half-byte
        coeff = CoeffFromHalfByte(z1);
        if(coeff != INVALID_COEFFICIENT) {
            poly[j++] = coeff;
            if(j == 256) break;
        }
    }
}
static const int32_t zetas[N] = {
         0,    25847, -2608894,  -518909,   237124,  -777960,  -876248,   466468,
   1826347,  2353451,  -359251, -2091905,  3119733, -2884855,  3111497,  2680103,
   2725464,  1024112, -1079900,  3585928,  -549488, -1119584,  2619752, -2108549,
  -2118186, -3859737, -1399561, -3277672,  1757237,   -19422,  4010497,   280005,
   2706023,    95776,  3077325,  3530437, -1661693, -3592148, -2537516,  3915439,
  -3861115, -3043716,  3574422, -2867647,  3539968,  -300467,  2348700,  -539299,
  -1699267, -1643818,  3505694, -3821735,  3507263, -2140649, -1600420,  3699596,
    811944,   531354,   954230,  3881043,  3900724, -2556880,  2071892, -2797779,
  -3930395, -1528703, -3677745, -3041255, -1452451,  3475950,  2176455, -1585221,
  -1257611,  1939314, -4083598, -1000202, -3190144, -3157330, -3632928,   126922,
   3412210,  -983419,  2147896,  2715295, -2967645, -3693493,  -411027, -2477047,
   -671102, -1228525,   -22981, -1308169,  -381987,  1349076,  1852771, -1430430,
  -3343383,   264944,   508951,  3097992,    44288, -1100098,   904516,  3958618,
  -3724342,    -8578,  1653064, -3249728,  2389356,  -210977,   759969, -1316856,
    189548, -3553272,  3159746, -1851402, -2409325,  -177440,  1315589,  1341330,
   1285669, -1584928,  -812732, -1439742, -3019102, -3881060, -3628969,  3839961,
   2091667,  3407706,  2316500,  3817976, -3342478,  2244091, -2446433, -3562462,
    266997,  2434439, -1235728,  3513181, -3520352, -3759364, -1197226, -3193378,
    900702,  1859098,   909542,   819034,   495491, -1613174,   -43260,  -522500,
   -655327, -3122442,  2031748,  3207046, -3556995,  -525098,  -768622, -3595838,
    342297,   286988, -2437823,  4108315,  3437287, -3342277,  1735879,   203044,
   2842341,  2691481, -2590150,  1265009,  4055324,  1247620,  2486353,  1595974,
  -3767016,  1250494,  2635921, -3548272, -2994039,  1869119,  1903435, -1050970,
  -1333058,  1237275, -3318210, -1430225,  -451100,  1312455,  3306115, -1962642,
  -1279661,  1917081, -2546312, -1374803,  1500165,   777191,  2235880,  3406031,
   -542412, -2831860, -1671176, -1846953, -2584293, -3724270,   594136, -3776993,
  -2013608,  2432395,  2454455,  -164721,  1957272,  3369112,   185531, -1207385,
  -3183426,   162844,  1616392,  3014001,   810149,  1652634, -3694233, -1799107,
  -3038916,  3523897,  3866901,   269760,  2213111,  -975884,  1717735,   472078,
   -426683,  1723600, -1803090,  1910376, -1667432, -1104333,  -260646, -3833893,
  -2939036, -2235985,  -420899, -2286327,   183443,  -976891,  1612842, -3545687,
   -554416,  3919660,   -48306, -1362209,  3937738,  1400424,  -846154,  1976782
};
// Montgomery reduction helper function
int32_t montgomery_reduce(int64_t a) {
    int32_t t;
    t = (int32_t)a * 58728449; // Q^-1 mod 2^32
    t = (a - (int64_t)t * Q) >> 32;
    return (int32_t)t;
}
// Improved NTT with correct butterfly operations
void NTT(int32_t a[256]) {
    int len, start, j, k;
    int32_t zeta, t;

    k = 1;
    for (len = 128; len >= 1; len >>= 1) {
        for (start = 0; start < 256; start = j + len) {
            zeta = zetas[k++];
            for (j = start; j < start + len; j++) {
                t = montgomery_reduce((int64_t)zeta * a[j + len]);
                a[j + len] = a[j] - t;
                a[j] = a[j] + t;
                
                // Modular reduction
                a[j] = (a[j] < 0) ? a[j] + Q : (a[j] >= Q) ? a[j] - Q : a[j];
                a[j + len] = (a[j + len] < 0) ? a[j + len] + Q : 
                                (a[j + len] >= Q) ? a[j + len] - Q : a[j + len];
            }
        }
    }
}

// Improved Inverse NTT
void NTT_Inverse(int32_t a[256]) {
    int len, start, j, k;
    int32_t zeta, t;

    k = 256;
    for (len = 1; len < 256; len <<= 1) {
        for (start = 0; start < 256; start = j + len) {
            zeta = zetas[--k];
            for (j = start; j < start + len; j++) {
                t = a[j];
                a[j] = t + a[j + len];
                a[j + len] = t - a[j + len];
                a[j + len] = montgomery_reduce((int64_t)zeta * a[j + len]);
                
                // Modular reduction
                a[j] = (a[j] < 0) ? a[j] + Q : (a[j] >= Q) ? a[j] - Q : a[j];
            }
        }
    }

    // Final scaling by N^-1
    const int32_t f = 8347681; // 256^-1 mod Q
    for (j = 0; j < 256; j++) {
        a[j] = montgomery_reduce((int64_t)a[j] * f);
    }
}



int bitlen(int b) {
    return (b == 0) ? 1 : (int)log2(b) + 1;
}


int SimpleBitPack(uint8_t *output, int32_t w[256], int b)
{
    int bitlen = 0;
    int temp = b;
    while (temp > 0)
    {
        bitlen++;
        temp >>= 1;
    }

    uint8_t bits[256 * bitlen];
    for (int i = 0; i < 256; i++)
    {
        for (int j = 0; j < bitlen; j++)
        {
            bits[i * bitlen + j] = (w[i] >> j) & 1;
        }
    }
    BitsToBytes(output, bits, 256 * bitlen);
    return 0;
}

// Algorithm 35: Power2Round
void Power2Round(int32_t r, int32_t *r1, int32_t *r0) {
    int32_t r_plus = r % Q;
    *r0 = r_plus % (1 << D);
    if (*r0 > (1 << (D-1))) {
        *r0 -= (1 << D);
    }
    *r1 = (r_plus - *r0) / (1 << D);
}


void pkEncode(uint8_t *pk, uint8_t rho[SEEDBYTES], int32_t t1[ML_K][256])
{
    memcpy(pk, rho, SEEDBYTES);
    int offset = SEEDBYTES;

    for (int i = 0; i < ML_K; i++)
    {
        SimpleBitPack(pk + offset, t1[i], (1 << (BITLEN_Q1 - D)) - 1);
        offset += 32 * (BITLEN_Q1 - D);
    }
}


// Modified BitPack for signed coefficients (used in skEncode)
void SignedBitPack(uint8_t *output, const int32_t w[256], int eta, int max_coeff) {
    int bit_length = bitlen(2 * eta);
    uint8_t bit_string[256 * bit_length];
    memset(bit_string, 0, sizeof(bit_string));
    
    for (int i = 0; i < 256; i++) {
        int32_t coeff = w[i];
        // Center lift to [0, 2*eta] or [0, 2^(d-1)]
        if (coeff < 0) coeff += max_coeff + 1;
        
        // Convert to bits
        for (int j = 0; j < bit_length; j++) {
            bit_string[i * bit_length + j] = (coeff >> j) & 1;
        }
    }
    
    // Convert bits to bytes
    int byte_length = (256 * bit_length + 7) / 8;
    memset(output, 0, byte_length);
    for (int i = 0; i < 256 * bit_length; i++) {
        if (bit_string[i]) {
            output[i / 8] |= (1 << (i % 8));
        }
    }
}

// skEncode algorithm
void skEncode(uint8_t *sk, uint8_t rho[SEEDBYTES], uint8_t K[SEEDBYTES], 
              uint8_t tr[CRHBYTES], int32_t s1[ML_L][256], 
              int32_t s2[ML_K][256], int32_t t0[ML_K][256]) {
    // Copy ρ, K, tr
    memcpy(sk, rho, SEEDBYTES);
    memcpy(sk + SEEDBYTES, K, SEEDBYTES);
    memcpy(sk + 2*SEEDBYTES, tr, CRHBYTES);
    size_t offset = 2*SEEDBYTES + CRHBYTES;
    
    // Pack s1 (ℓ polynomials)
    for (int i = 0; i < ML_L; i++) {
        SignedBitPack(sk + offset, s1[i], ETA, 2*ETA);
        offset += 32 * bitlen(2*ETA);
    }
    
    // Pack s2 (k polynomials)
    for (int i = 0; i < ML_K; i++) {
        SignedBitPack(sk + offset, s2[i], ETA, 2*ETA);
        offset += 32 * bitlen(2*ETA);
    }
    
    // Pack t0 (k polynomials)
    for (int i = 0; i < ML_K; i++) {
        SignedBitPack(sk + offset, t0[i], (1 << (D-1)) - 1, 1 << D);
        offset += 32 * D;
    }
}