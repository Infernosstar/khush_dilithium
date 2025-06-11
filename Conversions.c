#include <stdint.h>
#include <stdio.h>
#include<string.h>
/*******************************|
 * | Integer to bits program algo |
 * 
 *  x' ← x
 *  for i from 0 to α − 1 do
 *    y[i] ← x' mod 2
 *    x' ← ⌊x' / 2⌋
 *  end for
 *  return y
 * 
 ********************************/
void IntegerToBits(int x, int alpha, int *y) {
    int x_dash = x;
    for (int i = 0; i < alpha; i++) {
        y[i] = x_dash & 1;
        x_dash = x_dash >> 1;
    }
}

/*******************************
 * Algorithm 10: BitsToInteger
 * Computes the integer value expressed by a bit string using little-endian order.
 * Input: bit string y of length alpha
 * Output: integer x
 ********************************/
int BitsToInteger(const int *y, int alpha) {
    int x = 0;
    for (int i = 1; i <= alpha; i++) {
        x = (x << 1) | y[alpha - i];
    }
    return x;
}

/*******************************
 * Algorithm 11: IntegerToBytes
 * Computes a base-256 representation of x mod 256^alpha using little-endian order.
 * Input: integer x, positive integer alpha
 * Output: byte array y of length alpha
 ********************************/
void IntegerToBytes(int x, int alpha, uint8_t *y) {
    int x_dash = x;
    for (int i = 0; i < alpha; i++) {
        y[i] = x_dash & 0xFF;
        x_dash = x_dash >> 8;
    }
}

/*******************************
 * Algorithm 12: BitsToBytes
 * Converts a bit string y of length alpha to a byte string z of length alpha/8 using little-endian order.
 * Input: bit string y of length alpha
 * Output: byte string z of length alpha/8
 ********************************/

void BitsToBytes(uint8_t *output, const uint8_t *bits, int bitlen) {
    memset(output, 0, bitlen/8 + (bitlen%8 != 0));
    for (int i = 0; i < bitlen; i++) {
        output[i/8] |= (bits[i] & 1) << (i%8);
    }
}

/*******************************
 * Algorithm 13: BytesToBits
 * Converts a byte string z of length alpha to a bit string y of length 8*alpha using little-endian order.
 * Input: byte string z of length alpha
 * Output: bit string y of length 8*alpha
 ********************************/
void BytesToBits(const uint8_t *z, int alpha, uint8_t *y) {
    for (int i = 0; i < alpha; i++) {
        for (int j = 0; j < 8; j++) {
            y[8 * i + j] = (z[i] >> j) & 1;
        }
    }
}


