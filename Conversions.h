#ifndef CONVERSIONS_H
#define CONVERSIONS_H

#include <stdint.h>

void IntegerToBits(int x, int alpha, int *y);
int BitsToInteger(const int *y, int alpha);
void IntegerToBytes(int x, int alpha, uint8_t *y);
void BitsToBytes(uint8_t *output, const uint8_t *bits, int bitlen);
void BytesToBits(const uint8_t *z, int alpha, uint8_t *y);

#endif