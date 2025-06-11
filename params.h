// params.h
#ifndef PARAMS_H
#define PARAMS_H

#include <stdint.h>

// Security parameters
#define ML_K 4
#define ML_L 4
#define ETA 2
#define D 13
#define GAMMA1 (1 << 17)
#define GAMMA2 ((8380417-1)/88)  // (Q-1)/88
#define BETA 78
#define OMEGA 80
#define TU 64

// Bit lengths
#define BITLEN_2ETA 3
#define BITLEN_GAMMA1 20
#define BITLEN_Q1 23

// Sizes
#define SEEDBYTES 32
#define CRHBYTES 64
#define PKBYTES 1312
#define SKBYTES 2560
#define SIGBYTES 2420

// Constants
#define Q 8380417
#define N 256
#define INVALID_COEFFICIENT -1

#endif