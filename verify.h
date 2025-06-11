#ifndef VERIFY_H
#define VERIFY_H

#include <stdint.h>
#include "fips202.h"
#include "params.h"

int ml_dsa_verify(const uint8_t *pk, const uint8_t *M, size_t M_len,
                 keccak_state *state, const uint8_t *sigma);

int ml_dsa_verify_internal(const uint8_t *pk, const uint8_t *M_prime, 
                          size_t M_prime_len, keccak_state *state,
                          const uint8_t *sigma);

#endif