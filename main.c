#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "randombytes.h"
#include "fips202.h"
#include "externalfunctions.h"
#include "params.h"
#include "key_gen.h"
#include "key_sign.h"
/*
void print_Hex(const uint8_t *data, size_t size)
{
    for (size_t i = 0; i < size; i++)
    {
        printf("%02X", data[i]);
    }
    printf("\n");
}*/

int main() {
    // Test parameters
    const char *message = "Hello There!";
    const char *context = "Test context";
    size_t message_len = strlen(message);
    size_t context_len = strlen(context);

    // Allocate buffers
    uint8_t pk[PKBYTES];
    uint8_t sk[SKBYTES];
    uint8_t sigma[SIGBYTES];

    // Generate test keypair
    keygen(pk, sk);
    
    printf("=== ML-DSA Signing Test ===\n");
    printf("security key:\n");
    print_Hex( sk, SKBYTES);
    printf("\npublic key:\n");
    print_Hex( pk, PKBYTES);

    // Sign the message
    int ret = ml_dsa_sign(sk, (const uint8_t *)message, message_len, 
                         (const uint8_t *)context, context_len, sigma);
    
    if (ret != 0) {
        printf("Signing failed with error %d\n", ret);
        return 1;
    }

    printf("Signing successful!\n");
    printf("Signature: \n");
    print_Hex( sigma, SIGBYTES);

  
    int all_zero = 1;
    for (size_t i = 0; i < SIGBYTES; i++) {
        if (sigma[i] != 0) {
            all_zero = 0;
            break;
        }
    }

    if (all_zero) {
        printf("ERROR: Signature is all zeros!\n");
        return 1;
    }

    printf("Signature appears valid (basic checks passed)\n");
    return 0;
}