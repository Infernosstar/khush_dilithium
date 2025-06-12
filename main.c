#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "randombytes.h"
#include "fips202.h"
#include "externalfunctions.h"
#include "params.h"
#include "key_gen.h"
#include "key_sign.h"
#include "verify.h"
/*
void print_Hex(const uint8_t *data, size_t size) {
    for (size_t i = 0; i < size; i++) {
        printf("%02X", data[i]);
    }
    printf("\n");
}*/

int main() {
    // Test parameters
    const char *message = "Hello There!";
    size_t message_len = strlen(message);
    
    // Keccak state for hashing
    keccak_state ctx;
    
    // Allocate buffers
    uint8_t pk[PKBYTES];
    uint8_t sk[SKBYTES];
    uint8_t sigma[SIGBYTES];
    memset(sigma, 0, SIGBYTES);

    // Generate test keypair
    printf("=== Generating ML-DSA Key Pair ===\n");
    keygen(pk, sk);
    
    printf("Secret key (%d bytes):\n", SKBYTES);
    print_Hex(sk, SKBYTES);
    printf("\nPublic key (%d bytes):\n", PKBYTES);
    print_Hex(pk, PKBYTES);

    // Sign the message
    printf("\n=== Signing Message ===\n");
    printf("Message: '%s'\n", message);
    
    // Initialize the hash states
    shake256_init(&ctx);
    
    int ret = ml_dsa_sign(sk, (const uint8_t *)message, message_len, 
                        &ctx, sigma);
    
    if (ret != 0) {
        printf("Signing failed with error %d\n", ret);
        return 1;
    }

    printf("\nSignature (%d bytes):\n", SIGBYTES);
    print_Hex(sigma, SIGBYTES);

    // Basic signature validation
    printf("\n=== Basic Signature Validation ===\n");
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

    // Full verification (reinitialize state for verification)
   /* printf("\n=== Verifying Signature ===\n");
    shake256_init(&ctx);
    int verify_result = ml_dsa_verify(pk, (const uint8_t *)message, message_len, &ctx, sigma);
    
    if (verify_result) {
        printf("Signature verification SUCCESS\n");
    } else {
        printf("Signature verification FAILED\n");
        return 1;
    }
*/
    return 0;
}