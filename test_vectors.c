/* This is free and unencumbered software released into the public domain.
 *
 * Anyone is free to copy, modify, publish, use, compile, sell, or
 * distribute this software, either in source code form or as a compiled
 * binary, for any purpose, commercial or non-commercial, and by any
 * means.
 *
 * See LICENSE for complete information.
 */


#include "test_vectors.h"

int test_vector(test_vec_t vec){
  int success = 1;

  /*check that public key is derived from signing key*/
  public_key_t pk_test;
  gen_pk(&pk_test, vec.sk);
  for(int i = 0; i < N; i++){
    if((pk_test.t1)[i] != (vec.pk.t1)[i] || (pk_test.t2)[i] != (vec.pk.t2)[i]){
      printf("failed to derive public key from signing key\n");
      success = 0;
      break;
    }
  }
  
  /*check that supplied signature verifies correctly*/
  if(!verify(vec.sig, vec.pk, (unsigned char*)(vec.message), strlen(vec.message))
     ){
      printf("failed to verify test vector with supplied signatures!\n");
      success = 0;
    }

  /*check that computed signature matches provided one*/
  signature_t sig_test;
  deterministic_sign(&sig_test, vec.y, vec.sk, (unsigned char *)vec.message, strlen(vec.message));
  for(int i = 0; i < N; i++){
    if((sig_test.z)[i] != (vec.sig.z)[i]){
      printf("failed to produce signature from test vector!\n");
      success = 0;
      break;
    }
    
  }
  for(int i = 0; i < TESLA_DIGEST_LENGTH; i++){
    if((sig_test.c_dash)[i] != (vec.sig.c_dash)[i]){
      printf("failed to produce signature from test vector!\n");
      success = 0;
      break;
    }
    
  }

  return success;

}

