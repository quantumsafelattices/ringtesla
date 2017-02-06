/* This is free and unencumbered software released into the public domain.
 *
 * Anyone is free to copy, modify, publish, use, compile, sell, or
 * distribute this software, either in source code form or as a compiled
 * binary, for any purpose, commercial or non-commercial, and by any
 * means.
 *
 * See LICENSE for complete information.
 */


#include "tesla.h"
#include "tesla_utils.h"
#include "test_vectors.h"
#define SIGN_TRIALS 1000


int main(){

 signing_key_t sk;
  public_key_t pk;
  char *message = "testtest";
  signature_t sig;

  /*print a single example*/
  printf("example signature");
  printf("\nmessage:\n");
  printf("%s\n",message);
  gen_sk(&sk);
  gen_pk(&pk, sk);
  printf("\nsigning key:\n");
  print_sk(sk);
  printf("\npublic key:\n");
  print_pk(pk);
  sign(&sig, sk,(unsigned char *)message,strlen(message));
  printf("\nsignature:\n");
  print_sig(sig);
  
  
  printf("******************************************************************************************************\n");


  /*test a lot of verifications*/
  printf("trying %d independent sign/verifies\n", SIGN_TRIALS);
  for(uint32_t i=0; i < SIGN_TRIALS; i++){
    gen_sk(&sk);
    gen_pk(&pk, sk);
    sign(&sig, sk,(unsigned char *)message,strlen(message));
    if(!verify(sig, pk,(unsigned char *)message,strlen(message))){
      printf("verification failure round %d!\n",i);
      return 1;
    }
    if(!(i % 100)){
      printf("passed trial %d\n",i);
    }
  }
  printf("signature scheme validates across %d independent trials\n", SIGN_TRIALS);
  printf("******************************************************************************************************\n");

  /*run test vectors*/
  printf("running %d test vectors\n",N_TESLA_TEST_VECS);
  for(uint32_t i=0;i<N_TESLA_TEST_VECS; i++){
    printf("vector %d of %d...",i+1,N_TESLA_TEST_VECS);
    if(test_vector(tesla_test_vecs[i])) printf("success\n");
    else printf("failure\n");
  }


  
  return 0;
}
