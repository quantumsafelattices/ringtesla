/* This is free and unencumbered software released into the public domain.
 *
 * Anyone is free to copy, modify, publish, use, compile, sell, or
 * distribute this software, either in source code form or as a compiled
 * binary, for any purpose, commercial or non-commercial, and by any
 * means.
 *
 * See LICENSE for complete information.
 */

#ifndef HEADER_TESLA_H
#define HEADER_TESLA_H

#include <stdio.h>
#include <inttypes.h>
#include <openssl/sha.h>
#include "tesla_consts.h"
#include "tesla_rand.h"

#define TESLA_DIGEST_LENGTH 32
  
#ifndef RINGELT
#define RINGELT uint32_t
#endif


/*signature structs */
struct tesla_public_key_st{
  RINGELT t1[N];
  RINGELT t2[N];
};

struct tesla_signing_key_st{
  RINGELT s[N];
  RINGELT e1[N];
  RINGELT e2[N];
};

struct tesla_signature_st{
  RINGELT z[N];
  unsigned char c_dash[TESLA_DIGEST_LENGTH];
};


typedef struct tesla_public_key_st public_key_t;
typedef struct tesla_signing_key_st signing_key_t;
typedef struct tesla_signature_st signature_t; 

  /*macros for working mod q */
#define SIGN(A) (0 == (A))? 0 : ((2*(A) <= Q) ? 1 : -1) 
#define ABS(A) (2*(A) <= Q ? (A) : Q - (A))
#define NEG(A) ((SIGN(A) >= 0)? (A) : (Q - (A)))

/* Internal function prototypes */
void gen_sk(signing_key_t *sk);
void gen_pk(public_key_t *pk, signing_key_t sk);
void sign(signature_t *sig,
	  const signing_key_t sk,  
	  const unsigned char *message_digest,
	  const size_t dgst_len);
int verify(const signature_t sig, 
	   const public_key_t pk, 
	   const unsigned char *message_digest,
	   const size_t dgst_len);
int deterministic_sign(signature_t *signature, 
		       const RINGELT y[N],
		       const signing_key_t sk, 
		       const unsigned char *message_digest,
		       const size_t dgst_len);


#endif /* HEADER_TESLA_H */
