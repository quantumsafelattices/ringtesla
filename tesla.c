/* This is free and unencumbered software released into the public domain.
 *
 * Anyone is free to copy, modify, publish, use, compile, sell, or
 * distribute this software, either in source code form or as a compiled
 * binary, for any purpose, commercial or non-commercial, and by any
 * means.
 *
 * See LICENSE for complete information.
 */

#include <string.h>

#include "tesla.h"
#include "tesla_rand_openssl_aes.h"
#include "tesla_utils.h"

/* Ring-TESLA primitive functions */


/***********************************************SIGNING KEY GENERATION********************************/

/*generates signing key (s,e1,e2), stored in physical form */
void gen_sk(signing_key_t *sk){
    
  /*loop through sampling s,e1,e2 until e1,e2 pass the relevant checks*/
  while(1){
    sample_gaussian_ct(sk->e1);
    sample_gaussian_ct(sk->e2);
    if(checkE(sk->e1) && (checkE(sk->e2))){
      sample_gaussian_ct(sk->s);
      break;
    }
  }
}


/***********************************************PUBLIC KEY GENERATION********************************/

/*takes a signing key stored in physical space and computes the public key in physical space */
/*points a1, a2 are stored in FFT space */
void gen_pk(public_key_t *pk, signing_key_t sk){
  FFT_FORWARD(sk.s);
  FFT_FORWARD(sk.e1);
  FFT_FORWARD(sk.e2);
  POINTWISE_MUL_ADD(pk->t1, a1, sk.s, sk.e1, N, Q);
  POINTWISE_MUL_ADD(pk->t2, a2, sk.s, sk.e2, N, Q);  
  FFT_BACKWARD(sk.s);
  FFT_BACKWARD(sk.e1);
  FFT_BACKWARD(sk.e2);
  FFT_BACKWARD(pk->t1);
  FFT_BACKWARD(pk->t2);
}


/**********************************************SIGNING*********************************************/

#if NISPOWEROFTWO
#define N_LIMIT N
#else
#define N_LIMIT N-1
#endif

/*signs a message as (z,c_dash) where z is a ring elt in physical form, and c_dash is a hash output */
void sign(signature_t *sig,
          const signing_key_t sk,  
          const unsigned char *message_digest,
	  const size_t dgst_len){
  RINGELT y[N]={0};
  int success = 0;
  size_t i;
 
  /*initialise random sampling for y*/
  RANDOM_VARS
  /*sample y randomly, and repeat until it passes rejection sampling*/
  while(!success){
    for(i = 0; i < N_LIMIT; i++){
#if !NISPOWEROFTWO
    y[N-1] = 0;
#endif
      while(1){
        y[i] = RANDOM32; /*get 32 bits of random */
        y[i] &= ~(~0 << (B_BITS + 1)); /*take bottom (B_BITS + 1) bits */
        if (y[i] <= 2*B + 1) break;
      }
      y[i] = (y[i] <= B)? y[i] : Q - (y[i] - B); 
    }

    success = deterministic_sign(sig, y, sk, message_digest, dgst_len);
  }
}


/*signs a message for a fixed choice of ephemeral secret y in physcial space
 returns 1 or 0 according to success or failure in doing so (due to rejection sampling)*/
int deterministic_sign(signature_t *signature, 
                       const RINGELT y[N],
		       const signing_key_t sk, 
		       const unsigned char *message_digest,
		       const size_t dgst_len){
  RINGELT v1[N], v2[N], y_fft[N];
  uint16_t i;
  uint16_t c_sparse[OMEGA];
  uint16_t c_sparse_signs[OMEGA];
  
  copy_poly(y_fft,y);
  FFT_FORWARD(y_fft);

 /* v_i = a_i*y in FFT space*/
  POINTWISE_MUL(v1,a1,y_fft,N,Q);
  POINTWISE_MUL(v2,a2,y_fft,N,Q);
  FFT_BACKWARD(v1);
  FFT_BACKWARD(v2);

  MAPTOCYCLOTOMIC(v1, N, Q);
  MAPTOCYCLOTOMIC(v2, N, Q);

  /*compute hash in physical space*/
  if (hash(signature->c_dash, v1, v2, message_digest, dgst_len) == 0)
    return 0;

  /*put c_dash through encoding function to obtain sparse polynomial c in physical space*/
  if (encode_sparse(c_sparse, c_sparse_signs, signature->c_dash) == 0)
    return 0;
  /*compute z = y + s*c in physical space*/
  sparse_mul(signature->z, sk.s, c_sparse, c_sparse_signs);
  POINTWISE_ADD(signature->z, signature->z,y,N,Q);
  MAPTOCYCLOTOMIC(signature->z, N, Q);
  /*rejection sampling for security: check that coefficients of z have absolute value bounded by B-U*/
  for(i = 0; i < N_LIMIT; i++ ){
     if(ABS((signature->z)[i]) > (B-U))
     return 0; 
  }
  
  /*rejection sampling for correctness*/
  /*w_i = v_i - e_i c */
  RINGELT w1[N], w2[N], tmp[N];
  sparse_mul(tmp,sk.e1,c_sparse,c_sparse_signs);
  POINTWISE_SUB(w1,v1,tmp,N,Q);
  sparse_mul(tmp,sk.e2,c_sparse,c_sparse_signs);
  POINTWISE_SUB(w2,v2,tmp,N,Q);

  MAPTOCYCLOTOMIC(w1,N,Q);
  MAPTOCYCLOTOMIC(w2,N,Q);
  for(i = 0; i < N_LIMIT; i++ ){
    if(abs_low_bits(w1[i]) > (1<<(D-1))-L )
      return 0;
    if(abs_low_bits(w2[i]) > (1<<(D-1))-L )
      return 0;
  }

  /*all rejection sampling passed*/
#if TEST_VEC_VERBOSE
  printf("\ny used in signature:\n");
  print_poly(y);
  printf("\n");
  printf("\ny_fft used in signature:\n");
  print_poly(y_fft);
  printf("\n");
#endif

  return 1;

}



/**********************************************VERIFICATION*********************************************/

int verify(const signature_t sig, 
           const public_key_t pk, 
           const unsigned char *message_digest,
	   const size_t dgst_len){
  
  RINGELT w1[N], w2[N], temp[N], z_copy[N];
  uint16_t c_sparse[OMEGA];
  uint16_t c_sparse_signs[OMEGA];
  uint16_t i;
  unsigned char c_dash[TESLA_DIGEST_LENGTH];

  if (encode_sparse(c_sparse, c_sparse_signs, sig.c_dash) == 0) 
    return 0;
  
  copy_poly(z_copy, sig.z);
  FFT_FORWARD(z_copy);
  POINTWISE_MUL(w1, a1, z_copy, N, Q);
  POINTWISE_MUL(w2, a2, z_copy, N, Q);
  FFT_BACKWARD(w1);
  FFT_BACKWARD(w2);
  sparse_mul(temp,pk.t1,c_sparse,c_sparse_signs);  
  POINTWISE_SUB(w1,w1,temp,N,Q);
  sparse_mul(temp,pk.t2,c_sparse,c_sparse_signs); 
  POINTWISE_SUB(w2,w2,temp,N,Q);
  MAPTOCYCLOTOMIC(w1, N, Q);
  MAPTOCYCLOTOMIC(w2, N, Q);

  if (hash(c_dash, w1, w2, message_digest, dgst_len) == 0)
    return 0;

  for(i = 0; i < TESLA_DIGEST_LENGTH; i++){
    if(c_dash[i] != sig.c_dash[i]) return 0;
  }
  
  return 1;
}

