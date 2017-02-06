/* This is free and unencumbered software released into the public domain.
 *
 * Anyone is free to copy, modify, publish, use, compile, sell, or
 * distribute this software, either in source code form or as a compiled
 * binary, for any purpose, commercial or non-commercial, and by any
 * means.
 *
 * See LICENSE for complete information.
 */


#include "tesla_utils.h"
#include "tesla_rand_openssl_aes.h"


/**************************************************************PRINTING AND COPYING******************************************************/

void copy_poly(RINGELT f[N], const RINGELT g[N]){
  size_t i;
  for(i = 0; i < N; i++) f[i] = g[i];
}

void print_poly(const RINGELT f[N]){
  size_t i;
  for(i = 0; i < N ; i++){
    printf("%ld ", 2*f[i] < q_ringelt ? f[i] : f[i] - q_ringelt );
  }
  printf("\n");
}

void print_sk(const signing_key_t sk){
  printf("s:");
  print_poly(sk.s);
  printf("\n");
  printf("e1:");
  print_poly(sk.e1);
  printf("\n");
  printf("e2:");
  print_poly(sk.e2);
}

void print_pk(const public_key_t pk){
  printf("t1:");
  print_poly(pk.t1);
  printf("\n");
  printf("t2:");
  print_poly(pk.t2);
}

void print_hash(const unsigned char hash[TESLA_DIGEST_LENGTH]){
  size_t i;
  for(i = 0; i < TESLA_DIGEST_LENGTH; i++){
    printf("%02x",hash[i]);
  }
   printf("\n");
}

void print_sig(const signature_t sig){
  printf("z:");
  print_poly(sig.z);
  printf("\n");
  printf("c_dash:");
  print_hash(sig.c_dash);
  printf("\n");
}


/****************************************************ROUNDING FUNCTIONS**********************************/


void rounded_poly2bytes(unsigned char *y, RINGELT f_rd[N]){
  RINGELT x;
  uint16_t i,j;
  for(i = 0; i < N; i++){
    x = f_rd[i];
#if !NISPOWEROFTWO 
    if (i == N-1)
      x = 0;
#endif
    for(j = 0; j < ROUNDED_POLY_BYTES; j++){
      y[ROUNDED_POLY_BYTES*i + j] = x & 0xff;
      x >>= 8;
    }
  }
}

void round_poly(RINGELT f_rd[N], RINGELT f[N]){
  uint16_t i;
  for(i = 0; i < N; i++)
    f_rd[i] = f[i] >> D;
}


/*pass in x modulo q in [0,q) */
/*let r equal x modulo 2^d in (-2^(d-1), 2^(d-1)] */
/*return absolute value of r*/
RINGELT abs_low_bits(RINGELT x){
  int64_t out = x & ~(~((int64_t)0) << D);
  out -= (1<<(D-1));
  out += 1;
  return (RINGELT)abs(out);
}


/***************************************************************HASH ROUTINE*************************************/

/*hash function */
/*input: two polynomials, mu (usually itself a message digest),
         and the length of mu in bytes*/
/*output: a 256-bit hash */

int hash(unsigned char hash_output[TESLA_DIGEST_LENGTH],
	  RINGELT u[N],
	  RINGELT v[N],
	  const unsigned char* mu,
	  const size_t mulen){
  RINGELT u_rd[N], v_rd[N];
  size_t bytesPerRoundedPoly = ROUNDED_POLY_BYTES * TESLA_N;
  int ret;

  /*round u and v coefficients */
  round_poly(u_rd,u);
  round_poly(v_rd,v);

  uint64_t hash_input_bytes = 2*bytesPerRoundedPoly + mulen;
  size_t i;

  /*set hash input*/
  unsigned char *hash_input;
  if ((hash_input = (unsigned char *)malloc(hash_input_bytes)) == NULL) {
    return 0;
  }    

  rounded_poly2bytes(hash_input,u_rd);
  rounded_poly2bytes(&(hash_input[bytesPerRoundedPoly]),v_rd); 
  for(i = 0; i < mulen; i++) {
    hash_input[i + 2*bytesPerRoundedPoly] = mu[i];
  }


  /*compute hash*/
  SHA256_CTX sha256;
  ret = 1;
  ret &= SHA256_Init(&sha256);
  ret &= SHA256_Update(&sha256, hash_input, hash_input_bytes);
  ret &= SHA256_Final(hash_output, &sha256);
  
  /*clean up*/
  free(hash_input);
  return ret;

}

  
/***************************************************************SPARSE POLY MULTIPLICATION*************************/

/*input: a (regular),b (sparse, i.e. list of non-zero co-ordinates) in physical space*/
/*sets v = a*b,  */
 void sparse_mul(RINGELT v[N], const RINGELT a[N], const uint16_t b_sparse[OMEGA], const uint16_t b_sparse_signs[OMEGA]){
  RINGELT v_aux[2*N];
  uint16_t i,j;
  /*zero the output*/
  for(i = 0; i < 2*N; i++) v_aux[i] = 0;

  /*multiply in Z[x]*/
  for(i = 0; i < OMEGA; i++) {
    for(j = 0; j < N; j++) {
      if(b_sparse_signs[i] == 1) {
        ADD_MOD(v_aux[b_sparse[i] + j], v_aux[b_sparse[i] + j], a[j],Q);
      }
      else {
        SUB_MOD(v_aux[b_sparse[i] + j], v_aux[b_sparse[i] + j], a[j],Q);   
      }
    }
  }
  
  
#if NISPOWEROFTWO
  /*reduce mod x^n + 1*/
  for(i = 0; i < N; i++) SUB_MOD(v[i], v_aux[i], v_aux[i+N], Q);
#else
  /*reduce mod x^n - 1*/
  for(i = 0; i < N; i++) ADD_MOD(v[i], v_aux[i], v_aux[i+N], Q);
#endif
}

/**************************************************************ENCODE ROUTINES*************************************/


 int encode_sparse(uint16_t encode_output[OMEGA],
		    uint16_t encode_signs[OMEGA],
		    const unsigned char hash_output[TESLA_DIGEST_LENGTH]){ 
  /*key AES on hash output*/
  AES_KEY aes_key;
    
  if (AES_set_encrypt_key(hash_output,128,&aes_key) < 0)
    return 0;

  /*initialise AES */
  unsigned char aes_ivec[AES_BLOCK_SIZE];
  memset(aes_ivec, 0, AES_BLOCK_SIZE); 
  unsigned char aes_ecount_buf[AES_BLOCK_SIZE]; 
  memset(aes_ecount_buf, 0, AES_BLOCK_SIZE); 
  unsigned int aes_num = 0; 
  unsigned char aes_in[AES_BLOCK_SIZE]; 
  memset(aes_in, 0, AES_BLOCK_SIZE);

  /*get OMEGA values in [0,n), each with a 0 or 1 to indicate sign*/
  uint16_t pos,sign,i,j;
  for(i=0; i< OMEGA; i++){
    while(1){
      AES_ctr128_encrypt(aes_in, (unsigned char *) &pos, 8, &aes_key,
			 aes_ivec, aes_ecount_buf, &aes_num);
      pos &= ~((~0)<<(NBITS + 1));
      sign = pos&1;
      pos <<= 1;

#if NISPOWEROFTWO
      if(pos < N)
#else
      if(pos < N-1)
#endif
        {
          /*check we're not using this position already */
          int success=1;
          for(j = 0; j < i; j++){
            if(pos == j) success = 0;
          }
          if(success) break;
        }

    }
    encode_signs[i] = sign;
    encode_output[i] = pos;
  } 
  return 1;
}


/*compares two ring elements - called in qsort in checkE */
int cmp_checkE(const void *a, const void *b){
  if(*(RINGELT*)a == *(RINGELT*)b) return 0;
  if(*(RINGELT*)a >  *(RINGELT*)b) return -1;
  return 1;
}

/*checks e, which needs to be recorded in physical space */
int checkE(RINGELT e[N]){
  RINGELT f[N];
  uint64_t sum = 0;
  copy_poly(f,e);
  size_t i;
  for(i=0; i < N; i++){
    f[i] = ABS(f[i]);
  }
  qsort(&(f[0]),N,sizeof(*f),cmp_checkE);
  for(i = 0; i < OMEGA ; i++) sum += f[i];
  return (sum <= L);
}

