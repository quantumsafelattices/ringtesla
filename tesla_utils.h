/* This is free and unencumbered software released into the public domain.
 *
 * Anyone is free to copy, modify, publish, use, compile, sell, or
 * distribute this software, either in source code form or as a compiled
 * binary, for any purpose, commercial or non-commercial, and by any
 * means.
 *
 * See LICENSE for complete information.
 */

#ifndef TESLA_UTILS_H
#define TESLA_UTILS_H

#include "tesla_rand.h"

void copy_poly(RINGELT f[N],
	       const RINGELT g[N]);
void print_poly(const RINGELT f[N]);
void print_pk(const public_key_t pk);
void print_sk(const signing_key_t sk);
void print_hash(const unsigned char hash[TESLA_DIGEST_LENGTH]);
void print_sig(const signature_t sig);
void rounded_poly2bytes(unsigned char *y, RINGELT f_rd[N]);
void round_poly(RINGELT f_rd[N], RINGELT f[N]);
RINGELT abs_low_bits(RINGELT x);

int hash(unsigned char hash_output[TESLA_DIGEST_LENGTH],
	 RINGELT u[N],
	 RINGELT v[N],
	 const unsigned char* mu,
	 const size_t mulen);
void sparse_mul(RINGELT v[N],
		  const RINGELT a[N],
		  const uint16_t b_sparse[OMEGA],
		  const uint16_t b_sparse_signs[OMEGA]);
int encode_sparse(uint16_t encode_output[OMEGA],
		  uint16_t encode_output_signs[OMEGA],
		  const unsigned char hash_output[TESLA_DIGEST_LENGTH]);
void encode(RINGELT encode_output[N],
	    const unsigned char hash_output[TESLA_DIGEST_LENGTH]);
int cmp_checkE(const void *a, const void *b);
int checkE(RINGELT e[N]);

#endif
