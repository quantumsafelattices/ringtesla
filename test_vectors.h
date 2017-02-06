/* This is free and unencumbered software released into the public domain.
 *
 * Anyone is free to copy, modify, publish, use, compile, sell, or
 * distribute this software, either in source code form or as a compiled
 * binary, for any purpose, commercial or non-commercial, and by any
 * means.
 *
 * See LICENSE for complete information.
 */

#ifndef TESLA_TEST_VECTORS_H
#define TESLA_TEST_VECTORS_H

#include "tesla.h"

struct test_vec_st{
  public_key_t pk;
  signing_key_t sk;
  RINGELT y[N];
  signature_t sig;
  char *message;
};

typedef struct test_vec_st test_vec_t;


#if TESLA_N == 1024
#include "test_vectors_1024.h"
#endif

#if TESLA_N == 709
#include "test_vectors_709.h"
#endif

#if TESLA_N == 827
#include "test_vectors_827.h"
#endif


int test_vector(test_vec_t vec);

#endif /* TESLA_TEST_VECTORS_H */
