/* This is free and unencumbered software released into the public domain.
 *
 * Anyone is free to copy, modify, publish, use, compile, sell, or
 * distribute this software, either in source code form or as a compiled
 * binary, for any purpose, commercial or non-commercial, and by any
 * means.
 *
 * See LICENSE for complete information.
 */

#ifndef HEADER_TESLA_CONSTS_H
#define HEADER_TESLA_CONSTS_H

#ifndef RINGELT
#define RINGELT uint_fast32_t
#endif

#define RINGELT_BYTES 4

#undef VALID_N

#if TESLA_N == 709
#include "tesla_consts_709.h"
#define VALID_N
#endif

#if TESLA_N == 827
#include "tesla_consts_827.h"
#define VALID_N
#endif

#if TESLA_N == 1024
#include "tesla_consts_1024.h"
#define VALID_N
#endif

#endif /* HEADER_TESLA_CONSTS_H */
