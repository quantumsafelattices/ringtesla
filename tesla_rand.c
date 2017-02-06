/* This is free and unencumbered software released into the public domain.
 *
 * Anyone is free to copy, modify, publish, use, compile, sell, or
 * distribute this software, either in source code form or as a compiled
 * binary, for any purpose, commercial or non-commercial, and by any
 * means.
 *
 * See LICENSE for complete information.
 */


#include "tesla_rand.h"
#include "tesla_rand_openssl_aes.h"

#if TESLA_N == 709
#include "table_stdev28.h"
#endif

#if TESLA_N == 827
#include "table_stdev23.h"
#endif

#if TESLA_N == 1024
#include "table_stdev2.h"
#endif

#ifndef LOWER_INDEX_START
#error tesla_rand: invalid TESLA_N
#endif


#define RANDOM192(c) c[0] = RANDOM64; c[1] = RANDOM64; c[2] = RANDOM64

static uint64_t ct_isnonzero_u64(uint64_t x);
static uint64_t ct_lt_u64(uint64_t x, uint64_t y);
static uint64_t ct_mask_u64(uint64_t bit);
static uint64_t ct_ne_u64(uint64_t x, uint64_t y);
static uint64_t ct_eq_u64(uint64_t x, uint64_t y);
static int cmplt_ct(uint64_t *a, uint64_t *b);
static uint64_t ct_select_u64(uint64_t x, uint64_t y, uint64_t bit);

static RINGELT single_sample_ct(uint64_t *in);



/*
 * Returns 1 if x != 0
 * Returns 0 if x == 0
 * x and y are arbitrary unsigned 64-bit integers
 */
static uint64_t ct_isnonzero_u64(uint64_t x) {
    return (x|-x) >> 63;
}


/* Returns 1 if x < y
 * Returns 0 if x >= y
 * x and y are arbitrary unsigned 64-bit integers
 */
static uint64_t ct_lt_u64(uint64_t x, uint64_t y) {
    return (x^((x^y)|((x-y)^y))) >> 63;
}

/* Returns 0xFFFF..FFFF if bit != 0
 * Returns            0 if bit == 0
 */
static uint64_t ct_mask_u64(uint64_t bit)
{
    return 0 - (uint64_t)ct_isnonzero_u64(bit);
}

/*
 * Returns 1 if x != y
 * Returns 0 if x == y
 * x and y are arbitrary unsigned 64-bit integers
 */
static uint64_t ct_ne_u64(uint64_t x, uint64_t y) {
    return ((x-y)|(y-x)) >> 63;
}

/*
 * Returns 1 if x == y
 * Returns 0 if x != y
 * x and y are arbitrary unsigned 64-bit integers
 */
static uint64_t ct_eq_u64(uint64_t x, uint64_t y) {
    return 1 ^ ct_ne_u64(x, y);
}



/* Returns 0 if a >= b
 * Returns 1 if a < b
 * Where a and b are both 3-limb 64-bit integers.
 * This function runs in constant time.
 */
static int cmplt_ct(uint64_t *a, uint64_t *b) {

	uint64_t r = 0; /* result */
	uint64_t _m = 0; /* mask   */
	int i;
	for(i = 2; i >= 0; i--) {
		r |= ct_lt_u64(a[i], b[i]) & ~_m;
		_m |= ct_mask_u64(ct_ne_u64(a[i], b[i])); /* stop when a[i] != b[i] */
	}
	return r & 1;
}

/* Conditionally return x or y depending on whether bit is set
 * Equivalent to: return bit ? x : y
 * x and y are arbitrary 64-bit unsigned integers
 * bit must be either 0 or 1.
 */
static uint64_t ct_select_u64(uint64_t x, uint64_t y, uint64_t bit) {
    uint64_t _m = ct_mask_u64(bit);
    return (x&_m) | (y&~_m);
}


/***********************************************DISCRETE GAUSSIAN SAMPLING FROM  TABLE ****************************/

static RINGELT single_sample_ct(uint64_t *in) {
	uint16_t lower_index = LOWER_INDEX_START, this_index = THIS_INDEX_START, upper_index = UPPER_INDEX_START;
	size_t i;
	for (i = 0; i <= UPPER_INDEX_START_LOG2; i++) {
		if (cmplt_ct(in, tesla_table[this_index])) {
			upper_index = this_index;
		} else {
			lower_index = this_index;
		}
		this_index = (lower_index + upper_index) / 2;		
	}
	return upper_index;
}

void sample_gaussian_ct(RINGELT f[N]){
  	RANDOM_VARS
	  uint64_t r=0, rnd[3]={0,0,0};
	size_t i;
#if NISPOWEROFTWO
	for(i = 0; i < N; i++) {
#else
        f[N-1] = 0;
        for(i = 0; i < N -1 ; i++) {
#endif
          if (!( i &0x3F)) r = RANDOM64;
          RANDOM192(rnd);		
          f[i] = single_sample_ct(rnd);
          RINGELT t = Q - f[i];
          f[i] = ct_select_u64(t, f[i], ct_eq_u64(r&1, 1) & ct_ne_u64(f[i], 0));
          r >>= 1;
	}

}
