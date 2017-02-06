
/* Code to compute a Number Theoretic Transform for multiplication in the ring 
F_q[x] / <x^n+1>.
	n = 709, q = 10486111	*/

#include <stdio.h>
#include "FFT_includes.h"
#include "FFT_709_10486111_constants.h"
#include "FFT_2048_475404289_constants.h"
#include "FFT_2048_950808577_constants.h"

/*
We use Bluestein's trick and integer convoution by CRT.
*/
void FFT_forward_709_10486111(FFTSHORT x[709]) {
        const FFTSHORT n = 709;
	const FFTSHORT q = 10486111;
	const FFTSHORT N = 2048;
  
	FFTSHORT x0[2048], x1[2048];
	FFTSHORT i;
	FFTLONG x_crt[2048];
	const FFTSHORT q0 = 475404289, Ninvq0 = 475172158;	
	const FFTSHORT q1 = 950808577, Ninvq1 = 950344315;
	const FFTLONG h0 = 950808578UL, h1 = 950808577UL;
	const FFTLONG q0q1 = 452018475523786753UL;
	FFTLONG tmp_long;

	/*Setup Bluestein's method*/
	for (i = 0; i < n; ++i) {
		MUL_MOD(x0[i], x[i], Bluestein_mul_709_10486111[i], q);
	}
	memset((void *) (x0+n), 0, (N-n)*sizeof(FFTSHORT)); /*Pad with 0's*/
	memcpy((void *) x1, (void *) x0, N*sizeof(FFTSHORT)); /*Copy x0 into x1*/
	
	/*Cyclic convolution*/
	FFT_forward_2048_475404289(x0);
	FFT_forward_2048_950808577(x1);
	
	for (i = 0; i < N; ++i) {
		MUL_MOD(x0[i], x0[i], Bluestein_roots_fft_2048_475404289[i], q0);
		MUL_MOD(x1[i], x1[i], Bluestein_roots_fft_2048_950808577[i], q1);		
	}
	FFT_backward_2048_475404289(x0);
	FFT_backward_2048_950808577(x1);
	
	/*Apply the CRT*/
	for (i = 0; i < N; ++i) {
		MUL_MOD(x0[i], x0[i], Ninvq0, q0); /*Scaling for the convolution*/
		MUL_MOD(x1[i], x1[i], Ninvq1, q1); /*Scaling for the convolution*/
                MUL_MOD_LONG(x_crt[i], h1, x0[i], q0q1);
	        x_crt[i] = q0q1 - x_crt[i];
                MUL_MOD_LONG(tmp_long, h0, x1[i], q0q1);
                ADD_MOD_LONG(x_crt[i], x_crt[i], tmp_long, q0q1);
		x_crt[i] = x_crt[i] % q0q1; //Will now be the integer convolution
		x_crt[i] = x_crt[i] % q; 
	}
		
	/*Complete Bluestein's trick*/
	x[0] = (FFTSHORT) x_crt[(N>>1)-1];
	for (i = 0; i < n-1; ++i) {
		MUL_MOD(x[i+1], x_crt[(N>>1)+i], Bluestein_mul_709_10486111[i], q);
	}		
	
}


void FFT_backward_709_10486111(FFTSHORT x[709]) {
        const FFTSHORT n = 709;
	const FFTSHORT q = 10486111;
	const FFTSHORT N = 2048;
	
	FFTSHORT x0[2048], x1[2048];
	FFTSHORT i;
	FFTLONG x_crt[2048];
	const FFTSHORT q0 = 475404289, Ninvq0 = 475172158;	
	const FFTSHORT q1 = 950808577, Ninvq1 = 950344315;
	const FFTLONG h0 = 950808578UL, h1 = 950808577UL;
	const FFTLONG q0q1 = 452018475523786753UL;
	
	/*Setup Bluestein's method*/
	for (i = 0; i < n; ++i) {
		MUL_MOD(x0[i], x[i], Bluestein_mul_inv_709_10486111[i], q);
	}
	memset((void *) (x0+n), 0, (N-n)*sizeof(FFTSHORT)); /*Pad with 0's*/
	memcpy((void *) x1, (void *) x0, N*sizeof(FFTSHORT)); /*Copy x0 into x1*/
		
	/*Cyclic convolution*/
	FFT_forward_2048_475404289(x0);
	FFT_forward_2048_950808577(x1);
	for (i = 0; i < N; ++i) {
		MUL_MOD(x0[i], x0[i], Bluestein_roots_inv_fft_2048_475404289[i], q0);
		MUL_MOD(x1[i], x1[i], Bluestein_roots_inv_fft_2048_950808577[i], q1);		
	}
	FFT_backward_2048_475404289(x0);
	FFT_backward_2048_950808577(x1);
	
	/*Apply the CRT*/
        FFTLONG tmp_long;
	for (i = 0; i < N; ++i) {
		MUL_MOD(x0[i], x0[i], Ninvq0, q0); /*Scaling for the convolution*/
		MUL_MOD(x1[i], x1[i], Ninvq1, q1); /*Scaling for the convolution*/
		/*x_crt[i] = h1*x0[i] % q0q1;*/
                MUL_MOD_LONG(x_crt[i], h1,x0[i],q0q1);
		x_crt[i] = q0q1 - x_crt[i];
		/*x_crt[i] += h0*x1[i];*/
                MUL_MOD_LONG(tmp_long,h0,x1[i],q0q1);
                ADD_MOD_LONG(x_crt[i],x_crt[i],tmp_long,q0q1);
		/*x_crt[i] = x_crt[i] % q0q1;*/ //Will now be the integer convolution
		x_crt[i] = x_crt[i] % q; 
	}
			
	/*Complete Bluestein's trick*/
	x[0] = (FFTSHORT) x_crt[(N>>1)-1];	
	for (i = 0; i < n-1; ++i) {
		MUL_MOD(x[i+1], x_crt[(N>>1)+i], Bluestein_mul_inv_709_10486111[i], q);		
	}
				
}


void _FFT_forward_709_10486111(FFTSHORT *x) {
  FFT_forward_709_10486111(x);
}

void _FFT_backward_709_10486111(FFTSHORT *x) {
  int i;
  FFT_backward_709_10486111(x);
  for (i=0; i<709; ++i)
    MUL_MOD(x[i], x[i], 10471321, 10486111);
}


/* Code to compute a Number Theoretic Transform for multiplication in the ring 
F_q[x] / <x^n+1>.
	n = 2048, q = 475404289	*/


/*
We use Gentleman-Sande, decimation-in-frequency FFT, for the forward FFT.
Note that we will not perform the usual scambling / bit-reversal procedure here because we will invert 
the fourier transform using decimation-in-time.
*/
void FFT_forward_2048_475404289(FFTSHORT x[2048]) {
        const FFTSHORT n = 2048;
	const FFTSHORT q = 475404289;
	FFTSHORT index, step;
	FFTSHORT i,j,m;
	FFTSHORT t0,t1;

	step = 1;
	for (m = n>>1; m >= 1; m=m>>1) {
		index = 0;
		for (j = 0 ; j < m; ++j) {
			for (i = j; i < n; i += (m<<1)) {
				ADD_MOD(t0, x[i], x[i+m], q);
				ADD(t1, x[i], q - x[i+m]);
				MUL_MOD(x[i+m], t1, W_2048_475404289[index], q);				
				x[i] = t0;				
			}
			SUB_MODn(index, index, step, n);
		}
		step = step << 1;
	}	 
}

/*
We use Cooley-Tukey, decimation-in-time FFT, for the inverse FFT.
Note that we will not perform the usual scambling / bit-reversal procedure here because we will the forward
fourier transform is using decimation-in-frequency.
*/
void FFT_backward_2048_475404289(FFTSHORT x[2048]) {
        const FFTSHORT n = 2048;
	const FFTSHORT q = 475404289;
	FFTSHORT index, step;
	FFTSHORT i,j,m;
	FFTSHORT t0,t1;

	step = n>>1;
	for (m = 1; m < n; m=m<<1) {
		index = 0;
		for (j = 0 ; j < m; ++j) {
			for (i = j; i < n; i += (m<<1)) {							
				t0 = x[i];
				t0 -= (t0 >= q) ? q : 0;
				MUL_MOD(t1, x[i+m], W_rev_2048_475404289[index], q);				
				ADD(x[i], t0, t1);
				ADD(x[i+m], t0, q - t1);
				
			}
			SUB_MODn(index, index, step, n);
		}
		step = step >> 1;
	}	
	for (i = 0; i < n; ++i) {
		x[i] -= (x[i] >= q) ? q : 0;
	}
}



/* Code to compute a Number Theoretic Transform for multiplication in the ring 
F_q[x] / <x^n+1>.
	n = 2048, q = 950808577	*/

/*
We use Gentleman-Sande, decimation-in-frequency FFT, for the forward FFT.
Note that we will not perform the usual scambling / bit-reversal procedure here because we will invert 
the fourier transform using decimation-in-time.
*/
void FFT_forward_2048_950808577(FFTSHORT x[2048]) {
        const FFTSHORT n = 2048;
	const FFTSHORT q = 950808577;
	FFTSHORT index, step;
	FFTSHORT i,j,m;
	FFTSHORT t0,t1;

	step = 1;
	for (m = n>>1; m >= 1; m=m>>1) {
		index = 0;
		for (j = 0 ; j < m; ++j) {
			for (i = j; i < n; i += (m<<1)) {
				ADD_MOD(t0, x[i], x[i+m], q);
				ADD(t1, x[i], q - x[i+m]);
				MUL_MOD(x[i+m], t1, W_2048_950808577[index], q);				
				x[i] = t0;				
			}
			SUB_MODn(index, index, step, n);
		}
		step = step << 1;
	}	 
}

/*
We use Cooley-Tukey, decimation-in-time FFT, for the inverse FFT.
Note that we will not perform the usual scambling / bit-reversal procedure here because we will the forward
fourier transform is using decimation-in-frequency.
*/
void FFT_backward_2048_950808577(FFTSHORT x[2048]) {
        const FFTSHORT n = 2048;
	const FFTSHORT q = 950808577;
	FFTSHORT index, step;
	FFTSHORT i,j,m;
	FFTSHORT t0,t1;

	step = n>>1;
	for (m = 1; m < n; m=m<<1) {
		index = 0;
		for (j = 0 ; j < m; ++j) {
			for (i = j; i < n; i += (m<<1)) {							
				t0 = x[i];
				t0 -= (t0 >= q) ? q : 0;
				MUL_MOD(t1, x[i+m], W_rev_2048_950808577[index], q);				
				ADD(x[i], t0, t1);
				ADD(x[i+m], t0, q - t1);
				
			}
			SUB_MODn(index, index, step, n);
		}
		step = step >> 1;
	}	
	for (i = 0; i < n; ++i) {
		x[i] -= (x[i] >= q) ? q : 0;
	}
}

