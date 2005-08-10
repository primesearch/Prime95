/**************************************************************
 *
 *  gwdbldbl.cpp
 *
 *  This file contains all the gwnum initialization routines that require
 *  extended precision floating point.  We want to initialize our sin/cos
 *  and FFT weight arrays with doubles that are as accurate as possible.
 *
 *  This is the only C++ routine in the gwnum library.  Since gwnum is
 *  a C based library, we declare all routines here as extern "C".
 * 
 *  Copyright 2005 Just For Fun Software, Inc.
 *  All Rights Reserved.
 *
 **************************************************************/

/* Include files */

#include "gwdbldbl.h"

/* Pick which doubledouble package we will use. */

#define QD
//#define KEITH_BRIGGS

/* Turn on #define that will disable extended precision floating point */
/* registers, as they wreak havoc with double-double library routines. */

#ifndef X86_64
#define x86
#endif

/* Use Hida, Li & Bailey's QD doubledouble C++ package. */

#ifdef QD
#include "dd.cc"
#endif

/* Use Keith Briggs' doubledouble C++ package.  I find the QD package */
/* a better choice. */

#ifdef KEITH_BRIGGS
#define DD_INLINE
#include "doubledouble.cc"
#include "math.cc"
#define	dd_real doubledouble
#define	_2pi TwoPi
#define	_log2 Log2
#endif


/* Now write all the routines that use the dd_real package. */


/* Utility routine to compute many of the constants used by the */
/* assembly language code */

extern "C"
void gwasm_constants (
	double	k,
	signed long c,
	int	sse2,
	int	zero_pad,
	unsigned long fftlen,
	unsigned long bits_in_small_word,
	double	*asm_values)
{
	dd_real arg, sine1, cosine1, sine2, cosine2, sine3, cosine3;
	double	small_word, big_word, temp;
#define MINUS_CARG			asm_values[0]
#define BIGVAL				((float *) &asm_values[1])[0]
#define BIGBIGVAL			((float *) &asm_values[1])[1]
#define XMM_BIGVAL			asm_values[2]
#define XMM_BIGBIGVAL			asm_values[3]
#define SQRTHALF			asm_values[4]
#define ttmp_ff_inv			asm_values[5]
#define XMM_NORM012_FF			asm_values[6]
#define XMM_K_HI			asm_values[7]
#define XMM_K_LO			asm_values[8]
#define XMM_K_HI_2			asm_values[9]
#define XMM_K_HI_1			asm_values[10]
#define XMM_LIMIT_INVERSE		(asm_values+11)
#define XMM_LIMIT_BIGMAX		(asm_values+19)
#define XMM_LIMIT_BIGMAX_NEG		(asm_values+27)
#define P951				asm_values[35]
#define P618				asm_values[36]
#define P309				asm_values[37]
#define M262				asm_values[38]
#define P588				asm_values[39]
#define M162				asm_values[40]
#define M809				asm_values[41]
#define M382				asm_values[42]
#define P866				asm_values[43]
#define P433				asm_values[44]
#define P577				asm_values[45]
#define P975				asm_values[46]
#define P445				asm_values[47]
#define P180				asm_values[48]
#define P623				asm_values[49]
#define M358				asm_values[50]
#define P404				asm_values[51]
#define M223				asm_values[52]
#define M901				asm_values[53]
#define M691				asm_values[54]

/* Do some initial setup */
	
	x86_FIX
	small_word = (double) (1 << bits_in_small_word);
	big_word = (double) (1 << (bits_in_small_word + 1));
	
/* Negate c and store as a double */

	MINUS_CARG = (double) -c;

/* Compute the x87 (64-bit) rounding constants */

	BIGVAL = (float) (3.0 * pow (2.0, 62.0));
	BIGBIGVAL = (float) (big_word * BIGVAL);

/* Compute the SSE2 (53-bit) rounding constants */

	XMM_BIGVAL = 3.0 * pow (2.0, 51.0);
	XMM_BIGBIGVAL = big_word * XMM_BIGVAL;

/* Compute square root of 0.5 */

	SQRTHALF = sqrt (0.5);

/* The sumout value is FFTLEN/2 times larger than it should be.  Create an */
/* inverse to properly calculate the sumout when a multiplication ends. */

	ttmp_ff_inv = 2.0 / (double) fftlen;

/* Compute constant that converts fft_weight_over_fftlen found in the */
/* two-to-minus-phi tables into the true fft_weight value.  This is usually */
/* FFTLEN / 2, but when doing a non-zero-padded FFT this is FFTLEN / 2k. */

	XMM_NORM012_FF = (double) fftlen / 2.0;
	if (!zero_pad) XMM_NORM012_FF /= k;

/* Split k for zero-padded FFTs emulating modulo k*2^n+c */

	XMM_K_HI = floor (k / big_word) * big_word;
	XMM_K_LO = k - XMM_K_HI;
	XMM_K_HI_2 = floor (k / big_word / big_word) * big_word * big_word;
	XMM_K_HI_1 = XMM_K_HI - XMM_K_HI_2;

/* Compute the normalization constants indexed by biglit array entries */

	temp = 1.0 / small_word;	/* Compute lower limit inverse */
	XMM_LIMIT_INVERSE[0] =
	XMM_LIMIT_INVERSE[1] =
	XMM_LIMIT_INVERSE[3] =
	XMM_LIMIT_INVERSE[4] = temp;

					/* Compute lower limit bigmax */
	if (sse2) temp = small_word * XMM_BIGVAL - XMM_BIGVAL;
	else temp = small_word * BIGVAL - BIGVAL;
	XMM_LIMIT_BIGMAX[0] =
	XMM_LIMIT_BIGMAX[1] =
	XMM_LIMIT_BIGMAX[3] =
	XMM_LIMIT_BIGMAX[4] = temp;

	temp = -temp;			/* Negative lower limit bigmax */
	XMM_LIMIT_BIGMAX_NEG[0] =
	XMM_LIMIT_BIGMAX_NEG[1] =
	XMM_LIMIT_BIGMAX_NEG[3] =
	XMM_LIMIT_BIGMAX_NEG[4] = temp;

	temp = 1.0 / big_word;		/* Compute upper limit inverse */
	XMM_LIMIT_INVERSE[2] =
	XMM_LIMIT_INVERSE[5] =
	XMM_LIMIT_INVERSE[6] =
	XMM_LIMIT_INVERSE[7] = temp;

					/* Compute upper limit bigmax */
	if (sse2) temp = big_word * XMM_BIGVAL - XMM_BIGVAL;
	else temp = big_word * BIGVAL - BIGVAL;
	XMM_LIMIT_BIGMAX[2] =
	XMM_LIMIT_BIGMAX[5] =
	XMM_LIMIT_BIGMAX[6] =
	XMM_LIMIT_BIGMAX[7] = temp;

	temp = -temp;			/* Negative upper limit bigmax */
	XMM_LIMIT_BIGMAX_NEG[2] =
	XMM_LIMIT_BIGMAX_NEG[5] =
	XMM_LIMIT_BIGMAX_NEG[6] =
	XMM_LIMIT_BIGMAX_NEG[7] = temp;

/* Initialize the five_reals sine-cosine data. */
/* NOTE: When computing cosine / sine, divide by the 64-bit sine not the */
/* extra precision sine since macros will multiply by the 64-bit sine. */

	arg = dd_real::_2pi / 5.0;		// 2*PI * 1 / 5
	sincos (arg, sine1, cosine1);

	arg = arg * 2.0;			// 2*PI * 2 / 5
	sincos (arg, sine2, cosine2);

	P951 = double (sine1);
	P618 = double (sine2 / P951);		// 0.588 / 0.951

	P309 = double (cosine1);
	M262 = double (cosine2 / P309);		// -0.809 / 0.309

	P588 = double (sine2);
	M162 = double (-sine1 / P588);		// -0.951 / 0.588

	M809 = double (cosine2);
	M382 = double (cosine1 / M809);		// 0.309 / -0.809

/* Initialize the six_reals sine-cosine data. */

	arg = dd_real::_2pi / 3.0;		// 2*PI / 3
	sine1 = sin (arg);			// Compute sine (0.866)

	P866 = double (sine1);
	P433 = double (sine1 * 0.5);		// 0.5 * P866
	P577 = double (dd_real (0.5) / sine1);	// 0.5 / 0.866

/* Initialize the seven_reals sine-cosine data. */

	arg = dd_real::_2pi / 7.0;		// 2*PI * 1 / 7
	sincos (arg, sine1, cosine1);		// cosine (0.623), sine (0.782)

	arg = arg * 2.0;			// 2*PI * 2 / 7
	sincos (arg, sine2, cosine2);		// cosine (-.223), sine (0.975)

	arg = arg * 1.5;			// 2*PI * 3 / 7
	sincos (arg, sine3, cosine3);		// cosine (-.901), sine (0.434)
		
	P975 = double (sine2);
	P445 = double (sine3 / P975);		// 0.434 / 0.975
	P180 = double (sine1 / P975 / P445);	// 0.782 / (0.975 * 0.445)

	P623 = double (cosine1);
	M358 = double (cosine2 / P623);		// -0.223 / 0.623
	P404 = double (cosine3 / P623 / M358);	// -.901 / (.623 * -.358)

	M223 = double (cosine2);
	M901 = double (cosine3);
	M691 = double (cosine1 / M901);		// 0.623 / -0.901

	END_x86_FIX
}

// Utility routine to compute a sin/cos premultiplier or a set of 3
// sine-cosine values.
// This is used during setup, formerly written in assembly language to take
// advantage of the extra precision in the FPU's 80-bit registers. 
// NOTE: When computing cosine / sine, divide by the 64-bit sine
// not the 80-bit sine since macros will multiply by the 64-bit sine.

extern "C"
void gwsincos (
	unsigned long x,
	unsigned long N,
	double	*results)
{
	dd_real arg, sine, cosine;

	x86_FIX
	arg = dd_real::_2pi * (double) x / (double) N;
	sincos (arg, sine, cosine);
	sine += 1E-200;			/* Protect against divide by zero */
	results[0] = sine;
	results[1] = cosine / results[0];
	END_x86_FIX
}

extern "C"
void gwsincos3 (
	unsigned long x,
	unsigned long N,
	double	*results)
{
	dd_real arg1, arg2, arg3, sine, cosine;

	x86_FIX
	arg1 = dd_real::_2pi * (double) x / (double) N;
	sincos (arg1, sine, cosine);
	sine += 1E-200;			/* Protect against divide by zero */
	results[0] = sine;
	results[1] = cosine / results[0];
	arg2 = arg1 + arg1;
	sincos (arg2, sine, cosine);
	sine += 1E-200;			/* Protect against divide by zero */
	results[2] = sine;
	results[3] = cosine / results[2];
	arg3 = arg2 + arg1;
	sincos (arg3, sine, cosine);
	sine += 1E-200;			/* Protect against divide by zero */
	results[4] = sine;
	results[5] = cosine / results[4];
	END_x86_FIX
}


//
// Utility routines to compute fft weights
//
// The FFT weight for the j-th FFT word doing a 2^q+c weighted transform is
//	2 ^ (ceil (j*q/FFTLEN) - j*q/FFTLEN)   *    abs(c) ^ j/FFTLEN
//

static	dd_real gw__bits_per_word;
static	int gw__c_is_one;
static	dd_real gw__log2_abs_c_div_fftlen;
static	dd_real gw__fftlen_inverse;
static	dd_real gw__over_fftlen;
static	double gwdbl__bits_per_word;
static	double gwdbl__log2_abs_c_div_fftlen;

extern "C"
void gwfft_weight_setup (
	int	zero_pad,
	double	k,
	unsigned long n,
	signed long c,
	unsigned long fftlen)
{
	x86_FIX
	gw__fftlen_inverse = dd_real (1.0) / dd_real ((double) fftlen);
	if (zero_pad) {
		gw__bits_per_word =
			dd_real ((double) (n + n)) * gw__fftlen_inverse;
		gw__c_is_one = 1;
		gw__over_fftlen = dd_real (2.0) * gw__fftlen_inverse;
	} else {
		gw__bits_per_word =
			(dd_real ((double) n) +
			 log (dd_real (k)) / dd_real::_log2) *
			gw__fftlen_inverse;
		gw__c_is_one = (abs ((int) c) == 1);
		gw__log2_abs_c_div_fftlen =
			log (dd_real (abs ((int) c))) / dd_real::_log2 *
			gw__fftlen_inverse;
		gw__over_fftlen = dd_real (k * 2.0) * gw__fftlen_inverse;
	}
	gwdbl__bits_per_word = (double) gw__bits_per_word;
	gwdbl__log2_abs_c_div_fftlen = (double) gw__log2_abs_c_div_fftlen;
	END_x86_FIX
}

extern "C"
double gwfft_weight (
	unsigned long j)
{
	dd_real temp, twopow, result;

	x86_FIX
	temp = dd_real ((double) j) * gw__bits_per_word;
	twopow = ceil (temp) - temp;
	if (! gw__c_is_one)
		twopow += gw__log2_abs_c_div_fftlen * dd_real ((double) j);
	result = exp (dd_real::_log2 * twopow);
	END_x86_FIX
	return (double (result));
}

// Like the above, but faster and does not guarantee quite as much accuracy.

extern "C"
double gwfft_weight_sloppy (
	unsigned long j)
{
	dd_real temp, twopow;

	x86_FIX
	temp = dd_real ((double) j) * gw__bits_per_word;
	twopow = ceil (temp) - temp;
	if (! gw__c_is_one)
		twopow += gw__log2_abs_c_div_fftlen * dd_real ((double) j);
	END_x86_FIX
	return (pow (2.0, double (twopow)));
}

// Compute the inverse of the fft weight

extern "C"
double gwfft_weight_inverse (
	unsigned long j)
{
	dd_real temp, twopow, result;

	x86_FIX
	temp = dd_real ((double) j) * gw__bits_per_word;
	twopow = ceil (temp) - temp;
	if (! gw__c_is_one)
		twopow += gw__log2_abs_c_div_fftlen * dd_real ((double) j);
	result = exp (dd_real::_log2 * -twopow);
	END_x86_FIX
	return (double (result));
}

// Like the above, but faster and does not guarantee quite as much accuracy.

extern "C"
double gwfft_weight_inverse_sloppy (
	unsigned long j)
{
	double	tempdbl, twopowdbl;

	tempdbl = (double) j * gwdbl__bits_per_word;
	twopowdbl = ceil (tempdbl) - tempdbl;
	if (twopowdbl < 0.001 || twopowdbl > 0.999) {
		dd_real temp;

		x86_FIX
		temp = dd_real ((double) j) * gw__bits_per_word;
		twopowdbl = (double) (ceil (temp) - temp);
		END_x86_FIX
	}

	if (! gw__c_is_one)
		twopowdbl += gwdbl__log2_abs_c_div_fftlen * (double) j;
	return (pow (2.0, - double (twopowdbl)));
}

// This computes the inverse FFT weight multiplied by the appropriate constant
// to produce an integer during an FFT multiply's normalize stage.  This
// constant is 2/FFTLEN for a zero-padded FFT and k*2/FFTLEN for a
// non-zero-padded FFT.

extern "C"
double gwfft_weight_inverse_over_fftlen (
	unsigned long j)
{
	dd_real temp, twopow, result;

	x86_FIX
	temp = dd_real ((double) j) * gw__bits_per_word;
	twopow = ceil (temp) - temp;
	if (! gw__c_is_one)
		twopow += gw__log2_abs_c_div_fftlen * dd_real ((double) j);
	result = exp (dd_real::_log2 * -twopow) * gw__over_fftlen;
	END_x86_FIX
	return (double (result));
}

// This computes the three FFT weights in one call.  It is faster than
// calling the above individually.

extern "C"
void gwfft_weights3 (
	unsigned long j,
	double	*fft_weight,
	double	*fft_weight_inverse,
	double	*fft_weight_inverse_over_fftlen)
{
	dd_real temp, twopow, weight;

	x86_FIX
	temp = dd_real ((double) j) * gw__bits_per_word;
	twopow = ceil (temp) - temp;
	if (! gw__c_is_one)
		twopow += gw__log2_abs_c_div_fftlen * dd_real ((double) j);
	weight = exp (dd_real::_log2 * twopow);
	*fft_weight = double (weight);
	weight = dd_real (1.0) / weight;
	if (fft_weight_inverse != NULL)
		*fft_weight_inverse = double (weight);
	if (fft_weight_inverse_over_fftlen != NULL)
		*fft_weight_inverse_over_fftlen = double (weight * gw__over_fftlen);
	END_x86_FIX
}

// Returns log2(fft_weight).  This is used in determining the FFT weight
// fudge factor in two-pass FFTs.  This is much faster than computing the
// fft_weight because it eliminates a call to the double-double exp routine.

extern "C"
double gwfft_weight_exponent (
	unsigned long j)
{
	double	tempdbl, twopowdbl;
	dd_real temp, twopow;

// For speed, try this with plain old doubles first

	if (j == 0) return (0);
	tempdbl = (double) j * gwdbl__bits_per_word;
	twopowdbl = ceil (tempdbl) - tempdbl;
	if (twopowdbl > 0.001 && twopowdbl < 0.999) return (twopowdbl);

// If at all uncertain of the result, use doubledoubles to do the calculation

	x86_FIX
	temp = dd_real ((double) j) * gw__bits_per_word;
	twopow = ceil (temp) - temp;
	END_x86_FIX
	return (double (twopow));
}

//
// Utility routine to compute fft base for j-th fft word
//
// The FFT base for the j-th FFT word doing a 2^q+c weighted transform is
//	ceil (j*q/FFTLEN)
// This routine returns ceil (j*q/FFTLEN) taking great care to return a
// value accurate to 53 bits.  This is important when j*q is really close to
// a multiple of FFTLEN (admittedly quite rare).  It would be really bad if
// rounding differences caused this routine to compute ceil (j*q/FFTLEN)
// differently than the weighting functions.
//

extern "C"
unsigned long gwfft_base (
	unsigned long j)
{
	double	tempdbl, ceildbl, diffdbl;
	dd_real temp;
	unsigned long twopow;

// For speed, try this with plain old doubles first

	if (j == 0) return (0);
	tempdbl = (double) j * gwdbl__bits_per_word;
	ceildbl = ceil (tempdbl);
	diffdbl = ceildbl - tempdbl;
	if (diffdbl > 0.001 && diffdbl < 0.999)
		return ((unsigned long) ceildbl);

// If at all uncertain of the result, use doubledoubles to do the calculation

	x86_FIX
	temp = dd_real ((double) j) * gw__bits_per_word;
	twopow = (int) ceil (temp);
	END_x86_FIX
	return (twopow);
}
