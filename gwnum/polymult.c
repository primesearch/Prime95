/*----------------------------------------------------------------------
| polymult.c
|
| This file contains the C routines for fast polynomial multiplication
| 
|  Copyright 2021 Mersenne Research, Inc.  All rights reserved.
+---------------------------------------------------------------------*/

/* Include files */

#include <stdlib.h>
#include <math.h>
#include <memory.h>
#include "cpuid.h"
#include "gwnum.h"
#include "gwutil.h"
#include "polymult.h"
#if defined (SSE2) || defined (AVX) || defined (FMA) || defined (AVX512)
#include <immintrin.h>
#endif

// TO DO:
// 1) Support non-power-of-two FFTs below 32?  Support FFT size requiring a radix-2 step?
// 2) Zero padded gwnums.  Test all FFT sizes and FFT types -- especially funky PFA ffts.
// 7) Variable GW_HEADER_SIZE.
// 9) Is a fast convert-gwnum-to-larger-FFT feasible? necessary?
// 10) Tune for Brute / Karatsuba / FFT crossovers.
// 12) Tune routine that returns needed safety margin needed for given poly sizes
// 14) Compressed gwnums.  Does it pay off?
// 15) Add polymult to gwtest or as a standalone program
// 16) Handle out-of-memory gracefully
// 17) Right collection of PFA FFT sizes?
// 19) For large ffts, use Pavel's idea of copying the complex values to a gwnum and use an unweighted FFT with a special no-carry-propagation norm routine.
//     Perhaps only do this for most common architectures, such as AVX-512 and FMA 64-bit.
// 21) Should polymult be its own library or a part of the gwnum library?
// 22) Can we add asserts that polymult isn't violating EXTRA_BITS?
// 25) Optimize.  Read/write several cache lines at a time (fill up L2 cache?).  Prefetch future cache lines.
// 28) Optimize brute_force to do less compute when outvec entries are NULL.
// 28) Optimize brute_force to do less compute when one or both invec entries are RLPs.
// 29) monic polynomials are tricky if FFT(1) is not all ones.  Use GW_FFT1 and treat as non-monic?
// 30) option to negate outputs??
// 32) If monic input and polymult is doing an FFT and there is zero padding going on, then pad with one (or GW_FFT1) rather than use monic_add.  Be careful
//	that this does not cause any unwanted wrapping of data!
// 33) Are there cases with brute/Karatsuba would benefit from using a bigger poly rather than using monic_add?
// 34) Support GWMUL_BY_CONST, GWMUL_ADD_CONST option in gwunfft2.
// 35) A monic 1024 times monic 1024 using length 2048 FFT could apply the two ones in the FFT vectors and patch outvec[0] (possibly with gwsetaddin(-1)).
//	A monic 1023 times monic 1023 using length 2048 FFT could apply the two ones in the FFT vectors with no patching rather than monic  post-processing.
//	That should be a little faster.  Other monic combinations have similar optimizations.
// 36) For cached FFTs, is there a way to avoid the transpose(write_line/read_line) step?
// 37) If transpose is expensive, can we offer a way to "pre-transpose" invec1 for use multiple times?  With compression of exponent field in double??
//	Compression still needs semi-random/sequential access to invec1 for multi-threading.
// 38)  For poly #2, is there a quick block compression we can do where the extra read/write cost will be more than compensated by the few more extra
//	gwnums we can add to poly #2?
// 39)  Think more on unbalanced poly mults.  Is a 100 x 200 mult better done as two 100 x 100 multiplies?  Answer: Never for brute force, maybe for FFT.
//	If we can save the smaller RLP in FFTed state using smaller FFT size that might influence our decision.
// 40) Is invec1_negate code OK when GW_FFT1 is not all ones?  I would think so.
// 41) Add callback routine to check for ^C/ESC
// 42) Transparent / or easy way to store gwnums on disk as part of input and/or output vectors.  Picture using 12GB memory and 50GB disk for really huge polys.
// 43) Read_line/write_line that reads/writes data with stride=1?
// 44) Cache several sincos tables for different FFT sizes.
// 45) If twiddle_stride gets large, build a smaller twiddle table for smaller strides.
// 46) Look for more FMA opportunities.
// 48) Test if hyperthreading is beneficial
// 49) Problems if FFT(1)!=1!!! The pointwise mul code ought to divide by FFT(1) OR have a special asm entry point for unfft that deals with data that
//	has already had a pointwise mul.  Furthermore, the monic_adjustments code can't simply add data - it must add FFT(1)*data.
// 51) Scratch area (like prime95) to avoid large strides,  Read_partial_line could populate the scratch area (makes monic_adjustment difficult)?
// 52) Half-sized DCT for RLPs?  The DJB FFTs may make reordering DCT outputs for an inverse FFT tricky.
// 53) Can P-1's circular monic-poly times monic-RLP be safely done without monic_adjustments allowing us to save the FFT'ed RLP that is half-sized due to DCT?
// 54) Does monolithic two-forward-FFTs-pointmul-and-one-inverse-FFT routine suggest different preferred pass sizes (use half the L2 cache)?
// 55) For VFMA machines, store twiddles as cos/sin, sin to allow a few more FMA opportunities
// 56) Move mul by 1/FFTLEN to last unfft level (perhaps into the twiddle table itself).  Saves 65%+ of the multiplies if VFMA is available.
// 57) Option to save memory by processing AVX-512 data using FMA3.  (reduces sizeof(CVDT)).

// Forward declarations required for the SSE2, AVX, FMA, AVX512 compiles
int next_line (pmhandle *pmdata);
unsigned long cache_line_offset (gwhandle *gwdata, unsigned long i);
void generate_sincos (pmhandle *pmdata, const unsigned int size);
void polymult_line_sse2 (pmhandle *pmdata, gwnum *gw_invec1, int invec1_size, gwnum *gw_invec2, int invec2_size, gwnum *gw_outvec, int outvec_size, int options);
void polymult_line_avx (pmhandle *pmdata, gwnum *gw_invec1, int invec1_size, gwnum *gw_invec2, int invec2_size, gwnum *gw_outvec, int outvec_size, int options);
void polymult_line_fma (pmhandle *pmdata, gwnum *gw_invec1, int invec1_size, gwnum *gw_invec2, int invec2_size, gwnum *gw_outvec, int outvec_size, int options);
void polymult_line_avx512 (pmhandle *pmdata, gwnum *gw_invec1, int invec1_size, gwnum *gw_invec2, int invec2_size, gwnum *gw_outvec, int outvec_size, int options);

#define INVEC1_OPTIONS	(POLYMULT_INVEC1_MONIC | POLYMULT_INVEC1_RLP | POLYMULT_INVEC1_NEGATE | POLYMULT_INVEC1_FFTED)
#define INVEC2_OPTIONS	(POLYMULT_INVEC2_MONIC | POLYMULT_INVEC2_RLP | POLYMULT_INVEC2_NEGATE | POLYMULT_INVEC2_FFTED)

#if defined (AVX512)
#define VLEN		512			// Vector is 256 bits wide
#define VDT		__m512d			// Vector data type
#define VFMA		1			// FMA is supported
#define	read_line	read_line_avx512	// Naming extension for this option
#define	write_line	write_line_avx512
#define	brute_line	brute_line_avx512
#define	karatsuba_line	karatsuba_line_avx512
#define	fft_line	fft_line_avx512
#define pick_pass_sizes	pick_pass_sizes_avx512
#define	monic_line_add_hi monic_line_add_hi_avx512
#define	monic_line_add_lo monic_line_add_lo_avx512
#define	monic_line_adjustments monic_line_adjustments_avx512
#define	polymult_line	polymult_line_avx512
#define broadcastsd	_mm512_set1_pd
#define addpd		_mm512_add_pd
#define subpd		_mm512_sub_pd
#define mulpd		_mm512_mul_pd
#define fmaddpd		_mm512_fmadd_pd
#define fmsubpd		_mm512_fmsub_pd
#define fnmaddpd	_mm512_fnmadd_pd
#elif defined (FMA)
#define VLEN		256			// Vector is 256 bits wide
#define VDT		__m256d			// Vector data type
#define VFMA		1			// FMA is supported
#define	read_line	read_line_fma		// Naming extension for this option
#define	write_line	write_line_fma
#define	brute_line	brute_line_fma
#define	karatsuba_line	karatsuba_line_fma
#define	fft_line	fft_line_fma
#define pick_pass_sizes	pick_pass_sizes_fma
#define	monic_line_add_hi monic_line_add_hi_fma
#define	monic_line_add_lo monic_line_add_lo_fma
#define	monic_line_adjustments monic_line_adjustments_fma
#define	polymult_line	polymult_line_fma
#define broadcastsd	_mm256_set1_pd
#define addpd		_mm256_add_pd
#define subpd		_mm256_sub_pd
#define mulpd		_mm256_mul_pd
#define fmaddpd		_mm256_fmadd_pd
#define fmsubpd		_mm256_fmsub_pd
#define fnmaddpd	_mm256_fnmadd_pd
#elif defined (AVX)
#define VLEN		256			// Vector is 256 bits wide
#define VDT		__m256d			// Vector data type
#define VFMA		0			// FMA not supported
#define	read_line	read_line_avx		// Naming extension for this option
#define	write_line	write_line_avx
#define	brute_line	brute_line_avx
#define	karatsuba_line	karatsuba_line_avx
#define	fft_line	fft_line_avx
#define pick_pass_sizes	pick_pass_sizes_avx
#define	monic_line_add_hi monic_line_add_hi_avx
#define	monic_line_add_lo monic_line_add_lo_avx
#define	monic_line_adjustments monic_line_adjustments_avx
#define	polymult_line	polymult_line_avx
#define broadcastsd	_mm256_set1_pd
#define addpd		_mm256_add_pd
#define subpd		_mm256_sub_pd
#define mulpd		_mm256_mul_pd
#elif defined (SSE2)
#define VLEN		128			// Vector is 256 bits wide
#define VDT		__m128d			// Vector data type
#define VFMA		0			// FMA not supported
#define	read_line	read_line_sse2		// Naming extension for this option
#define	write_line	write_line_sse2
#define	brute_line	brute_line_sse2
#define	karatsuba_line	karatsuba_line_sse2
#define	fft_line	fft_line_sse2
#define pick_pass_sizes	pick_pass_sizes_sse2
#define	monic_line_add_hi monic_line_add_hi_sse2
#define	monic_line_add_lo monic_line_add_lo_sse2
#define	monic_line_adjustments monic_line_adjustments_sse2
#define	polymult_line	polymult_line_sse2
#define broadcastsd	_mm_set1_pd
#define addpd		_mm_add_pd
#define subpd		_mm_sub_pd
#define mulpd		_mm_mul_pd
#else
#define VLEN		64			// Vector is 256 bits wide
#define VDT		double			// Vector data type
#define VFMA		0			// FMA not supported
#define	read_line	read_line_dbl		// Naming extension for this option
#define	write_line	write_line_dbl
#define	brute_line	brute_line_dbl
#define	karatsuba_line	karatsuba_line_dbl
#define	fft_line	fft_line_dbl
#define pick_pass_sizes	pick_pass_sizes_dbl
#define	monic_line_add_hi monic_line_add_hi_dbl
#define	monic_line_add_lo monic_line_add_lo_dbl
#define	monic_line_adjustments monic_line_adjustments_dbl
#define	polymult_line	polymult_line_dbl
#define broadcastsd	_mm64_set1_pd
#define addpd		_mm64_add_pd
#define subpd		_mm64_sub_pd
#define mulpd		_mm64_mul_pd
// Emulate vector length 1 using Intel-style instrinsic naming conventions
#define __m64d			double
#define _mm64_set1_pd(a)	(a)
#define _mm64_add_pd(a,b)	((a)+(b))
#define _mm64_sub_pd(a,b)	((a)-(b))
#define _mm64_mul_pd(a,b)	((a)*(b))
#endif

// Emulate FMA on CPUs that do not have an FMA instruction
#if !VFMA
#define fmaddpd(a,b,c)	addpd(mulpd(a,b),c)
#define fmsubpd(a,b,c)	subpd(mulpd(a,b),c)
#define fnmaddpd(a,b,c)	subpd(c,mulpd(a,b))
#endif

// Structure for complex VDTs
typedef struct {
	VDT	real;
	VDT	imag;
} CVDT;
#define cvload(a,b,c)		(a).real = broadcastsd (b), (a).imag = broadcastsd (c)
#define cvadd(a,b,c)		(a).real = addpd ((b).real, (c).real), (a).imag = addpd ((b).imag, (c).imag)				// a = b + c
#define cvsub(a,b,c)		(a).real = subpd ((b).real, (c).real), (a).imag = subpd ((b).imag, (c).imag)				// a = b - c
#define cviadd(a,b,c)		(a).real = subpd ((b).real, (c).imag), (a).imag = addpd ((b).imag, (c).real)				// a = b + i*c
#define cvisub(a,b,c)		(a).real = addpd ((b).real, (c).imag), (a).imag = subpd ((b).imag, (c).real)				// a = b - i*c
#define cvaddfm(a,b,const_d,c)	(a).real = fmaddpd (const_d, (c).real, (b).real), (a).imag = fmaddpd (const_d, (c).imag, (b).imag)	// a = b + const_d*c
#define cvsubfm(a,b,const_d,c)	(a).real = fnmaddpd (const_d, (c).real, (b).real), (a).imag = fnmaddpd (const_d, (c).imag, (b).imag)	// a = b - const_d*c
#define cvmul(a,b,c)		(a).real = fmsubpd ((b).real, (c).real, mulpd ((b).imag, (c).imag)), \
				(a).imag = fmaddpd ((b).imag, (c).real, mulpd ((b).real, (c).imag))					// a = b * c
#define cvconjmul(a,b,c)	(a).real = fmaddpd ((b).real, (c).real, mulpd ((b).imag, (c).imag)), \
				(a).imag = fmsubpd ((b).imag, (c).real, mulpd ((b).real, (c).imag))					// a = b * conjugate(c)

#if !defined (SSE2) && !defined (AVX) && !defined (FMA) && !defined (AVX512)

// Get the needed safety_margin required for an invec1_size by invec2_size polymult
double polymult_safety_margin (int invec1_size, int invec2_size) {
	int	n = (invec1_size < invec2_size) ? invec1_size : invec2_size;
	// The number of EXTRA_BITS (output FFT bits) consumed by a polymult ~= log2(n)/2.  However, gwinit's safety margin parameter refers to
	// input FFT bits.  Thus, divide by another 2.
	return (log2 ((double) n) / 4.0);
}

// Get the FFT size that will be used for an n = invec1_size + invec2_size polymult
int polymult_fft_size (int n)
{
//GW: These FFT sizes haven't been timed.  Some may be slower than a larger FFT or more radix-3/5 options might be profitable
	int	pfas[] = {40, 45, 48, 50, 54, 60, 64, 72, 75};		// Provide a variety of radix-3 and radix-5 FFT sizes

	// Power of two when n is small -- brute force is probably faster anyway
	if (n <= 32) for (int i = 4; ; i *= 2) if (i >= n) return (i);

	// Increase n by a power of two until it
	for (int two_multiplier = 1; ; two_multiplier *= 2) {
		if (n > two_multiplier * 75) continue;
		for (int i = 0; i < sizeof (pfas) / sizeof (int); i++) {
			int fft_size = two_multiplier * pfas[i];
			if (fft_size >= n && (fft_size % 4 == 0 || fft_size % 2 == 1)) return (fft_size);  // Weed out fft sizes that would require a radix-2 step
		}
	}
}

// Initialize a polymult handle
void polymult_init (
	pmhandle *pmdata,		// Handle for polymult library that needs initializing
	gwhandle *gwdata)		// Handle for gwnum FFT library
{
	memset (pmdata, 0, sizeof (pmhandle));
	pmdata->gwdata = gwdata;

	// Initialize default number of threads
	pmdata->num_threads = gwget_num_threads (gwdata);

	// Initialize default breakeven points
	if (pmdata->gwdata->cpu_flags & CPU_AVX) {
		pmdata->KARAT_BREAK = 32;			//GW: Fix me
		pmdata->FFT_BREAK = 64;				//GW: Fix me
	}
	else {
		pmdata->KARAT_BREAK = 32;			//GW: Fix me
		pmdata->FFT_BREAK = 64;				//GW: Fix me
	}

	// Set default cache size to optimize FFTs for
	pmdata->L2_CACHE_SIZE = 256;
}

// Terminate use of a polymult handle.  Free up memory.
void polymult_done (
	pmhandle *pmdata)		// Handle for polymult library
{
	if (pmdata->thread_ids != NULL) {
		pmdata->termination_in_progress = TRUE;
		gwevent_signal (&pmdata->work_to_do);
		for (int i = 1; i < pmdata->num_threads; i++) gwthread_wait_for_exit (&pmdata->thread_ids[i]);
		free (pmdata->thread_ids);
		pmdata->thread_ids = NULL;
		gwmutex_destroy (&pmdata->poly_mutex);
		gwevent_destroy (&pmdata->work_to_do);
		gwevent_destroy (&pmdata->helpers_done);
		gwevent_destroy (&pmdata->main_can_wakeup);
	}
	aligned_free (pmdata->twiddles1), pmdata->twiddles1 = NULL;
	aligned_free (pmdata->twiddles2), pmdata->twiddles2 = NULL;
}


// Return the next line that needs processing
int next_line (
	pmhandle *pmdata)		// Handle for polymult library
{
	int	return_value;

	// Use mutex to guard next_line counter, return next available line
	if (pmdata->num_threads > 1) gwmutex_lock (&pmdata->poly_mutex);
	return_value = pmdata->next_line++;
	if (pmdata->num_threads > 1) gwmutex_unlock (&pmdata->poly_mutex);
	return (return_value);
}

//GW -- move to gwnum.c?

/* Get offset to nth "cache line" from an FFTed gwnum.  Where a cache line is 64 bytes. */

unsigned long cache_line_offset (
	gwhandle *gwdata,	/* Handle initialized by gwsetup */
	unsigned long i)
{
	unsigned long addr, grps;

/* Memory layouts for AVX-512 machines */

	if (gwdata->cpu_flags & CPU_AVX512F) {
		if (gwdata->PASS1_SIZE == 0) {
			addr = i * 64;
			/* Now optionally add 64 pad bytes every 4KB */
			if (gwdata->FOURKBGAPSIZE) addr = addr + (addr / (gwdata->FOURKBGAPSIZE << 6)) * 64;
		}
		else {
			grps = i / (gwdata->PASS2_SIZE / 4);
			addr = i * 64;
			/* Now optionally add bytes for every 4KB page and every pass 2 block */
			addr = addr + (addr >> 12) * gwdata->FOURKBGAPSIZE + grps * gwdata->PASS2GAPSIZE;
		}
	}

/* Memory layouts for AVX machines */

	else if (gwdata->cpu_flags & CPU_AVX) {
		if (gwdata->PASS2_SIZE == 0) {
			addr = i * 64;
			/* Now optionally add 64 pad bytes every 1KB, 2KB or 4KB */
			if (gwdata->FOURKBGAPSIZE) addr = addr + (addr / (gwdata->FOURKBGAPSIZE << 6)) * 64;
		}
		else {
			grps = i / (gwdata->PASS2_SIZE / 4);
			addr = i * 64;
			/* Now add 64 bytes every 4KB and one pass2gapsize for every pass 2 block. */
			addr = addr + (addr >> 12) * gwdata->FOURKBGAPSIZE + grps * gwdata->PASS2GAPSIZE;
		}
	}

/* Memory layouts for SSE2 machines */

	else if (gwdata->cpu_flags & CPU_SSE2) {

/* PFA-style FFTs are a little tricker.  See assembly code for example.	*/
		if (gwdata->PASS2_SIZE == 0) {
			addr = i * 64;
		}
//GW: BUG!?
/* and PFA layouts are even funkier.					*/
		else if (gwdata->FFT_TYPE == FFT_TYPE_HOME_GROWN) {
			grps = i / (gwdata->PASS2_SIZE / 2);
			addr = i * 64;
			/* Now add 128 bytes every 8KB and one pass2gapsize for every pass 2 block. */
			addr = addr + (addr >> 13) * 128 + grps * gwdata->PASS2GAPSIZE;
		}
/* Newer traditional radix-4 large FFTs use don't have a special layout for PFA. */
		else {
			grps = i / (gwdata->PASS2_SIZE / 2);
			addr = i * 64;
			/* Now add 128 bytes every 8KB and one pass2gapsize for every pass 2 block. */
			addr = addr + (addr >> 13) * 128 + grps * gwdata->PASS2GAPSIZE;
		}
	}

/* One pass x87 FFTs use a near flat memory model. */

//BUG --- complex vals are 8 away from reals (at least in integer space) ---- not an issue if polymult is only for x64
	else if (gwdata->PASS2_SIZE == 0) {
		addr = i * 64;
	}
	else {
		grps = (i * 8) / gwdata->PASS2_SIZE;
		addr = i * 64 + grps * 64;
	}

/* Return the offset */

	return (addr);
}


// Helper routines for FFT code

// Generate radix-3/4/5 twiddle factors
void generate_sincos (
	pmhandle *pmdata,		// Handle for polymult library
	const unsigned int fft_size)
{
	const long double twopi = 2.0L * 3.1415926535897932384626433L;
	unsigned int size = fft_size;

	// Return quickly if twiddles have already been initialized
	if (pmdata->twiddles_initialized && pmdata->twiddles_initialized == size) return;

	// Make sure only one thread initializes the twiddles
	if (pmdata->num_threads > 1) gwmutex_lock (&pmdata->poly_mutex);

	// If some other thread beat us to the initialization, then return
	if (pmdata->twiddles_initialized && pmdata->twiddles_initialized == size) {
		if (pmdata->num_threads > 1) gwmutex_unlock (&pmdata->poly_mutex);
		return;
	}

	// Free the previous twiddles
	if (pmdata->twiddles_initialized) {
//GW: accept size < pmdata->twiddles_initialized (at least for power of 2)
//GW: cache several different table sizes?
		aligned_free (pmdata->twiddles1), pmdata->twiddles1 = NULL;
		aligned_free (pmdata->twiddles2), pmdata->twiddles2 = NULL;
	}

	// Generate twiddles for radix-3 butterflies
	if (size % 3 == 0) {
		long double twopi_over_size = twopi / (long double) size;
		pmdata->twiddles1 = (double *) aligned_malloc (size/3 * 2 * sizeof (double), 64);
		for (unsigned int i = 0; i < size/3; i++) {
			long double angle = twopi_over_size * (long double) i;
			pmdata->twiddles1[2*i] = (double) cosl (angle);
			pmdata->twiddles1[2*i+1] = (double) sinl (angle);
		}
		while (size % 3 == 0) size /= 3;
	}

	// Generate twiddles for radix-4/5 butterflies
	long double twopi_over_size = twopi / (long double) size;
	pmdata->twiddles2 = (double *) aligned_malloc (size/4 * 2 * 2 * sizeof (double), 64);
	for (unsigned int i = 0; i < size/4; i++) {
		double long angle = twopi_over_size * (long double) i;
		pmdata->twiddles2[2*(2*i)] = (double) cosl (angle);
		pmdata->twiddles2[2*(2*i)+1] = (double) sinl (angle);
		pmdata->twiddles2[2*(2*i+1)] = (double) cosl (2.0L*angle);
		pmdata->twiddles2[2*(2*i+1)+1] = (double) sinl (2.0L*angle);
	}

	pmdata->twiddles_initialized = fft_size;
	if (pmdata->num_threads > 1) gwmutex_unlock (&pmdata->poly_mutex);
}

#endif

//--------------------------------------------------------------------
//		Read/write cache lines from/to a gwnum
//    Requires intimate knowledge of layout of gwnums in memory
//--------------------------------------------------------------------

#if defined (AVX512)

bool read_line_avx512 (gwhandle *gwdata, gwnum *data, int size, int index, CVDT *vec, int options, bool post_process_monics)
{
	// If not post-processing monics and the input is a monic RLP, then output a 1+0i
	if (!post_process_monics && (options & (POLYMULT_INVEC1_MONIC | POLYMULT_INVEC2_MONIC)) && (options & (POLYMULT_INVEC1_RLP | POLYMULT_INVEC2_RLP))) {
//GW: check if GW_FFT1 should be used instead
		vec[0].real = broadcastsd (1.0);
		vec[0].imag = broadcastsd (0.0);
		vec++;
	}

	// If input is an RLP then read the data into the high half of vec.  We will duplicate the data as the last step.
	if (options & (POLYMULT_INVEC1_RLP | POLYMULT_INVEC2_RLP)) vec += (size-1);

	// Transformed real data ends up with 2 real values and FFTLEN/2-1 complex values.  Transformed complex data has FFTLEN/2 complex values.
	// Adjust index to handle the extra data point that must be returned for transformed real data.
	if (!gwdata->ALL_COMPLEX_FFT) index--;

	// Zero padded FFTs must also return for processing the 7 real values stored in the header
	if (gwdata->ZERO_PADDED_FFT) index -= 7;
//BUG-Zero-padded not handled yet.  Need bigger gwnum header (if we allow storing ffted polys)
//BUG-current use of header location -10 conflicts with zero-padded use of location -10.

	// If we've read all the gwnum data, return FALSE
	if (index * 16 >= (int) gwdata->FFTLEN) return (FALSE);

	// AVX512 implementation
	if (index == -1) {			// First real-only value
		for (int j = 0; j < size; j++) {
			if (data[j] == NULL) { memset (&vec[j], 0, sizeof (CVDT)); continue; }
			vec[j].real = _mm512_set_pd (0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, data[j][0]);
//bug? return 0.0 instead of -11
//			vec[j].imag = _mm512_set_pd (0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, data[j][-11]);
			vec[j].imag = _mm512_set_pd (0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0);
		}
	}
	else if (index == 0 && !gwdata->ALL_COMPLEX_FFT) {			// Second real-only value
		for (int j = 0; j < size; j++) {
			if (data[j] == NULL) { memset (&vec[j], 0, sizeof (CVDT)); continue; }
			vec[j].real = _mm512_set_pd (data[j][7], data[j][6], data[j][5], data[j][4], data[j][3], data[j][2], data[j][1], data[j][8]);
//bug? return 0.0 instead of -10
//			vec[j].imag = _mm512_set_pd (data[j][7], data[j][6], data[j][5], data[j][-10]);
			vec[j].imag = _mm512_set_pd (data[j][15], data[j][14], data[j][13], data[j][12], data[j][11], data[j][10], data[j][9], 0.0);
		}
	}
	else {				// Process complex values where imaginary part is in the high half of the cache line
		unsigned long offset = cache_line_offset (gwdata, index * 2) / sizeof (VDT);
		for (int j = 0; j < size; j++) {
			if (data[j] == NULL) { memset (&vec[j], 0, sizeof(CVDT)); continue; }
			ASSERTG (FFT_state (data[j]) == FULLY_FFTed);
			vec[j].real = ((VDT *) data[j])[offset];
			vec[j].imag = ((VDT *) data[j])[offset+1];
		}
	}

	// Negate coefficients if requested
	if (options & (POLYMULT_INVEC1_NEGATE | POLYMULT_INVEC2_NEGATE)) {
		VDT zero = broadcastsd (0.0);
		for (int j = 0; j < size; j++) {
			vec[j].real = subpd (zero, vec[j].real);
			vec[j].imag = subpd (zero, vec[j].imag);
		}
	}

	// Duplicate RLP data in reverse order
	if (options & (POLYMULT_INVEC1_RLP | POLYMULT_INVEC2_RLP)) {
		for (int j = 1; j < size; j++) vec[-j] = vec[j];
	}

	// If not post-processing monics and the input is monic, then output a 1+0i
	if (!post_process_monics && (options & (POLYMULT_INVEC1_MONIC | POLYMULT_INVEC2_MONIC))) {
//GW: check if GW_FFT1 should be used instead
		vec[size].real = broadcastsd (1.0);
		vec[size].imag = broadcastsd (0.0);
	}

	// Return success
	return (TRUE);
}

void write_line_avx512 (gwhandle *gwdata, gwnum *data, int size, int index, CVDT *vec, int options)
{
	// If either input was a monic RLP, back up the output vector pointer to the newly created entry by monic_line_adjustments
	if (((options & POLYMULT_INVEC1_MONIC) && (options & POLYMULT_INVEC1_RLP)) || ((options & POLYMULT_INVEC2_MONIC) && (options & POLYMULT_INVEC2_RLP))) vec--;

	// If both inputs are RLPs, then the output is an RLP too.  Skip writing half of the data.
	if ((options & POLYMULT_INVEC1_RLP) && (options & POLYMULT_INVEC2_RLP)) vec += (size-1);

	// Transformed real data ends up with 2 real values and FFTLEN/2-1 complex values.  Transformed complex data has FFTLEN/2 complex values.
	// Adjust index to handle the extra data point that must be returned for transformed real data.
	if (!gwdata->ALL_COMPLEX_FFT) index--;

	// Zero padded FFTs must also return for processing the 7 real values stored in the header
	if (gwdata->ZERO_PADDED_FFT) index -= 7;
//BUG-Zero-padded not handled yet.  Need bigger gwnum header (if we allow storing ffted polys)
//BUG-current use of header location -10 conflicts with zero-padded use of location -10.

	// AVX512 implementation
	union {
		VDT	a;
		double	b[8];
	} x;
	if (index == -1) {			// First real-only value
		for (int j = 0; j < size; j++) {
			if (data[j] == NULL) continue;
			x.a = vec[j].real, data[j][0] = x.b[0];
//bug? return 0.0 instead of -11
//			x.a = vec[j].imag, data[j][-11] = x.b[0];
		}
	}
	else if (index == 0 && !gwdata->ALL_COMPLEX_FFT) {			// Second real-only value
		for (int j = 0; j < size; j++) {
			if (data[j] == NULL) continue;
			x.a = vec[j].real, data[j][7] = x.b[7], data[j][6] = x.b[6], data[j][5] = x.b[5], data[j][4] = x.b[4],
					   data[j][3] = x.b[3], data[j][2] = x.b[2], data[j][1] = x.b[1], data[j][8] = x.b[0];
//bug? return 0.0 instead of -10
			x.a = vec[j].imag, data[j][15] = x.b[7], data[j][14] = x.b[6], data[j][13] = x.b[5], data[j][12] = x.b[4],
					   data[j][11] = x.b[3], data[j][10] = x.b[2], data[j][9] = x.b[1];
		}
	}
	else {				// Process complex values where imaginary part is in the high half of the cache line
		unsigned long offset = cache_line_offset (gwdata, index * 2) / sizeof (VDT);
		for (int j = 0; j < size; j++) {
			if (data[j] == NULL) continue;
			((VDT *) data[j])[offset] = vec[j].real;
			((VDT *) data[j])[offset+1] = vec[j].imag;
		}
	}
}

#elif defined (FMA)				// FMA and AVX FFTs use the same memory layout

#define read_line_fma		read_line_avx
#define write_line_fma		write_line_avx
#define pick_pass_sizes_fma	pick_pass_sizes_avx
bool read_line_avx (gwhandle *gwdata, gwnum *data, int size, int index, CVDT *vec, int options, bool post_process_monics);
void write_line_avx (gwhandle *gwdata, gwnum *data, int size, int index, CVDT *vec, int options);
void pick_pass_sizes_avx (pmhandle *pmdata, unsigned int fftsize, unsigned int *p1size, unsigned int *p2size);

#elif defined (AVX)

bool read_line_avx (gwhandle *gwdata, gwnum *data, int size, int index, CVDT *vec, int options, bool post_process_monics)
{
	// If not post-processing monics and the input is a monic RLP, then output a 1+0i
	if (!post_process_monics && (options & (POLYMULT_INVEC1_MONIC | POLYMULT_INVEC2_MONIC)) && (options & (POLYMULT_INVEC1_RLP | POLYMULT_INVEC2_RLP))) {
//GW: check if GW_FFT1 should be used instead
		vec[0].real = broadcastsd (1.0);
		vec[0].imag = broadcastsd (0.0);
		vec++;
	}

	// If input is an RLP then read the data into the high half of vec.  We will duplicate the data as the last step.
	if (options & (POLYMULT_INVEC1_RLP | POLYMULT_INVEC2_RLP)) vec += (size-1);

	// Transformed real data ends up with 2 real values and FFTLEN/2-1 complex values.  Transformed complex data has FFTLEN/2 complex values.
	// Adjust index to handle the extra data point that must be returned for transformed real data.
	if (!gwdata->ALL_COMPLEX_FFT) index--;

	// Zero padded FFTs must also return for processing the 7 real values stored in the header
	if (gwdata->ZERO_PADDED_FFT) index -= 7;
//BUG-Zero-padded not handled yet.  Need bigger gwnum header (if we allow storing ffted polys)
//BUG-current use of header location -10 conflicts with zero-padded use of location -10.

	// If we've read all the gwnum data, return FALSE
	if (index * 8 >= (int) gwdata->FFTLEN) return (FALSE);

	// AVX implementation
	if (index == -1) {			// First real-only value
		for (int j = 0; j < size; j++) {
			if (data[j] == NULL) { memset (&vec[j], 0, sizeof (CVDT)); continue; }
			vec[j].real = _mm256_set_pd (0.0, 0.0, 0.0, data[j][0]);
//bug? return 0.0 instead of -11
//			vec[j].imag = _mm256_set_pd (0.0, 0.0, 0.0, data[j][-11]);
			vec[j].imag = _mm256_set_pd (0.0, 0.0, 0.0, 0.0);
		}
	}
	else if (index == 0 && !gwdata->ALL_COMPLEX_FFT) {			// Second real-only value
		for (int j = 0; j < size; j++) {
			if (data[j] == NULL) { memset (&vec[j], 0, sizeof (CVDT)); continue; }
			vec[j].real = _mm256_set_pd (data[j][3], data[j][2], data[j][1], data[j][4]);
//bug? return 0.0 instead of -10
//			vec[j].imag = _mm256_set_pd (data[j][7], data[j][6], data[j][5], data[j][-10]);
			vec[j].imag = _mm256_set_pd (data[j][7], data[j][6], data[j][5], 0.0);
		}
	}
	else {				// Process complex values where imaginary part is in the high half of the cache line
		unsigned long offset = cache_line_offset (gwdata, index) / sizeof (VDT);
		for (int j = 0; j < size; j++) {
			if (data[j] == NULL) { memset (&vec[j], 0, sizeof (CVDT)); continue; }
			ASSERTG (FFT_state (data[j]) == FULLY_FFTed);
			vec[j].real = ((VDT *) data[j])[offset];
			vec[j].imag = ((VDT *) data[j])[offset+1];
		}
	}

	// Negate coefficients if requested
	if (options & (POLYMULT_INVEC1_NEGATE | POLYMULT_INVEC2_NEGATE)) {
		VDT zero = broadcastsd (0.0);
		for (int j = 0; j < size; j++) {
			vec[j].real = subpd (zero, vec[j].real);
			vec[j].imag = subpd (zero, vec[j].imag);
		}
	}

	// Duplicate RLP data in reverse order
	if (options & (POLYMULT_INVEC1_RLP | POLYMULT_INVEC2_RLP)) {
		for (int j = 1; j < size; j++) vec[-j] = vec[j];
	}

	// If not post-processing monics and the input is monic, then output a 1+0i
	if (!post_process_monics && (options & (POLYMULT_INVEC1_MONIC | POLYMULT_INVEC2_MONIC))) {
//GW: check if GW_FFT1 should be used instead
		vec[size].real = broadcastsd (1.0);
		vec[size].imag = broadcastsd (0.0);
	}

	// Return success
	return (TRUE);
}

void write_line_avx (gwhandle *gwdata, gwnum *data, int size, int index, CVDT *vec, int options)
{
	// If either input was a monic RLP, back up the output vector pointer to the newly created entry by monic_line_adjustments
	if (((options & POLYMULT_INVEC1_MONIC) && (options & POLYMULT_INVEC1_RLP)) || ((options & POLYMULT_INVEC2_MONIC) && (options & POLYMULT_INVEC2_RLP))) vec--;

	// If both inputs are RLPs, then the output is an RLP too.  Skip writing half of the data.
	if ((options & POLYMULT_INVEC1_RLP) && (options & POLYMULT_INVEC2_RLP)) vec += (size-1);

	// Transformed real data ends up with 2 real values and FFTLEN/2-1 complex values.  Transformed complex data has FFTLEN/2 complex values.
	// Adjust index to handle the extra data point that must be returned for transformed real data.
	if (!gwdata->ALL_COMPLEX_FFT) index--;

	// Zero padded FFTs must also return for processing the 7 real values stored in the header
	if (gwdata->ZERO_PADDED_FFT) index -= 7;
//BUG-Zero-padded not handled yet.  Need bigger gwnum header (if we allow storing ffted polys)
//BUG-current use of header location -10 conflicts with zero-padded use of location -10.

	// AVX implementation
	union {
		VDT	a;
		double	b[4];
	} x;
	if (index == -1) {			// First real-only value
		for (int j = 0; j < size; j++) {
			if (data[j] == NULL) continue;
			x.a = vec[j].real, data[j][0] = x.b[0];
//bug? return 0.0 instead of -11
//			x.a = vec[j].imag, data[j][-11] = x.b[0];
		}
	}
	else if (index == 0 && !gwdata->ALL_COMPLEX_FFT) {			// Second real-only value
		for (int j = 0; j < size; j++) {
			if (data[j] == NULL) continue;
			x.a = vec[j].real, data[j][3] = x.b[3], data[j][2] = x.b[2], data[j][1] = x.b[1], data[j][4] = x.b[0];
//bug? return 0.0 instead of -10
			x.a = vec[j].imag, data[j][7] = x.b[3], data[j][6] = x.b[2], data[j][5] = x.b[1] /*, data[j][-10] = x.b[0]*/ ;
		}
	}
	else {				// Process complex values where imaginary part is in the high half of the cache line
		unsigned long offset = cache_line_offset (gwdata, index) / sizeof (VDT);
		for (int j = 0; j < size; j++) {
			if (data[j] == NULL) continue;
			((VDT *) data[j])[offset] = vec[j].real;
			((VDT *) data[j])[offset+1] = vec[j].imag;
		}
	}
}

#elif defined (SSE2)

bool read_line_sse2 (gwhandle *gwdata, gwnum *data, int size, int index, CVDT *vec, int options, bool post_process_monics)
{
	// If not post-processing monics and the input is a monic RLP, then output a 1+0i
	if (!post_process_monics && (options & (POLYMULT_INVEC1_MONIC | POLYMULT_INVEC2_MONIC)) && (options & (POLYMULT_INVEC1_RLP | POLYMULT_INVEC2_RLP))) {
//GW: check if GW_FFT1 should be used instead
		vec[0].real = broadcastsd (1.0);
		vec[0].imag = broadcastsd (0.0);
		vec++;
	}

	// If input is an RLP then read the data into the high half of vec.  We will duplicate the data as the last step.
	if (options & (POLYMULT_INVEC1_RLP | POLYMULT_INVEC2_RLP)) vec += (size-1);

	// Transformed real data ends up with 2 real values and FFTLEN/2-1 complex values.  Transformed complex data has FFTLEN/2 complex values.
	// Adjust index to handle the extra data point that must be returned for transformed real data.
	if (!gwdata->ALL_COMPLEX_FFT) index--;

	// Zero padded FFTs must also return for processing the 7 real values stored in the header
	if (gwdata->ZERO_PADDED_FFT) index -= 7;
//BUG-Zero-padded not handled yet.  Need bigger gwnum header (if we allow storing ffted polys)
//BUG-current use of header location -10 conflicts with zero-padded use of location -10.

	// If we've read all the gwnum data, return FALSE
	if (index * 4 >= (int) gwdata->FFTLEN) return (FALSE);

	// Various SSE2 FFT implementations store the imaginary part of a complex number either 8, 16, or 32 bytes above the real part
	int	imag_dist;
	if (gwdata->FFT_TYPE > 0 &&					// r4, r4delay, r4dwpn with 8-complex final levels
	    (gwdata->PASS2_SIZE == 1536 || gwdata->PASS2_SIZE == 2048 || gwdata->PASS2_SIZE == 2560 ||
	     gwdata->PASS2_SIZE == 4608 || gwdata->PASS2_SIZE == 6144 || gwdata->PASS2_SIZE == 7680 ||
	     gwdata->PASS2_SIZE == 8192 || gwdata->PASS2_SIZE == 10240 || gwdata->PASS2_SIZE == 12800)) imag_dist = 32;
	else if (gwdata->ALL_COMPLEX_FFT) imag_dist = 16;
	else if (gwdata->PASS2_SIZE == 0 && index < 4) imag_dist = 8;
	else if (gwdata->PASS2_SIZE != 0 && gwdata->FFT_TYPE == 0 && (index <= 1 || index == 4 || index == 5)) imag_dist = 8;	// hg ffts
	else imag_dist = 16;

	// SSE2 instructions implementation
	if (index == -1) {			// First real-only value
		for (int j = 0; j < size; j++) {
			if (data[j] == NULL) { memset (&vec[j], 0, sizeof (CVDT)); continue; }
			vec[j].real = _mm_set_pd (0.0, data[j][0]);
//bug? return 0.0 instead of -11 ---- fix when we support FFTing invec
			vec[j].imag = _mm_set_pd (0.0, 0.0);
		}
	}
	else if (index == 0 && !gwdata->ALL_COMPLEX_FFT) {			// Second real-only value
		for (int j = 0; j < size; j++) {
			if (data[j] == NULL) { memset (&vec[j], 0, sizeof (CVDT)); continue; }
			if (imag_dist == 8) {
				vec[j].real = _mm_set_pd (data[j][2], data[j][1]);
				vec[j].imag = _mm_set_pd (data[j][3], 0.0); //data[j][-10];
			} else if (imag_dist == 16) {
				vec[j].real = _mm_set_pd (data[j][1], data[j][2]);
//bug? return 0.0 instead of -10.  How are we going to implement FFTed poly read?  Its own routine? It's own FFT state? Compressed?
				vec[j].imag = _mm_set_pd (data[j][3], 0.0); //data[j][-10];
			} else {
				vec[j].real = _mm_set_pd (data[j][1], data[j][4]);
				vec[j].imag = _mm_set_pd (data[j][5], 0.0); //data[j][-10];
			}
		}
	}
	else {				// Process complex values where imaginary part is 8, 16, or 32 bytes above the real part
		unsigned long offset = cache_line_offset (gwdata, index / 2) / sizeof (VDT);
		for (int j = 0; j < size; j++) {
			if (data[j] == NULL) { memset (&vec[j], 0, sizeof (CVDT)); continue; }
			ASSERTG (FFT_state (data[j]) == FULLY_FFTed);
			if (imag_dist == 8) {
				vec[j].real = _mm_set_pd (data[j][index*4+2], data[j][index*4]);
				vec[j].imag = _mm_set_pd (data[j][index*4+3], data[j][index*4+1]);
			} else if (imag_dist == 16) {
				vec[j].real = ((VDT *) data[j])[offset+(index&1)*2];
				vec[j].imag = ((VDT *) data[j])[offset+(index&1)*2+1];
			} else {
				vec[j].real = ((VDT *) data[j])[offset+(index&1)];
				vec[j].imag = ((VDT *) data[j])[offset+(index&1)+2];
			}
		}
	}

	// Negate coefficients if requested
	if (options & (POLYMULT_INVEC1_NEGATE | POLYMULT_INVEC2_NEGATE)) {
		VDT zero = broadcastsd (0.0);
		for (int j = 0; j < size; j++) {
			vec[j].real = subpd (zero, vec[j].real);
			vec[j].imag = subpd (zero, vec[j].imag);
		}
	}

	// Duplicate RLP data in reverse order
	if (options & (POLYMULT_INVEC1_RLP | POLYMULT_INVEC2_RLP)) {
		for (int j = 1; j < size; j++) vec[-j] = vec[j];
	}

	// If not post-processing monics and the input is monic, then output a 1+0i
	if (!post_process_monics && (options & (POLYMULT_INVEC1_MONIC | POLYMULT_INVEC2_MONIC))) {
//GW: check if GW_FFT1 should be used instead
		vec[size].real = broadcastsd (1.0);
		vec[size].imag = broadcastsd (0.0);
	}

	// Return success
	return (TRUE);
}

void write_line_sse2 (gwhandle *gwdata, gwnum *data, int size, int index, CVDT *vec, int options)
{
	// If either input was a monic RLP, back up the output vector pointer to the newly created entry by monic_line_adjustments
	if (((options & POLYMULT_INVEC1_MONIC) && (options & POLYMULT_INVEC1_RLP)) || ((options & POLYMULT_INVEC2_MONIC) && (options & POLYMULT_INVEC2_RLP))) vec--;

	// If both inputs are RLPs, then the output is an RLP too.  Skip writing half of the data.
	if (options & POLYMULT_INVEC1_RLP && options & POLYMULT_INVEC2_RLP) vec += (size-1);

	// Transformed real data ends up with 2 real values and FFTLEN/2-1 complex values.  Transformed complex data has FFTLEN/2 complex values.
	// Adjust index to handle the extra data point that must be returned for transformed real data.
	if (!gwdata->ALL_COMPLEX_FFT) index--;

	// Zero padded FFTs must also return for processing the 7 real values stored in the header
	if (gwdata->ZERO_PADDED_FFT) index -= 7;
//BUG-Zero-padded not handled yet.  Need bigger gwnum header (if we allow storing ffted polys)
//BUG-current use of header location -10 conflicts with zero-padded use of location -10.

	// Various SSE2 FFT implementations store the imaginary part of a complex number either 8, 16, or 32 bytes above the real part
	int	imag_dist;
	if (gwdata->FFT_TYPE > 0 &&					// r4, r4delay, r4dwpn with 8-complex final levels
	    (gwdata->PASS2_SIZE == 1536 || gwdata->PASS2_SIZE == 2048 || gwdata->PASS2_SIZE == 2560 ||
	     gwdata->PASS2_SIZE == 4608 || gwdata->PASS2_SIZE == 6144 || gwdata->PASS2_SIZE == 7680 ||
	     gwdata->PASS2_SIZE == 8192 || gwdata->PASS2_SIZE == 10240 || gwdata->PASS2_SIZE == 12800)) imag_dist = 32;
	else if (gwdata->ALL_COMPLEX_FFT) imag_dist = 16;
	else if (gwdata->PASS2_SIZE == 0 && index < 4) imag_dist = 8;
	else if (gwdata->PASS2_SIZE != 0 && gwdata->FFT_TYPE == 0 && (index <= 1 || index == 4 || index == 5)) imag_dist = 8;	// hg ffts
	else imag_dist = 16;

	// SSE2 instructions implementation
	union {
		VDT	a;
		double	b[2];
	} x;
	if (index == -1) {			// First real-only value
		for (int j = 0; j < size; j++) {
			if (data[j] == NULL) continue;
			x.a = vec[j].real, data[j][0] = x.b[0];
//			data[j][-11] = vec[j].imag;  //GW:bug?
		}
	}
	else if (index == 0 && !gwdata->ALL_COMPLEX_FFT) {		// Second real-only value
		for (int j = 0; j < size; j++) {
			if (data[j] == NULL) continue;
			;
			if (imag_dist == 8) {
				x.a = vec[j].real, data[j][2] = x.b[1], data[j][1] = x.b[0];
				x.a = vec[j].imag, data[j][3] = x.b[1];
			} else if (imag_dist == 16) {
				x.a = vec[j].real, data[j][1] = x.b[1], data[j][2] = x.b[0];
				x.a = vec[j].imag, data[j][3] = x.b[1];
			} else {
				x.a = vec[j].real, data[j][1] = x.b[1], data[j][4] = x.b[0];
				x.a = vec[j].imag, data[j][5] = x.b[1];
			}
//bug? return 0.0 instead of -10
		}
	}
	else {				// Process complex values where imaginary part is 8 or 16 bytes above the real part
		unsigned long offset = cache_line_offset (gwdata, index / 2) / sizeof (VDT);
		for (int j = 0; j < size; j++) {
			if (data[j] == NULL) continue;
			if (imag_dist == 8) {
				x.a = vec[j].real, data[j][index*4+2] = x.b[1], data[j][index*4] = x.b[0];
				x.a = vec[j].imag, data[j][index*4+3] = x.b[1], data[j][index*4+1] = x.b[0];
			} else if (imag_dist == 16) {
				((VDT *) data[j])[offset+(index&1)*2] = vec[j].real;
				((VDT *) data[j])[offset+(index&1)*2+1] = vec[j].imag;
			} else {
				((VDT *) data[j])[offset+(index&1)] = vec[j].real;
				((VDT *) data[j])[offset+(index&1)+2] = vec[j].imag;
			}
		}
	}
}

#elif !defined (X86_64)				// x87 FFTs only supported in 32-bit mode

bool read_line_dbl (gwhandle *gwdata, gwnum *data, int size, int index, CVDT *vec, int options, bool post_process_monics)
{
	// If not post-processing monics and the input is a monic RLP, then output a 1+0i
	if (!post_process_monics && (options & (POLYMULT_INVEC1_MONIC | POLYMULT_INVEC2_MONIC)) && (options & (POLYMULT_INVEC1_RLP | POLYMULT_INVEC2_RLP))) {
//GW: check if GW_FFT1 should be used instead
		vec[0].real = 1.0;
		vec[0].imag = 0.0;
		vec++;
	}

	// If input is an RLP then read the data into the high half of vec.  We will duplicate the data as the last step.
	if (options & (POLYMULT_INVEC1_RLP | POLYMULT_INVEC2_RLP)) vec += (size-1);

	// Transformed real data ends up with 2 real values and FFTLEN/2-1 complex values.  Transformed complex data has FFTLEN/2 complex values.
	// Adjust index to handle the extra data point that must be returned for transformed real data.
	if (!gwdata->ALL_COMPLEX_FFT) index--;

	// Zero padded FFTs must also return for processing the 7 real values stored in the header
	if (gwdata->ZERO_PADDED_FFT) index -= 7;
//BUG-Zero-padded not handled yet.  Need bigger gwnum header (if we allow storing ffted polys)
//BUG-current use of header location -10 conflicts with zero-padded use of location -10.

	// If we've read all the gwnum data, return FALSE
	if (index * 2 >= (int) gwdata->FFTLEN) return (FALSE);

	// No special instructions implementation
	if (index == -1) {			// First real-only value
		for (int j = 0; j < size; j++) {
			if (data[j] == NULL) { memset (&vec[j], 0, sizeof (CVDT)); continue; }
			vec[j].real = data[j][0];
//bug? return 0.0 instead of -11 ---- fix when we support FFTing invec
			vec[j].imag = 0.0; //data[j][-11];
		}
	}
	else if (index == 0 && !gwdata->ALL_COMPLEX_FFT) {			// Second real-only value
		for (int j = 0; j < size; j++) {
			if (data[j] == NULL) { memset (&vec[j], 0, sizeof (CVDT)); continue; }
			if (gwdata->PASS2_SIZE == 0) vec[j].real = data[j][1];
			else vec[j].real = data[j][2];
//bug? return 0.0 instead of -10.  How are we going to implement FFTed poly read?  Its own routine? It's own FFT state? Compressed?
			vec[j].imag = 0.0; //data[j][-10];
		}
	}
	else {				// Process complex values where imaginary part is 8 or 16 bytes above the real part
		unsigned long offset = cache_line_offset (gwdata, index / 4) / sizeof (VDT);
		for (int j = 0; j < size; j++) {
			if (data[j] == NULL) { memset (&vec[j], 0, sizeof (CVDT)); continue; }
			ASSERTG (FFT_state (data[j]) == FULLY_FFTed);
// bug???  Index assuming what?  skip every other set of 4?  Examine ancient FFT sources!
			if (gwdata->PASS2_SIZE == 0 && index < 8 && !gwdata->ALL_COMPLEX_FFT) {
				vec[j].real = data[j][offset+(index&3)*2];
				vec[j].imag = data[j][offset+(index&3)*2+1];
			} else {
				vec[j].real = data[j][offset+(index&2)*2+(index&1)];
				vec[j].imag = data[j][offset+(index&2)*2+(index&1)+2];
			}
		}
	}

	// Negate coefficients if requested
	if (options & (POLYMULT_INVEC1_NEGATE | POLYMULT_INVEC2_NEGATE)) {
		for (int j = 0; j < size; j++) {
			vec[j].real = -vec[j].real;
			vec[j].imag = -vec[j].imag;
		}
	}

	// Duplicate RLP data in reverse order
	if (options & (POLYMULT_INVEC1_RLP | POLYMULT_INVEC2_RLP)) {
		for (int j = 1; j < size; j++) vec[-j] = vec[j];
	}

	// If not post-processing monics and the input is monic, then output a 1+0i
	if (!post_process_monics && (options & (POLYMULT_INVEC1_MONIC | POLYMULT_INVEC2_MONIC))) {
//GW: check if GW_FFT1 should be used instead
		vec[size].real = 1.0;
		vec[size].imag = 0.0;
	}

	// Return success
	return (TRUE);
}

void write_line_dbl (gwhandle *gwdata, gwnum *data, int size, int index, CVDT *vec, int options)
{
	// If either input was a monic RLP, back up the output vector pointer to the newly created entry by monic_line_adjustments
	if (((options & POLYMULT_INVEC1_MONIC) && (options & POLYMULT_INVEC1_RLP)) || ((options & POLYMULT_INVEC2_MONIC) && (options & POLYMULT_INVEC2_RLP))) vec--;

	// If both inputs are RLPs, then the output is an RLP too.  Skip writing half of the data.
	if (options & POLYMULT_INVEC1_RLP && options & POLYMULT_INVEC2_RLP) vec += (size-1);

	// Transformed real data ends up with 2 real values and FFTLEN/2-1 complex values.  Transformed complex data has FFTLEN/2 complex values.
	// Adjust index to handle the extra data point that must be returned for transformed real data.
	if (!gwdata->ALL_COMPLEX_FFT) index--;

	// Zero padded FFTs must also return for processing the 7 real values stored in the header
	if (gwdata->ZERO_PADDED_FFT) index -= 7;
//BUG-Zero-padded not handled yet.  Need bigger gwnum header (if we allow storing ffted polys)
//BUG-current use of header location -10 conflicts with zero-padded use of location -10.

	// No special instructions implementation
	if (index == -1) {			// First real-only value
		for (int j = 0; j < size; j++) {
			if (data[j] == NULL) continue;
			data[j][0] = vec[j].real;
//			data[j][-11] = vec[j].imag;  //GW:bug?
		}
	}
	else if (index == 0 && !gwdata->ALL_COMPLEX_FFT) {		// Second real-only value
		for (int j = 0; j < size; j++) {
			if (data[j] == NULL) continue;
			if (gwdata->PASS2_SIZE == 0) data[j][1] = vec[j].real;
			else data[j][2] = vec[j].real;
//			data[j][-10] = vec[j].imag;	//GW:bug?
		}
	}
	else {				// Process complex values where imaginary part is 8 or 16 bytes above the real part
		unsigned long offset = cache_line_offset (gwdata, index / 4) / sizeof (VDT);
		for (int j = 0; j < size; j++) {
			if (data[j] == NULL) continue;
			if (gwdata->PASS2_SIZE == 0 && index < 8 && !gwdata->ALL_COMPLEX_FFT) {
				data[j][offset+(index&3)*2] = vec[j].real;
				data[j][offset+(index&3)*2+1] = vec[j].imag;
			} else {
				data[j][offset+(index&2)*2+(index&1)] = vec[j].real;
				data[j][offset+(index&2)*2+(index&1)+2] = vec[j].imag;
			}
		}
	}
}

#endif


/*--------------------------------------------------------------------------
|	     Generalized routines (vector data type independent)
+-------------------------------------------------------------------------*/

#if defined (AVX512) || defined (FMA) || defined (AVX) || defined (SSE2) || !defined (X86_64)

void brute_line (
	CVDT	*invec1,		// First input poly
	int	invec1_size,		// Size of the first input polynomial
	CVDT	*invec2,		// Second input poly
	int	invec2_size,		// Size of the second input polynomial
	CVDT	*outvec,		// Output poly
	int	outvec_size,		// Either invec1_size+invec2_size-1 or size of gw_outvec if circular convolution
	int	options,		// Options passed to polymult
	gwnum	*gw_outvec)		// Gwnum array passed to polymult (so we can detect outputs that do not need to be computed)
{
	// AVX implementation.  Loop through output indices.
	for (int out = 0; out < invec1_size + invec2_size - 1; out++) {
		// Adjust output index for circular convolutions
		int adjusted_out = out;
		bool first_time = TRUE;
		if (out >= outvec_size) {
			adjusted_out = out - outvec_size;
			first_time = FALSE;
		}
		// See if this output needs to be computed
		if (gw_outvec != NULL && gw_outvec[adjusted_out] == NULL &&	//GW: Catch the RLP cases too - monic_adjust creates another output and RLP*RLP tosses half the data
		    ((options & POLYMULT_CIRCULAR) ||
		     ((!(options & POLYMULT_INVEC1_MONIC) || !(options & POLYMULT_INVEC1_RLP)) &&
		      (!(options & POLYMULT_INVEC2_MONIC) || !(options & POLYMULT_INVEC2_RLP))))) continue;
		// Get index into input vectors
		int in1 = out < invec2_size ? 0 : out - invec2_size + 1;	// Index into invec1
		int in2 = out - in1;						// Index into invec2
		// Loop through all the pairs to accumulate in this output
		for ( ; in1 < invec1_size && in2 >= 0; in1++, in2--) {
			CVDT x = invec1[in1];
			CVDT y = invec2[in2];
			CVDT tmp;
			tmp.real = fmsubpd (x.real, y.real, mulpd (x.imag, y.imag));
			tmp.imag = fmaddpd (x.real, y.imag, mulpd (x.imag, y.real));
			if (first_time) {
				outvec[adjusted_out] = tmp;
				first_time = FALSE;
			} else {
				cvadd (outvec[adjusted_out], outvec[adjusted_out], tmp);
			}
		}
	}
}

// A recursive Karatsuba O(N^1.565) implementation of polynomial multiplication
//	a	b	(input vector 1)
//	c	d	(input vector 2)
//-----------------   
//ac   ad+bc   bd
//
//(a+b)*(c+d) = ac + ad + bc + bd
// ad+bc = (a+b)*(c+d) - ac - bd
void karatsuba_line (
	pmhandle *pmdata,		// Handle for polymult library
	CVDT	*invec1,		// Coefficients of first input poly
	int	invec1_size,		// Size of the first input polynomial
	CVDT	*invec2,		// Coefficients of second input poly
	int	invec2_size,		// Size of the second input polynomial
	CVDT	*tmp,			// Temporary space
	CVDT	*outvec)		// Output poly coefficients
{
	int	i, a_size, b_size, c_size, d_size, a_plus_b_size, c_plus_d_size, bd_size, ac_size, a_plus_b_times_c_plus_d_size;
	CVDT	*a, *b, *c, *d, *ac, *bd, *a_plus_b, *c_plus_d, *a_plus_b_times_c_plus_d;

	// Handle end of recursion
	if (invec1_size == 1 || invec2_size == 1 || invec1_size + invec2_size < pmdata->KARAT_BREAK) {
		brute_line (invec1, invec1_size, invec2, invec2_size, outvec, invec1_size + invec2_size - 1, 0, NULL);
		return;
	}

	// To simplify code, always make invec2 the larger poly.
	if (invec1_size > invec2_size) {
		int tmpsize = invec1_size; invec1_size = invec2_size; invec2_size = tmpsize;
		CVDT *tmpvec = invec1; invec1 = invec2; invec2 = tmpvec;
	}

	// Split larger vector roughly in half.  For readability, create pointers to two halves.
	d_size = (invec2_size + 1) / 2;
	c_size = invec2_size - d_size;
	c = invec2 + d_size, d = invec2;

	// Handle case where vectors are of wildly unequal lengths
	if (d_size >= invec1_size) {
		// Compute bd
		bd = outvec;
		karatsuba_line (pmdata, invec1, invec1_size, d, d_size, tmp, bd);
		bd_size = invec1_size + d_size - 1;
		// Compute bc and store it in tmp
		CVDT *bc = tmp;
		int bc_size = invec1_size + c_size - 1;
		tmp += bc_size;
		karatsuba_line (pmdata, invec1, invec1_size, c, c_size, tmp, bc);
		// add/copy bc to output
		for (i = 0; i < invec1_size - 1; i++) cvadd (outvec[i + d_size], outvec[i + d_size], bc[i]);
		for ( ; i < bc_size; i++) outvec[i + d_size] = bc[i];
		return;
	}

	// Split smaller vector, the size of the two lower halves must be the same.  For readability, create pointers to two halves.
	b_size = d_size;
	a_size = invec1_size - b_size;
	a = invec1 + b_size, b = invec1;

	// Compute (a+b)
	a_plus_b = outvec;
	for (i = 0; i < a_size && i < b_size; i++) cvadd (a_plus_b[i], a[i], b[i]);
//	for ( ; i < a_size; i++) a_plus_b[i] = a[i];		// In case a_size > b_size
	for ( ; i < b_size; i++) a_plus_b[i] = b[i];		// In case b_size > a_size
	a_plus_b_size = i;
	// Compute (c+d)
	c_plus_d = &outvec[i];
	for (i = 0; i < c_size && i < d_size; i++) cvadd (c_plus_d[i], c[i], d[i]);
//	for ( ; i < c_size; i++) c_plus_d[i] = c[i];		// In case c_size > d_size
	for ( ; i < d_size; i++) c_plus_d[i] = d[i];		// In case d_size > c_size
	c_plus_d_size = i;
	// Compute (a+b) * (c+d) and store it in tmp
	a_plus_b_times_c_plus_d = tmp;
	a_plus_b_times_c_plus_d_size = a_plus_b_size + c_plus_d_size - 1;
	tmp += a_plus_b_times_c_plus_d_size;
	karatsuba_line (pmdata, a_plus_b, a_plus_b_size, c_plus_d, c_plus_d_size, tmp, a_plus_b_times_c_plus_d);

	// Compute bd
	bd = outvec;
	karatsuba_line (pmdata, b, b_size, d, d_size, tmp, bd);
	bd_size = b_size + d_size - 1;
	// Output the complex value 0+0i
	memset (&outvec[bd_size], 0, sizeof (CVDT));
	// Compute ac
	ac = &outvec[bd_size+1];
	karatsuba_line (pmdata, a, a_size, c, c_size, tmp, ac);
	ac_size = a_size + c_size - 1;

	// subtract ac from middle
	for (i = 0; i < ac_size; i++) cvsub (a_plus_b_times_c_plus_d[i], a_plus_b_times_c_plus_d[i], ac[i]);
	// subtract bd from middle
	for (i = 0; i < bd_size; i++) cvsub (a_plus_b_times_c_plus_d[i], a_plus_b_times_c_plus_d[i], bd[i]);
	// add middle to output
	for (i = 0; i < a_plus_b_times_c_plus_d_size; i++) cvadd (outvec[i + b_size], outvec[i + b_size], a_plus_b_times_c_plus_d[i]);
}

// Pick FFT pass sizes
#ifndef FMA
void pick_pass_sizes (
	pmhandle *pmdata,		// Handle for polymult library
	unsigned int fftsize,		// FFT size
	unsigned int *p1size,		// Returned pass 1 size
	unsigned int *p2size)		// Returned pass 2 size
{
	unsigned int powers_of_two;
	double	max_p2size, log_fftsize;

	// Figure out what power-of-two FFT size will fit in the L2 cache.  One complex number fits in one CVDT.
	max_p2size = log2 ((double) pmdata->L2_CACHE_SIZE * 1024.0 / (double) sizeof (CVDT));

	// If the FFT can be done in one pass, choose that option
	log_fftsize = log2 ((double) fftsize);
	if (log_fftsize <= max_p2size) { *p1size = 1; *p2size = fftsize; return; }

	// Make pass sizes roughly equal
	*p2size = (unsigned int) log_fftsize / 2;
	if (*p2size > (unsigned int) max_p2size) *p2size = (unsigned int) max_p2size;

	// Make sure any odd power of two is done in pass 2
	for (powers_of_two = 0; (fftsize & (1 << powers_of_two)) == 0; powers_of_two++);
	if (*p2size > powers_of_two) *p2size = powers_of_two;
	if ((powers_of_two & 1) == 0) {		// no radix-2/8 step
		if ((*p2size & 1) == 1) *p2size -= 1;
	} else {
		if ((*p2size & 1) == 0) *p2size -= 1;
	}
	*p1size = fftsize >> *p2size;
	*p2size = 1 << *p2size;
}
#endif

// Perform forward FFT, pointmul, inverse FFT in as few passes over memory as possible.
// This routine can be called several different ways:
#define	FORWARD_FFT_INVEC1	0x1	// Forward FFT invec1 required
#define	FORWARD_FFT_INVEC2	0x2	// Forward FFT invec2 required
#define FORWARD_FFT_ONLY	0x4	// No pointwise multiply and inverse FFT

void fft_line (
	pmhandle *pmdata,		// Handle for polymult library
	CVDT	*invec1,		// Array of complex values
	CVDT	*invec2,		// Array of complex values
	CVDT	*outvec,		// Array of complex values
	CVDT	*scratch,		// Scratch area for two-pass FFTs       
	unsigned int fft_size,		// The fft size
	int	options)
{
	unsigned int pass1_size, pass2_size;	// Parameters for creating a two-pass FFT
	unsigned int twiddle3_stride = 1;	// Stride through the radix-3 twiddles data
	unsigned int twiddle45_stride = 1;	// Stride through the radix-4/5 twiddles data
	unsigned int inv_twiddle3_stride;	// Starting stride for the inverse FFT radix-3 twiddles data
	unsigned int inv_twiddle45_stride;	// Starting stride for the inverse FFT radix-4/5 twiddles data
	unsigned int powers_of_two;
	VDT	inv_fft_size = broadcastsd (1.0 / (double) fft_size);

	generate_sincos (pmdata, fft_size);

	// Init some needed items
	for (powers_of_two = 0; (fft_size & (1 << powers_of_two)) == 0; powers_of_two++);		// Inverse FFT check for a radix-8 step
	for (inv_twiddle3_stride = 1, inv_twiddle45_stride = fft_size; inv_twiddle45_stride % 3 == 0; ) {// Calc first inverse FFT twiddle strides
		inv_twiddle3_stride *= 3;
		inv_twiddle45_stride /= 3;
	}

	// When the data set gets large it overflows the L2 cache.  This causes every radix-step to read & write data from main memory.  To combat
	// this we use two passes to complete the FFT, doing as much work as possible while FFT data is in the L2 cache.
	pick_pass_sizes (pmdata, fft_size, &pass1_size, &pass2_size);

	// Do two-pass forward FFT merged with a pointmul and a two-pass inverse FFT.
	// Pass 1 is first half of the forward FFT.  Pass 2 is second half of the forward FFT, pointwise multiply, and half of the inverse FFT.
	// Pass 3 is the remaining half of the inverse FFT.
	for (unsigned int pass = 1; pass <= 3; pass++) {
	    unsigned int pass_size;			// The size of this pass
	    unsigned int pass_stride;			// Distance between smallest elements in the pass
	    unsigned int num_interior_blks;		// Num blocks smaller than pass_stride (for 2-pass FFTs this equals pass2_size or 1)
	    unsigned int num_exterior_blks;		// Num blocks larger than pass_stride (for 2-pass FFTs this equals 1 or pass1_size)
	    unsigned int saved_twiddle3_stride = twiddle3_stride;
	    unsigned int saved_twiddle45_stride = twiddle45_stride;

	    // Skip pass 1 and 3 for small FFTs (perform a one-pass FFT)
	    if (pass1_size == 1 && pass != 2) continue;

	    // Init some pass-specific variables
	    if (pass == 1 || pass == 3) {
		pass_size = pass1_size;
		pass_stride = pass2_size;
		num_interior_blks = pass2_size;
		num_exterior_blks = 1;
	    } else {
		pass_size = pass2_size;
		pass_stride = 1;
		num_interior_blks = 1;
		num_exterior_blks = pass1_size;
	    }

	    // Do pass an exterior a block at a time, then do pass an interior a block at a time (we're prepared for future three-pass FFTs)
	    for (unsigned int blke = 0; blke < num_exterior_blks; blke++) {
	    for (unsigned int blki = 0; blki < num_interior_blks; blki++) {
		CVDT *src;				// FFT source data
		CVDT *dest;				// FFT destination data
		unsigned int size;			// The current group size to process
		unsigned int stride;			// Logical stride through the FFT (often the same as input and output stride through memory)
		unsigned int instride1;			// Unit stride through input FFT data in memory
		unsigned int outstride1;		// Unit stride through output FFT data in memory

// Forward FFT.  Loop over both input vectors.  Both may need a forward FFT.
		
		if (pass <= 2) {
		    for (int srcarg = 1; srcarg <= 2; srcarg++) {

			// Skip input vectors that do not need a forward FFT
			if (srcarg == 1 && !(options & FORWARD_FFT_INVEC1)) continue;
			if (srcarg == 2 && !(options & FORWARD_FFT_INVEC2)) continue;
//GW: If neither are to be FFTed, need to skip to pointmul

			// Pass 1 reads from an input vector into a scratch area to reduce large strides which CPU caches hate.  At end, result written
			// back to input vector.  Pass 2 reads/writes directly from/to an input vector.
			if (pass == 1) {
			    src = (srcarg == 1 ? invec1 : invec2) + blki;
			    dest = scratch;
			    instride1 = pass2_size;
			    outstride1 = 1;
			}
			else {
			    src = dest = (srcarg == 1 ? invec1 : invec2) + blke * pass2_size;
			    instride1 = outstride1 = 1;
			}

			// Forward FFT.  Start working with group size as the whole pass size, then work on progressively smaller group sizes.
			size = pass_size;

			// Reset twiddle strides each time we execute the pass
			twiddle3_stride = saved_twiddle3_stride;
			twiddle45_stride = saved_twiddle45_stride;

//GW: Assume all radix-3 and radix-5 work is in the first (or only) pass???

			// Radix-3 FFT
			// A 3-complex FFT is:
			// Res1:  (R1+R2+R3) + (I1+I2+I3)i
			// Res2:  (R1-.5R2-.5R3-.866I2+.866I3) + (I1-.5I2-.5I3+.866R2-.866R3)i
			// Res3:  (R1-.5R2-.5R3+.866I2-.866I3) + (I1-.5I2-.5I3-.866R2+.866R3)i
			for ( ; size % 3 == 0; size /= 3, twiddle3_stride *= 3, src = dest, instride1 = 1) {
			    stride = size / 3;
			    if (pass == 1 && stride == 1) dest = (srcarg == 1 ? invec1 : invec2) + blki, outstride1 = pass2_size;
			    VDT half = broadcastsd (0.5);
			    VDT _866 = broadcastsd (0.86602540378443864676372317075294);
			    for (unsigned int group = 0; group < pass_size; group += stride * 3) {
				for (unsigned int group_member = 0; group_member < stride; group_member++) {
					unsigned int twiddle_idx = ((group_member * pass_stride) + blki) * twiddle3_stride;
					CVDT a = src[instride1*(group + group_member)];
					CVDT b = src[instride1*(group + group_member + stride)];
					CVDT c = src[instride1*(group + group_member + 2*stride)];

					CVDT b_plus_c; cvadd (b_plus_c, b, c);				// b + c

#if !VFMA
					VDT real_tmp23a = fnmaddpd (half, b_plus_c.real, a.real);	// a - 0.5 * (b + c)
					VDT imag_tmp23a = fnmaddpd (half, b_plus_c.imag, a.imag);
					VDT real_tmp23b = mulpd (_866, subpd (b.imag, c.imag));		// 0.866 * (b - c)
					VDT imag_tmp23b = mulpd (_866, subpd (b.real, c.real));
					b.real = subpd (real_tmp23a, real_tmp23b);			// combine 23a and 23b
					b.imag = addpd (imag_tmp23a, imag_tmp23b);
					c.real = addpd (real_tmp23a, real_tmp23b);
					c.imag = subpd (imag_tmp23a, imag_tmp23b);
#else
					VDT real_tmp23a = fnmaddpd (half, b_plus_c.real, a.real);	// a - 0.5 * (b + c)
					VDT imag_tmp23a = fnmaddpd (half, b_plus_c.imag, a.imag);
					VDT real_tmp23b = subpd (b.imag, c.imag);			// 0.866 delayed * (b - c)
					VDT imag_tmp23b = subpd (b.real, c.real);
					b.real = fnmaddpd (_866, real_tmp23b, real_tmp23a);		// combine 23a and 23b
					b.imag = fmaddpd  (_866, imag_tmp23b, imag_tmp23a);
					c.real = fmaddpd  (_866, real_tmp23b, real_tmp23a);
					c.imag = fnmaddpd (_866, imag_tmp23b, imag_tmp23a);
#endif

					cvadd (a, a, b_plus_c);						// a + b + c
					dest[outstride1*(group + group_member)] = a;
					if (twiddle_idx == 0) {	// First group member's twiddles are 1+0i, skip twiddle multiply
						dest[outstride1*(group + group_member + stride)]   = b;
						dest[outstride1*(group + group_member + 2*stride)] = c;
					}
					else {			// Apply twiddles
						CVDT twiddle1; cvload (twiddle1, pmdata->twiddles1[2*twiddle_idx], pmdata->twiddles1[2*twiddle_idx+1]);
						cvmul     (dest[outstride1*(group + group_member + stride)], b, twiddle1);
						cvconjmul (dest[outstride1*(group + group_member + 2*stride)], c, twiddle1);
					}
				}
			    }
			}

			// Radix-5 FFT
			// R1= r1     +(r2+r5)     +(r3+r4)
			// R2= r1 +.309(r2+r5) -.809(r3+r4)    -.951(i2-i5) -.588(i3-i4)
			// R5= r1 +.309(r2+r5) -.809(r3+r4)    +.951(i2-i5) +.588(i3-i4)
			// R3= r1 -.809(r2+r5) +.309(r3+r4)    -.588(i2-i5) +.951(i3-i4)
			// R4= r1 -.809(r2+r5) +.309(r3+r4)    +.588(i2-i5) -.951(i3-i4)
			// I1= i1     +(i2+i5)     +(i3+i4)
			// I2= i1 +.309(i2+i5) -.809(i3+i4)    +.951(r2-r5) +.588(r3-r4)
			// I5= i1 +.309(i2+i5) -.809(i3+i4)    -.951(r2-r5) -.588(r3-r4)
			// I3= i1 -.809(i2+i5) +.309(i3+i4)    +.588(r2-r5) -.951(r3-r4)
			// I4= i1 -.809(i2+i5) +.309(i3+i4)    -.588(r2-r5) +.951(r3-r4)
			for ( ; size % 5 == 0; size /= 5, twiddle45_stride *= 5, src = dest, instride1 = 1) {
			    VDT _309 = broadcastsd (0.30901699437494742410229341718282);
			    VDT _809 = broadcastsd (0.80901699437494742410229341718282);
			    VDT _951 = broadcastsd (0.95105651629515357211643933337938);
			    VDT _588_951 = broadcastsd (0.61803398874989484820458683436564);
			    stride = size / 5;
			    if (pass == 1 && stride == 1) dest = (srcarg == 1 ? invec1 : invec2) + blki, outstride1 = pass2_size;
			    for (unsigned int group = 0; group < pass_size; group += stride * 5) {
				for (unsigned int group_member = 0; group_member < stride; group_member++) {
					unsigned int twiddle_idx = ((group_member * pass_stride) + blki) * twiddle45_stride;
					CVDT a = src[instride1*(group + group_member)];
					CVDT b = src[instride1*(group + group_member + stride)];
					CVDT c = src[instride1*(group + group_member + 2*stride)];
					CVDT d = src[instride1*(group + group_member + 3*stride)];
					CVDT e = src[instride1*(group + group_member + 4*stride)];

					CVDT b_plus_e; cvadd (b_plus_e, b, e);				// b + e
					CVDT c_plus_d; cvadd (c_plus_d, c, d);				// c + d
					CVDT b_minus_e; cvsub (b_minus_e, b, e);			// b - e
					CVDT c_minus_d; cvsub (c_minus_d, c, d);			// c - d

#if !VFMA
					VDT real_tmp25a = fmaddpd (_309, b_plus_e.real, fnmaddpd (_809, c_plus_d.real, a.real)); // a +.309*(b+e) -.809*(c+d)
					VDT imag_tmp25a = fmaddpd (_309, b_plus_e.imag, fnmaddpd (_809, c_plus_d.imag, a.imag));
					VDT real_tmp34a = fnmaddpd (_809, b_plus_e.real, fmaddpd (_309, c_plus_d.real, a.real)); // a -.809*(b+e) +.309*(c+d)
					VDT imag_tmp34a = fnmaddpd (_809, b_plus_e.imag, fmaddpd (_309, c_plus_d.imag, a.imag));
					VDT real_tmp25b = mulpd (_951, fmaddpd (_588_951, c_minus_d.imag, b_minus_e.imag));	// .951*(b-e) + 0.588*(c-d)
					VDT imag_tmp25b = mulpd (_951, fmaddpd (_588_951, c_minus_d.real, b_minus_e.real));
					VDT real_tmp34b = mulpd (_951, fmsubpd (_588_951, b_minus_e.imag, c_minus_d.imag));	// .588*(b-e) - 0.951*(c-d)
					VDT imag_tmp34b = mulpd (_951, fmsubpd (_588_951, b_minus_e.real, c_minus_d.real));
					b.real = subpd (real_tmp25a, real_tmp25b);				// combine 25a and 25b
					b.imag = addpd (imag_tmp25a, imag_tmp25b);
					e.real = addpd (real_tmp25a, real_tmp25b);
					e.imag = subpd (imag_tmp25a, imag_tmp25b);
					c.real = subpd (real_tmp34a, real_tmp34b);				// combine 34a and 34b
					c.imag = addpd (imag_tmp34a, imag_tmp34b);
					d.real = addpd (real_tmp34a, real_tmp34b);
					d.imag = subpd (imag_tmp34a, imag_tmp34b);
#else
					VDT real_tmp25a = fmaddpd (_309, b_plus_e.real, fnmaddpd (_809, c_plus_d.real, a.real)); // a +.309*(b+e) -.809*(c+d)
					VDT imag_tmp25a = fmaddpd (_309, b_plus_e.imag, fnmaddpd (_809, c_plus_d.imag, a.imag));
					VDT real_tmp34a = fnmaddpd (_809, b_plus_e.real, fmaddpd (_309, c_plus_d.real, a.real)); // a -.809*(b+e) +.309*(c+d)
					VDT imag_tmp34a = fnmaddpd (_809, b_plus_e.imag, fmaddpd (_309, c_plus_d.imag, a.imag));
					VDT real_tmp25b = fmaddpd (_588_951, c_minus_d.imag, b_minus_e.imag);			// delayed .951*(b-e) + 0.588*(c-d)
					VDT imag_tmp25b = fmaddpd (_588_951, c_minus_d.real, b_minus_e.real);
					VDT real_tmp34b = fmsubpd (_588_951, b_minus_e.imag, c_minus_d.imag);			// delayed .588*(b-e) - 0.951*(c-d)
					VDT imag_tmp34b = fmsubpd (_588_951, b_minus_e.real, c_minus_d.real);
					b.real = fnmaddpd (_951, real_tmp25b, real_tmp25a);				// combine 25a and 25b
					b.imag = fmaddpd  (_951, imag_tmp25b, imag_tmp25a);
					e.real = fmaddpd  (_951, real_tmp25b, real_tmp25a);
					e.imag = fnmaddpd (_951, imag_tmp25b, imag_tmp25a);
					c.real = fnmaddpd (_951, real_tmp34b, real_tmp34a);				// combine 34a and 34b
					c.imag = fmaddpd  (_951, imag_tmp34b, imag_tmp34a);
					d.real = fmaddpd  (_951, real_tmp34b, real_tmp34a);
					d.imag = fnmaddpd (_951, imag_tmp34b, imag_tmp34a);
#endif

					cvadd (a, a, b_plus_e);							// a + b + c + d + e
					cvadd (a, a, c_plus_d);
					dest[outstride1*(group + group_member)]   = a;
					if (twiddle_idx == 0) {	// First group member's twiddles are 1+0i, skip twiddle multiply
						dest[outstride1*(group + group_member + stride)]   = b;
						dest[outstride1*(group + group_member + 2*stride)] = c;
						dest[outstride1*(group + group_member + 3*stride)] = d;
						dest[outstride1*(group + group_member + 4*stride)] = e;
					}
					else {			// Apply twiddles
						CVDT twiddle1; cvload (twiddle1, pmdata->twiddles2[4*twiddle_idx], pmdata->twiddles2[4*twiddle_idx+1]);
						CVDT twiddle2; cvload (twiddle2, pmdata->twiddles2[4*twiddle_idx+2], pmdata->twiddles2[4*twiddle_idx+3]);
						cvmul     (dest[outstride1*(group + group_member + stride)], b, twiddle1);
						cvmul     (dest[outstride1*(group + group_member + 2*stride)], c, twiddle2);
						cvconjmul (dest[outstride1*(group + group_member + 3*stride)], d, twiddle2);
						cvconjmul (dest[outstride1*(group + group_member + 4*stride)], e, twiddle1);
					}
				}
			    }
			}

			// Radix-4 FFT
			for ( ; size != 8 && size % 4 == 0; size /= 4, twiddle45_stride *= 4, src = dest, instride1 = 1) {
			    stride = size / 4;
			    if (pass == 1 && stride == 1) dest = (srcarg == 1 ? invec1 : invec2) + blki, outstride1 = pass2_size;
			    for (unsigned int group = 0; group < pass_size; group += stride * 4) {
				for (unsigned int group_member = 0; group_member < stride; group_member++) {
					unsigned int twiddle_idx = ((group_member * pass_stride) + blki) * twiddle45_stride;
					CVDT a = src[instride1*(group + group_member)];
					CVDT b = src[instride1*(group + group_member + stride)];
					CVDT c = src[instride1*(group + group_member + 2*stride)];
					CVDT d = src[instride1*(group + group_member + 3*stride)];

					CVDT new_a; cvadd (new_a, a, c);			// a + c
					CVDT new_c; cvsub (new_c, a, c);			// a - c
					CVDT new_b; cvadd (new_b, b, d);			// b + d
					CVDT new_d; cvsub (new_d, b, d);			// b - d

					cvadd  (a, new_a, new_b);				// a + b
					cvsub  (b, new_a, new_b);				// a - b
					cviadd (c, new_c, new_d);				// c + id
					cvisub (d, new_c, new_d);				// c - id

					dest[outstride1*(group + group_member)] = a;
					if (twiddle_idx == 0) {	// First group member's twiddles are 1+0i, skip twiddle multiply
						dest[outstride1*(group + group_member + 1*stride)] = c;		// Do not bit-reverse outputs
						dest[outstride1*(group + group_member + 2*stride)] = b;
						dest[outstride1*(group + group_member + 3*stride)] = d;
					}
					else {			// Apply twiddles
						CVDT twiddle1; cvload (twiddle1, pmdata->twiddles2[4*twiddle_idx], pmdata->twiddles2[4*twiddle_idx+1]);
						CVDT twiddle2; cvload (twiddle2, pmdata->twiddles2[4*twiddle_idx+2], pmdata->twiddles2[4*twiddle_idx+3]);
						cvmul     (dest[outstride1*(group + group_member + 1*stride)], c, twiddle1);
						cvmul     (dest[outstride1*(group + group_member + 2*stride)], b, twiddle2);
						cvconjmul (dest[outstride1*(group + group_member + 3*stride)], d, twiddle1);
					}
				}
			    }
			}

			// Radix-8 FFT (if needed it is the last radix block done in pass 2)
			if (size == 8) {
			    VDT _707 = broadcastsd (0.70710678118654752440084436210485);
			    stride = 1;
			    for (unsigned int group = 0; group < pass_size; group += stride * 8) {
				for (unsigned int group_member = 0; group_member < stride; group_member++) {
					CVDT a = src[group + group_member];
					CVDT b = src[group + group_member + stride];
					CVDT c = src[group + group_member + 2*stride];
					CVDT d = src[group + group_member + 3*stride];
					CVDT e = src[group + group_member + 4*stride];
					CVDT f = src[group + group_member + 5*stride];
					CVDT g = src[group + group_member + 6*stride];
					CVDT h = src[group + group_member + 7*stride];

					CVDT new_a, new_b, new_c, new_d, new_e, new_f, new_g, new_h, twid_f, twid_h;
					cvadd (new_a, a, e);					// a + e
					cvsub (new_e, a, e);					// a - e
					cvadd (new_b, b, f);					// b + f
					cvsub (new_f, b, f);					// b - f
					cvadd (new_c, c, g);					// c + g
					cvsub (new_g, c, g);					// c - g
					cvadd (new_d, d, h);					// d + h
					cvsub (new_h, d, h);					// d - h

					cvadd  (a, new_a, new_c);				// a + c
					cvsub  (c, new_a, new_c);				// a - c
					cvadd  (b, new_b, new_d);				// b + d
					cvsub  (d, new_b, new_d);				// b - d
					cviadd (e, new_e, new_g);				// e + ig
					cvisub (g, new_e, new_g);				// e - ig
					cviadd (f, new_f, new_h);				// f + ih
					cvisub (h, new_f, new_h);				// f - ih

					cvadd  (new_a, a, b);					// a + b
					cvsub  (new_b, a, b);					// a - b
					cviadd (new_c, c, d);					// c + id
					cvisub (new_d, c, d);					// c - id
					twid_f.real = subpd (f.real, f.imag);			// twiddled f / .707
					twid_f.imag = addpd (f.real, f.imag);
					twid_h.real = addpd (h.real, h.imag);			// -twiddled h / .707
					twid_h.imag = subpd (h.imag, h.real);
#if !VFMA
					twid_f.real = mulpd (_707, twid_f.real);		// twiddled f
					twid_f.imag = mulpd (_707, twid_f.imag);
					twid_h.real = mulpd (_707, twid_h.real);		// -twiddled h
					twid_h.imag = mulpd (_707, twid_h.imag);
					cvadd (new_e, e, twid_f);				// e + (sqrthalf+sqrthalf*i)f
					cvsub (new_f, e, twid_f);				// e - (sqrthalf+sqrthalf*i)f
					cvsub (new_g, g, twid_h);				// g + (-sqrthalf+sqrthalf*i)h
					cvadd (new_h, g, twid_h);				// g - (-sqrthalf+sqrthalf*i)h
#else
					cvaddfm (new_e, e, _707, twid_f);			// e + (sqrthalf+sqrthalf*i)f
					cvsubfm (new_f, e, _707, twid_f);			// e - (sqrthalf+sqrthalf*i)f
					cvsubfm (new_g, g, _707, twid_h);			// g + (-sqrthalf+sqrthalf*i)h
					cvaddfm (new_h, g, _707, twid_h);			// g - (-sqrthalf+sqrthalf*i)h
#endif

					dest[group + group_member] = new_a;			// Do not bit reverse outputs
					dest[group + group_member + 4*stride] = new_b;
					dest[group + group_member + 2*stride] = new_c;
					dest[group + group_member + 6*stride] = new_d;
					dest[group + group_member + 1*stride] = new_e;
					dest[group + group_member + 5*stride] = new_f;
					dest[group + group_member + 3*stride] = new_g;
					dest[group + group_member + 7*stride] = new_h;
				}
			    }
			}
		    }
		}

// If pass 1 then there is no pointmul and inverse FFT to do.  If only FFTing the input(s) also skip pointmul and inverse FFT.

		if (pass == 1 || (options & FORWARD_FFT_ONLY)) continue;

// Pointwise multiply

		if (pass == 2) {
		    for (unsigned int j = 0; j < pass_size; j++) {
			CVDT x = invec1[blke * pass2_size + j];
			CVDT y = invec2[blke * pass2_size + j];
			outvec[blke * pass2_size + j].real = mulpd (fmsubpd (x.real, y.real, mulpd (x.imag, y.imag)), inv_fft_size);
			outvec[blke * pass2_size + j].imag = mulpd (fmaddpd (x.real, y.imag, mulpd (x.imag, y.real)), inv_fft_size);
		    }
		}

// Inverse FFT
		
		// Pass 2 reads/writes directly from/to output vector.
		// Pass 3 reads from output vector into a scratch area to reduce large strides.  At end, result written back to output vector.
		if (pass == 2) {
			src = dest = outvec + blke * pass2_size;
			instride1 = outstride1 = 1;
		} else {
			src = outvec + blki;
			dest = scratch;
			instride1 = pass2_size;
			outstride1 = 1;
		}

		// Inverse FFT starts with stride = 1 and moves to progressive smaller group sizes with larger strides.
		stride = 1;

		// In pass 2, use the starting twiddle strides.  In pass 3, restore the saved twiddle strides from the end of pass 2.
		if (pass == 2) {
			twiddle3_stride = inv_twiddle3_stride;
			twiddle45_stride = inv_twiddle45_stride;
		} else {
			twiddle3_stride = saved_twiddle3_stride;
			twiddle45_stride = saved_twiddle45_stride;
		}

		// Radix-8 FFT (if needed it is the last radix block done in pass 2)
		if (pass == 2 && powers_of_two & 1) {		// If odd number of powers of two do a radix-8 step
			VDT _707 = broadcastsd (0.70710678118654752440084436210485);
			twiddle45_stride /= 8;		// Next twiddle stride amount
			for (unsigned int group = 0; group < pass_size; group += stride * 8) {
				for (unsigned int group_member = 0; group_member < stride; group_member++) {
					CVDT a = src[group + group_member];
					CVDT b = src[group + group_member + stride];
					CVDT c = src[group + group_member + 2*stride];
					CVDT d = src[group + group_member + 3*stride];
					CVDT e = src[group + group_member + 4*stride];
					CVDT f = src[group + group_member + 5*stride];
					CVDT g = src[group + group_member + 6*stride];
					CVDT h = src[group + group_member + 7*stride];

					CVDT new_a, new_b, new_c, new_d, new_e, new_f, new_g, new_h, twid_f, twid_h;
					cvadd (new_a, a, e);					// a + e
					cvsub (new_e, a, e);					// a - e
					cvadd (new_b, b, f);					// b + f
					cvsub (new_f, b, f);					// b - f
					cvadd (new_c, c, g);					// c + g
					cvsub (new_g, c, g);					// c - g
					cvadd (new_d, d, h);					// d + h
					cvsub (new_h, d, h);					// d - h

					cvadd  (a, new_a, new_c);				// a + c
					cvsub  (c, new_a, new_c);				// a - c
					cvadd  (b, new_b, new_d);				// b + d
					cvsub  (d, new_b, new_d);				// b - d
					cvisub (e, new_e, new_g);				// e - ig
					cviadd (g, new_e, new_g);				// e + ig
					cvisub (f, new_f, new_h);				// f - ih
					cviadd (h, new_f, new_h);				// f + ih

					cvadd  (new_a, a, b);					// a + b
					cvsub  (new_b, a, b);					// a - b
					cvisub (new_c, c, d);					// c - id
					cviadd (new_d, c, d);					// c + id
					twid_f.real = addpd (f.real, f.imag);			// twiddled f / .707
					twid_f.imag = subpd (f.imag, f.real);
					twid_h.real = subpd (h.real, h.imag);			// -twiddled h / .707
					twid_h.imag = addpd (h.real, h.imag);
#if !VFMA
					twid_f.real = mulpd (_707, twid_f.real);		// twiddled f
					twid_f.imag = mulpd (_707, twid_f.imag);
					twid_h.real = mulpd (_707, twid_h.real);		// -twiddled h
					twid_h.imag = mulpd (_707, twid_h.imag);
					cvadd (new_e, e, twid_f);				// e + (sqrthalf-sqrthalf*i)f
					cvsub (new_f, e, twid_f);				// e - (sqrthalf-sqrthalf*i)f
					cvsub (new_g, g, twid_h);				// g + (-sqrthalf-sqrthalf*i)h
					cvadd (new_h, g, twid_h);				// g - (-sqrthalf-sqrthalf*i)h
#else
					cvaddfm (new_e, e, _707, twid_f);			// e + (sqrthalf-sqrthalf*i)f
					cvsubfm (new_f, e, _707, twid_f);			// e - (sqrthalf-sqrthalf*i)f
					cvsubfm (new_g, g, _707, twid_h);			// g + (-sqrthalf-sqrthalf*i)h
					cvaddfm (new_h, g, _707, twid_h);			// g - (-sqrthalf-sqrthalf*i)h
#endif

					dest[group + group_member] = new_a;			// Do not bit reverse outputs
					dest[group + group_member + 4*stride] = new_b;
					dest[group + group_member + 2*stride] = new_c;
					dest[group + group_member + 6*stride] = new_d;
					dest[group + group_member + 1*stride] = new_e;
					dest[group + group_member + 5*stride] = new_f;
					dest[group + group_member + 3*stride] = new_g;
					dest[group + group_member + 7*stride] = new_h;
				}
			}
			stride *= 8;			// Next stride amount
		}

		// Radix-4 FFT, one fewer complex multiply than radix-2, does more work while data is in registers
		for ( ; (size = pass_size / stride) % 4 == 0; stride *= 4, src = dest, instride1 = 1) {
			twiddle45_stride /= 4;
			if (pass == 3 && stride * 4 == pass_size) dest = outvec + blki, outstride1 = pass2_size;
			for (unsigned int group = 0; group < pass_size; group += stride * 4) {
				for (unsigned int group_member = 0; group_member < stride; group_member++) {
					unsigned int twiddle_idx = ((group_member * pass_stride) + blki) * twiddle45_stride;
					CVDT a = src[instride1*(group + group_member)];
					CVDT b = src[instride1*(group + group_member + stride)];
					CVDT c = src[instride1*(group + group_member + 2*stride)];
					CVDT d = src[instride1*(group + group_member + 3*stride)];

					// Special case first group member, no twiddles necessary
					CVDT twiddled_b, twiddled_c, twiddled_d;
					if (twiddle_idx == 0) {
						twiddled_b = b;
						twiddled_c = c;
						twiddled_d = d;
					}
					// Apply the CONJUGATE of the twiddle
					else {
						CVDT twiddle1; cvload (twiddle1, pmdata->twiddles2[4*twiddle_idx], pmdata->twiddles2[4*twiddle_idx+1]);
						CVDT twiddle2; cvload (twiddle2, pmdata->twiddles2[4*twiddle_idx+2], pmdata->twiddles2[4*twiddle_idx+3]);
						cvconjmul (twiddled_b, b, twiddle1);
						cvconjmul (twiddled_c, c, twiddle2);
						cvmul     (twiddled_d, d, twiddle1);
					}

					CVDT new_a; cvadd (new_a, a, twiddled_c);			// a + c
					CVDT new_c; cvsub (new_c, a, twiddled_c);			// a - c
					CVDT new_b; cvadd (new_b, twiddled_b, twiddled_d);		// b + d
					CVDT new_d; cvsub (new_d, twiddled_b, twiddled_d);		// b - d

					cvadd  (a, new_a, new_b);					// a + b
					cvsub  (b, new_a, new_b);					// a - b
					cvisub (c, new_c, new_d);					// c - id
					cviadd (d, new_c, new_d);					// c + id

					// Save results
					dest[outstride1*(group + group_member)]            = a;		// Do not bit reverse outputs
					dest[outstride1*(group + group_member + 2*stride)] = b;
					dest[outstride1*(group + group_member + 1*stride)] = c;
					dest[outstride1*(group + group_member + 3*stride)] = d;
				}
			}
		}

		// Radix-5 FFT
		// R1= r1     +(r2+r5)     +(r3+r4)
		// R2= r1 +.309(r2+r5) -.809(r3+r4)    +.951(i2-i5) +.588(i3-i4)
		// R5= r1 +.309(r2+r5) -.809(r3+r4)    -.951(i2-i5) -.588(i3-i4)
		// R3= r1 -.809(r2+r5) +.309(r3+r4)    +.588(i2-i5) -.951(i3-i4)
		// R4= r1 -.809(r2+r5) +.309(r3+r4)    -.588(i2-i5) +.951(i3-i4)
		// I1= i1     +(i2+i5)     +(i3+i4)
		// I2= i1 +.309(i2+i5) -.809(i3+i4)    -.951(r2-r5) -.588(r3-r4)
		// I5= i1 +.309(i2+i5) -.809(i3+i4)    +.951(r2-r5) +.588(r3-r4)
		// I3= i1 -.809(i2+i5) +.309(i3+i4)    -.588(r2-r5) +.951(r3-r4)
		// I4= i1 -.809(i2+i5) +.309(i3+i4)    +.588(r2-r5) -.951(r3-r4)
		for ( ; (size = pass_size / stride) % 5 == 0; stride *= 5, src = dest, instride1 = 1) {
			twiddle45_stride /= 5;
			if (pass == 3 && stride * 5 == pass_size) dest = outvec + blki, outstride1 = pass2_size;
			VDT _309 = broadcastsd (0.30901699437494742410229341718282);
			VDT _809 = broadcastsd (0.80901699437494742410229341718282);
			VDT _951 = broadcastsd (0.95105651629515357211643933337938);
			VDT _588_951 = broadcastsd (0.61803398874989484820458683436564);
			for (unsigned int group = 0; group < pass_size; group += stride * 5) {
				for (unsigned int group_member = 0; group_member < stride; group_member++) {
					unsigned int twiddle_idx = ((group_member * pass_stride) + blki) * twiddle45_stride;
					CVDT a = src[instride1*(group + group_member)];
					CVDT b = src[instride1*(group + group_member + stride)];
					CVDT c = src[instride1*(group + group_member + 2*stride)];
					CVDT d = src[instride1*(group + group_member + 3*stride)];
					CVDT e = src[instride1*(group + group_member + 4*stride)];

					// Special case first group member, no twiddles necessary
					CVDT twiddled_b, twiddled_c, twiddled_d, twiddled_e;
					if (twiddle_idx == 0) {
						twiddled_b = b;
						twiddled_c = c;
						twiddled_d = d;
						twiddled_e = e;
					}
					// Apply the CONJUGATE of the twiddle
					else {
						CVDT twiddle1; cvload (twiddle1, pmdata->twiddles2[4*twiddle_idx], pmdata->twiddles2[4*twiddle_idx+1]);
						CVDT twiddle2; cvload (twiddle2, pmdata->twiddles2[4*twiddle_idx+2], pmdata->twiddles2[4*twiddle_idx+3]);
						cvconjmul (twiddled_b, b, twiddle1);
						cvconjmul (twiddled_c, c, twiddle2);
						cvmul     (twiddled_d, d, twiddle2);
						cvmul     (twiddled_e, e, twiddle1);
					}

					CVDT b_plus_e; cvadd (b_plus_e, twiddled_b, twiddled_e);			// b + e
					CVDT c_plus_d; cvadd (c_plus_d, twiddled_c, twiddled_d);			// c + d
					CVDT b_minus_e; cvsub (b_minus_e, twiddled_b, twiddled_e);			// b - e
					CVDT c_minus_d; cvsub (c_minus_d, twiddled_c, twiddled_d);			// c - d

#if !VFMA
					VDT real_tmp25a = fmaddpd (_309, b_plus_e.real, fnmaddpd (_809, c_plus_d.real, a.real)); // a +.309*(b+e) -.809*(c+d)
					VDT imag_tmp25a = fmaddpd (_309, b_plus_e.imag, fnmaddpd (_809, c_plus_d.imag, a.imag));
					VDT real_tmp34a = fnmaddpd (_809, b_plus_e.real, fmaddpd (_309, c_plus_d.real, a.real)); // a -.809*(b+e) +.309*(c+d)
					VDT imag_tmp34a = fnmaddpd (_809, b_plus_e.imag, fmaddpd (_309, c_plus_d.imag, a.imag));
					VDT real_tmp25b = mulpd (_951, fmaddpd (_588_951, c_minus_d.imag, b_minus_e.imag));	// .951*(b-e) + 0.588*(c-d)
					VDT imag_tmp25b = mulpd (_951, fmaddpd (_588_951, c_minus_d.real, b_minus_e.real));
					VDT real_tmp34b = mulpd (_951, fmsubpd (_588_951, b_minus_e.imag, c_minus_d.imag));	// .588*(b-e) - 0.951*(c-d)
					VDT imag_tmp34b = mulpd (_951, fmsubpd (_588_951, b_minus_e.real, c_minus_d.real));
					b.real = addpd (real_tmp25a, real_tmp25b);				// combine 25a and 25b
					b.imag = subpd (imag_tmp25a, imag_tmp25b);
					e.real = subpd (real_tmp25a, real_tmp25b);
					e.imag = addpd (imag_tmp25a, imag_tmp25b);
					c.real = addpd (real_tmp34a, real_tmp34b);				// combine 34a and 34b
					c.imag = subpd (imag_tmp34a, imag_tmp34b);
					d.real = subpd (real_tmp34a, real_tmp34b);
					d.imag = addpd (imag_tmp34a, imag_tmp34b);
#else
					VDT real_tmp25a = fmaddpd (_309, b_plus_e.real, fnmaddpd (_809, c_plus_d.real, a.real)); // a +.309*(b+e) -.809*(c+d)
					VDT imag_tmp25a = fmaddpd (_309, b_plus_e.imag, fnmaddpd (_809, c_plus_d.imag, a.imag));
					VDT real_tmp34a = fnmaddpd (_809, b_plus_e.real, fmaddpd (_309, c_plus_d.real, a.real)); // a -.809*(b+e) +.309*(c+d)
					VDT imag_tmp34a = fnmaddpd (_809, b_plus_e.imag, fmaddpd (_309, c_plus_d.imag, a.imag));
					VDT real_tmp25b = fmaddpd (_588_951, c_minus_d.imag, b_minus_e.imag);			// delayed .951*(b-e) + 0.588*(c-d)
					VDT imag_tmp25b = fmaddpd (_588_951, c_minus_d.real, b_minus_e.real);
					VDT real_tmp34b = fmsubpd (_588_951, b_minus_e.imag, c_minus_d.imag);			// delayed .588*(b-e) - 0.951*(c-d)
					VDT imag_tmp34b = fmsubpd (_588_951, b_minus_e.real, c_minus_d.real);
					b.real = fmaddpd  (_951, real_tmp25b, real_tmp25a);				// combine 25a and 25b
					b.imag = fnmaddpd (_951, imag_tmp25b, imag_tmp25a);
					e.real = fnmaddpd (_951, real_tmp25b, real_tmp25a);
					e.imag = fmaddpd  (_951, imag_tmp25b, imag_tmp25a);
					c.real = fmaddpd  (_951, real_tmp34b, real_tmp34a);				// combine 34a and 34b
					c.imag = fnmaddpd (_951, imag_tmp34b, imag_tmp34a);
					d.real = fnmaddpd (_951, real_tmp34b, real_tmp34a);
					d.imag = fmaddpd  (_951, imag_tmp34b, imag_tmp34a);
#endif

					cvadd (a, a, b_plus_e);							// a + b + c + d + e
					cvadd (a, a, c_plus_d);

					// Save results
					dest[outstride1*(group + group_member)]            = a;
					dest[outstride1*(group + group_member + stride)]   = b;
					dest[outstride1*(group + group_member + 2*stride)] = c;
					dest[outstride1*(group + group_member + 3*stride)] = d;
					dest[outstride1*(group + group_member + 4*stride)] = e;
				}
			}
		}

		// Radix-3 FFT
		// A 3-complex inverse FFT is:
		// Res1: (R1+R2+R3) + (I1+I2+I3)i
		// Res2:  (R1-.5R2-.5R3+.866I2-.866I3) + (I1-.5I2-.5I3-.866R2+.866R3)i
		// Res3:  (R1-.5R2-.5R3-.866I2+.866I3) + (I1-.5I2-.5I3+.866R2-.866R3)i
		for ( ; (size = pass_size / stride) % 3 == 0; stride *= 3, src = dest, instride1 = 1) {
			twiddle3_stride /= 3;
			if (pass == 3 && stride * 3 == pass_size) dest = outvec + blki, outstride1 = pass2_size;
			VDT half = broadcastsd (0.5);
			VDT _866 = broadcastsd (0.86602540378443864676372317075294);
			for (unsigned int group = 0; group < pass_size; group += stride * 3) {
				for (unsigned int group_member = 0; group_member < stride; group_member++) {
					unsigned int twiddle_idx = ((group_member * pass_stride) + blki) * twiddle3_stride;
					CVDT a = src[instride1*(group + group_member)];
					CVDT b = src[instride1*(group + group_member + stride)];
					CVDT c = src[instride1*(group + group_member + 2*stride)];

					// Special case first group member, no twiddles necessary
					CVDT twiddled_b, twiddled_c;
					if (twiddle_idx == 0) {
						twiddled_b = b;
						twiddled_c = c;
					}
					// Apply the CONJUGATE of the twiddle
					else {
						CVDT twiddle1; cvload (twiddle1, pmdata->twiddles1[2*twiddle_idx], pmdata->twiddles1[2*twiddle_idx+1]);
						cvconjmul (twiddled_b, b, twiddle1);
						cvmul     (twiddled_c, c, twiddle1);
					}

					CVDT b_plus_c; cvadd (b_plus_c, twiddled_b, twiddled_c);			// b + c

#if !VFMA
					VDT real_tmp23a = fnmaddpd (half, b_plus_c.real, a.real);			// a - 0.5 * (b + c)
					VDT imag_tmp23a = fnmaddpd (half, b_plus_c.imag, a.imag);
					VDT real_tmp23b = mulpd (_866, subpd (twiddled_b.imag, twiddled_c.imag));	// 0.866 * (b - c)
					VDT imag_tmp23b = mulpd (_866, subpd (twiddled_b.real, twiddled_c.real));
					b.real = addpd (real_tmp23a, real_tmp23b);					// combine 23a and 23b
					b.imag = subpd (imag_tmp23a, imag_tmp23b);
					c.real = subpd (real_tmp23a, real_tmp23b);
					c.imag = addpd (imag_tmp23a, imag_tmp23b);
#else
					VDT real_tmp23a = fnmaddpd (half, b_plus_c.real, a.real);			// a - 0.5 * (b + c)
					VDT imag_tmp23a = fnmaddpd (half, b_plus_c.imag, a.imag);
					VDT real_tmp23b = subpd (twiddled_b.imag, twiddled_c.imag);			// 0.866 delayed * (b - c)
					VDT imag_tmp23b = subpd (twiddled_b.real, twiddled_c.real);
					b.real = fmaddpd  (_866, real_tmp23b, real_tmp23a);				// combine 23a and 23b
					b.imag = fnmaddpd (_866, imag_tmp23b, imag_tmp23a);
					c.real = fnmaddpd (_866, real_tmp23b, real_tmp23a);
					c.imag = fmaddpd  (_866, imag_tmp23b, imag_tmp23a);
#endif

					cvadd (a, a, b_plus_c);								// a + b + c
					dest[outstride1*(group + group_member)] = a;
					dest[outstride1*(group + group_member + stride)] = b;
					dest[outstride1*(group + group_member + 2*stride)] = c;
				}
			}
		}
	    }
	    }
	}
}

// outvec += input vector shifted the proper amount to create a monic polymult result
void monic_line_add_hi (CVDT *invec, int invec_size, CVDT *outvec, int outvec_size)
{
	// Add words starting from the highest word
	for (invec += invec_size-1, outvec += outvec_size-1; invec_size; invec--, outvec--, invec_size--) cvadd (outvec[0], outvec[0], invec[0]);
}

// outvec += input vector shifted to the low part of outvec
void monic_line_add_lo (CVDT *invec, int invec_size, CVDT *outvec, int outvec_size)
{
	for ( ; invec_size; invec++, outvec++, invec_size--) cvadd (outvec[0], outvec[0], invec[0]);
}

// Make adjustments to outvec for any monic input vectors
void monic_line_adjustments (CVDT *invec1, int invec1_size, CVDT *invec2, int invec2_size, CVDT *outvec, int outvec_size, int options)
{
	// If neither input is monic, no adjustments are necessary
	if (! (options & POLYMULT_INVEC1_MONIC) && ! (options & POLYMULT_INVEC2_MONIC)) return;

	// If either input is a monic RLP, create a new entry before the output vector calculated thusfar
	if (((options & POLYMULT_INVEC1_MONIC) && (options & POLYMULT_INVEC1_RLP)) ||
	    ((options & POLYMULT_INVEC2_MONIC) && (options & POLYMULT_INVEC2_RLP))) {
		outvec--;
		outvec_size++;
		memset (&outvec[0], 0, sizeof (CVDT));
	}

	// Create a new entry after the output vector calculated thusfar so that we can add to it
	outvec_size++;
	memset (&outvec[outvec_size-1], 0, sizeof (CVDT));

	// Adjust for an implied high coefficient of one -- add in the other poly in the high words
	if (options & POLYMULT_INVEC1_MONIC) monic_line_add_hi (invec2, invec2_size, outvec, outvec_size);
	if (options & POLYMULT_INVEC2_MONIC) monic_line_add_hi (invec1, invec1_size, outvec, outvec_size);

	// Adjust for an implied low coefficient of one -- add in the other poly to the low words
	if ((options & POLYMULT_INVEC1_MONIC) && (options & POLYMULT_INVEC1_RLP)) monic_line_add_lo (invec2, invec2_size, outvec, outvec_size);
	if ((options & POLYMULT_INVEC2_MONIC) && (options & POLYMULT_INVEC2_RLP)) monic_line_add_lo (invec1, invec1_size, outvec, outvec_size);
}

// Select from various polymult algorithms
void polymult_line (
	pmhandle *pmdata,		// Handle for polymult library
	gwnum	*gw_invec1,		// First input poly
	int	invec1_size,		// Size of the first input polynomial
	gwnum	*gw_invec2,		// Second input poly
	int	invec2_size,		// Size of the second input polynomial
	gwnum	*gw_outvec,		// Output poly of size invec1_size + invec2_size - 1
	int	outvec_size,		// Size of the output polynomial.  Should be invec1_size + invec2_size (less one if not monic).
	int	options)
{
	int	adjusted_invec1_size, adjusted_invec2_size, adjusted_outvec_size, fft_size;
	bool	must_brute, must_fft, post_process_monics;
	CVDT	*invec1, *invec2, *outvec, *tmpvec, *scratch;

	// Adjust sizes due to RLPs.  Calculate the output size assuming no monic inputs.
	adjusted_invec1_size = (options & POLYMULT_INVEC1_RLP) ? 2 * invec1_size - 1 : invec1_size;
	adjusted_invec2_size = (options & POLYMULT_INVEC2_RLP) ? 2 * invec2_size - 1 : invec2_size;
	adjusted_outvec_size = adjusted_invec1_size + adjusted_invec2_size - 1;

	// Sometimes we must (or it is beneficial) handle monic input polynomials as part of invec1/invec2 as opposed to a post-processing call to
	// monic_line_adjustments.  Find those cases and adjust vector sizes.
	if (options & POLYMULT_CIRCULAR) {
		// Warn about the dangerous 1*1 roundoff problem (one is NOT a random looking number).
		ASSERTG (!(options & POLYMULT_INVEC1_MONIC) || !(options & POLYMULT_INVEC2_MONIC) || gw_outvec[(adjusted_invec1_size + adjusted_invec2_size)%outvec_size] == NULL);
//GW: Find more cases where monic processing would be better done this way
		if (options & POLYMULT_INVEC1_MONIC) adjusted_invec1_size += (options & POLYMULT_INVEC1_RLP) ? 2 : 1;
		if (options & POLYMULT_INVEC2_MONIC) adjusted_invec2_size += (options & POLYMULT_INVEC2_RLP) ? 2 : 1;
		adjusted_outvec_size = outvec_size;
		post_process_monics = FALSE;
	} else
		post_process_monics = TRUE;

	// Polymult using AVX instructions.  Allocate memory.  For the output vector allocate two extra in case there are any monic RLP inputs.
	must_brute = (options & POLYMULT_CIRCULAR) && outvec_size != polymult_fft_size (outvec_size);
	must_fft = !must_brute && (options & (POLYMULT_FORWARD_FFT | POLYMULT_INVEC1_FFTED | POLYMULT_INVEC2_FFTED));
	if (must_brute || (!must_fft && adjusted_outvec_size < pmdata->KARAT_BREAK)) {
		invec1 = (CVDT *) aligned_malloc (adjusted_invec1_size * sizeof (CVDT), 64);
		invec2 = (CVDT *) aligned_malloc (adjusted_invec2_size * sizeof (CVDT), 64);
		outvec = (CVDT *) aligned_malloc ((adjusted_outvec_size + 2) * sizeof (CVDT), 64) + 1;
		tmpvec = NULL;
		scratch = NULL;
	}
	else if (!must_fft && !(options & POLYMULT_CIRCULAR) && adjusted_outvec_size < pmdata->FFT_BREAK) {
		invec1 = (CVDT *) aligned_malloc (adjusted_invec1_size * sizeof (CVDT), 64);
		invec2 = (CVDT *) aligned_malloc (adjusted_invec2_size * sizeof (CVDT), 64);
		outvec = (CVDT *) aligned_malloc ((adjusted_outvec_size + 2) * sizeof (CVDT), 64) + 1;
		tmpvec = (CVDT *) aligned_malloc ((adjusted_outvec_size + 32) * sizeof (CVDT), 64);	// Could need one extra for each recursion
		scratch = NULL;
	}
	else {
		unsigned int pass1_size, pass2_size;
		fft_size = (options & (POLYMULT_FORWARD_FFT | POLYMULT_CIRCULAR)) ? outvec_size :
			   (options & POLYMULT_INVEC1_FFTED) ? invec1_size :
			   (options & POLYMULT_INVEC2_FFTED) ? invec2_size : polymult_fft_size (adjusted_outvec_size);
		pick_pass_sizes (pmdata, fft_size, &pass1_size, &pass2_size);
		tmpvec = (CVDT *) aligned_malloc (fft_size * sizeof (CVDT), 64);		// FFT of invec1 goes here
		outvec = (CVDT *) aligned_malloc ((fft_size + 2) * sizeof (CVDT), 64) + 1;	// FFT of invec2 goes here as does outvec
		if (post_process_monics && (options & POLYMULT_INVEC2_MONIC)) invec1 = (CVDT *) aligned_malloc (fft_size * sizeof (CVDT), 64);
		else invec1 = tmpvec;
		if (post_process_monics && (options & POLYMULT_INVEC1_MONIC)) invec2 = (CVDT *) aligned_malloc (fft_size * sizeof (CVDT), 64);
		else invec2 = outvec;
		if (pass1_size > 1) scratch = (CVDT *) aligned_malloc (pass1_size * sizeof (CVDT), 64);
		else scratch = NULL;
	}

	// Grab four complex values from each gwnum coefficient.  Work that data set.  Then move onto remaining complex vals in gwnums.
	for ( ; ; ) {
		int line = next_line (pmdata);
		if (!read_line (pmdata->gwdata, gw_invec1, invec1_size, line, invec1, options & INVEC1_OPTIONS, post_process_monics)) break;
		read_line (pmdata->gwdata, gw_invec2, invec2_size, line, invec2, options & INVEC2_OPTIONS, post_process_monics);

		if (must_brute || (!must_fft && adjusted_outvec_size < pmdata->KARAT_BREAK)) {
			brute_line (invec1, adjusted_invec1_size, invec2, adjusted_invec2_size, outvec, adjusted_outvec_size, options, gw_outvec);
		}
		else if (!must_fft && !(options & POLYMULT_CIRCULAR) && adjusted_outvec_size < pmdata->FFT_BREAK) {
			karatsuba_line (pmdata, invec1, adjusted_invec1_size, invec2, adjusted_invec2_size, tmpvec, outvec);
		}
		else {
			ASSERTG (!(options & POLYMULT_INVEC1_FFTED));	// Not supported yet
			ASSERTG (!(options & POLYMULT_INVEC2_FFTED));	// Not supported yet
			if (invec1 != tmpvec) memcpy (tmpvec, invec1, adjusted_invec1_size * sizeof (CVDT));
			memset (tmpvec + adjusted_invec1_size, 0, (fft_size - adjusted_invec1_size) * sizeof (CVDT));
//			if (options & POLYMULT_FORWARD_FFT) {
//				write_line (pmdata->gwdata, gw_outvec, fft_size, line, tmpvec, 0);
//				continue;
//			}
			if (invec2 != outvec) memcpy (outvec, invec2, adjusted_invec2_size * sizeof (CVDT));
			memset (outvec + adjusted_invec2_size, 0, (fft_size - adjusted_invec2_size) * sizeof (CVDT));
			fft_line (pmdata, tmpvec, outvec, outvec, scratch, fft_size, FORWARD_FFT_INVEC1 | FORWARD_FFT_INVEC2);
		}
		if (post_process_monics) monic_line_adjustments (invec1, adjusted_invec1_size, invec2, adjusted_invec2_size, outvec, adjusted_outvec_size, options);
		write_line (pmdata->gwdata, gw_outvec, outvec_size, line, outvec, options & POLYMULT_CIRCULAR ? 0 : options);
	}

	// Free memory
	if (invec1 != tmpvec) aligned_free (invec1);
	if (invec2 != outvec) aligned_free (invec2);
	aligned_free (outvec - 1);
	aligned_free (tmpvec);
	aligned_free (scratch);
}

#endif


/*--------------------------------------------------------------------------
|			Polymult helper threads
+-------------------------------------------------------------------------*/

#if !defined (SSE2) && !defined (AVX) && !defined (FMA) && !defined (AVX512)

/* Entry point for helper thread to do part of the polymult */
void polymult_thread (void *arg)
{
	pmhandle *pmdata = (pmhandle *) arg;
	gwhandle cloned_gwdata;

/* Generate a unique helper thread number (so that callback routine to uniquely identify helper threads) */
/* Call gwnum's optional user provided callback routine so that the caller can set the thread's priority and affinity */

	gwmutex_lock (&pmdata->poly_mutex);
	int thread_num = ++pmdata->next_thread_num;
	gwmutex_unlock (&pmdata->poly_mutex);
	if (pmdata->gwdata->thread_callback != NULL)
		(*pmdata->gwdata->thread_callback) (thread_num, 20, pmdata->gwdata->thread_callback_data);

/* Clone the gwdata for user-defined callbacks to use */

	if (pmdata->helper_callback != NULL) gwclone (&cloned_gwdata, pmdata->gwdata);

/* Loop doing polymult work */

	while (!pmdata->termination_in_progress) {
		// After a user-defined callback completes, we must sync up all the stats and counters generated by the cloned gwdatas
		if (pmdata->helpers_sync_clone_stats) {
			gwmutex_lock (&pmdata->poly_mutex);
			gwclone_merge_stats (&cloned_gwdata);
			gwmutex_unlock (&pmdata->poly_mutex);
		}
		// Call user-defined routine with a cloned gwdata so the user can do multi-threaded gwnum work
		else if (!pmdata->helpers_doing_polymult)
			(*pmdata->helper_callback) (thread_num, &cloned_gwdata, pmdata->helper_callback_data);
		// Optimized routines for AVX512 instructions
		else if (pmdata->gwdata->cpu_flags & CPU_AVX512F)
			polymult_line_avx512 (pmdata, pmdata->invec1, pmdata->invec1_size, pmdata->invec2, pmdata->invec2_size,
					      pmdata->outvec, pmdata->outvec_size, pmdata->options);
		// Optimized routines for FMA3 instructions
		else if (pmdata->gwdata->cpu_flags & CPU_FMA3)
			polymult_line_fma (pmdata, pmdata->invec1, pmdata->invec1_size, pmdata->invec2, pmdata->invec2_size,
					   pmdata->outvec, pmdata->outvec_size, pmdata->options);
		// Optimized routines for AVX instructions
		else if (pmdata->gwdata->cpu_flags & CPU_AVX)
			polymult_line_avx (pmdata, pmdata->invec1, pmdata->invec1_size, pmdata->invec2, pmdata->invec2_size,
					   pmdata->outvec, pmdata->outvec_size, pmdata->options);
		// Optimized routines for SSE2 instructions
		else if (pmdata->gwdata->cpu_flags & CPU_SSE2)
			polymult_line_sse2 (pmdata, pmdata->invec1, pmdata->invec1_size, pmdata->invec2, pmdata->invec2_size,
					    pmdata->outvec, pmdata->outvec_size, pmdata->options);
		// Default routines using no vector instructions
#ifndef X86_64
		else
			polymult_line_dbl (pmdata, pmdata->invec1, pmdata->invec1_size, pmdata->invec2, pmdata->invec2_size,
					   pmdata->outvec, pmdata->outvec_size, pmdata->options);
#endif

		// When all helpers are finished, signal that all helpers are done and can now wait on the work_to_do event.
		// Otherwise wait for all helpers to finish
		gwmutex_lock (&pmdata->poly_mutex);
		int num_active = --pmdata->helpers_active;
		gwmutex_unlock (&pmdata->poly_mutex);
		if (num_active == 0) {
			gwevent_reset (&pmdata->work_to_do);
			gwevent_signal (&pmdata->helpers_done);
		} else
			gwevent_wait (&pmdata->helpers_done, 0);

		// All helpers are done.  When all helpers wake up it will be safe to signal the main thread.
		gwmutex_lock (&pmdata->poly_mutex);
		int num_waiting = ++pmdata->helpers_waiting_work;
		gwmutex_unlock (&pmdata->poly_mutex);
		if (num_waiting == pmdata->num_threads - 1) gwevent_signal (&pmdata->main_can_wakeup);

/* Wait for more work (or termination) */

		gwevent_wait (&pmdata->work_to_do, 0);
	}

/* Call gwnum's optional user provided callback routine so that the caller can do any necessary cleanup */

	if (pmdata->gwdata->thread_callback != NULL)
		(*pmdata->gwdata->thread_callback) (thread_num, 21, pmdata->gwdata->thread_callback_data);

/* Free the cloned gwdata that user-defined callbacks used */

	if (pmdata->helper_callback != NULL) gwdone (&cloned_gwdata);
}

/* This routine launches the polymult helper threads.  The polymult library uses this routine to do multithreading, users can too! */
void polymult_launch_helpers (
	pmhandle *pmdata)		// Handle for polymult library
{
	if (pmdata->num_threads == 1) return;
	// Set counters necessary for coordinating helper threads
	pmdata->helpers_active = pmdata->num_threads - 1;	// When this count reaches 0, it is safe to signal helpers_done
	pmdata->helpers_waiting_work = 0;			// When this reaches num_threads-1, it is safe for main thread to wake up
	// "Neuter" gwnum to single-threaded for user-definied callback
	if (!pmdata->helpers_doing_polymult && !pmdata->helpers_sync_clone_stats) {
		pmdata->saved_gwdata_num_threads = gwget_num_threads (pmdata->gwdata);
		gwset_num_threads (pmdata->gwdata, 1);
	}
	// Launch helpers for the first time
	if (pmdata->thread_ids == NULL) {
		pmdata->thread_ids = (gwthread *) malloc (pmdata->num_threads * sizeof (gwthread));
		gwmutex_init (&pmdata->poly_mutex);
		gwevent_init (&pmdata->work_to_do);
		gwevent_init (&pmdata->helpers_done);
		gwevent_init (&pmdata->main_can_wakeup);
		gwevent_signal (&pmdata->work_to_do);
		pmdata->next_thread_num = 0;
		for (int i = 1; i < pmdata->num_threads; i++) gwthread_create_waitable (&pmdata->thread_ids[i], &polymult_thread, (void *) pmdata);
	}
	// Reactivate helpers.  They're waiting for work to do.
	else {
		gwmutex_lock (&pmdata->poly_mutex);
		gwevent_reset (&pmdata->helpers_done);
		gwevent_reset (&pmdata->main_can_wakeup);
		gwevent_signal (&pmdata->work_to_do);		// Start all helper threads
		gwmutex_unlock (&pmdata->poly_mutex);
	}
}

/* This routine waits for the launched polymult helper threads to finish.  This is called by the main thread when it has run out of work to do. */
void polymult_wait_on_helpers (
	pmhandle *pmdata)		// Handle for polymult library
{
	if (pmdata->num_threads == 1) return;
	// Wait for helpers to end
	gwevent_wait (&pmdata->main_can_wakeup, 0);
	// User-defined helpers require syncing up cloned gwdata stats and restoring main thread's gwnum's multithreading capabilities
	if (!pmdata->helpers_doing_polymult) {
		pmdata->helpers_sync_clone_stats = TRUE;
		polymult_launch_helpers (pmdata);
		gwevent_wait (&pmdata->main_can_wakeup, 0);
		pmdata->helpers_sync_clone_stats = FALSE;
		gwset_num_threads (pmdata->gwdata, pmdata->saved_gwdata_num_threads);
	}
}

/*--------------------------------------------------------------------------
|			Polymult main routine
+-------------------------------------------------------------------------*/

/* Multiply two polynomials producing FFTed gwnums as output.  The output gwnums need to be normalized before use in further gwnum operations. */
/* It is safe to use input gwnums in the output vector.  For a normal polynomial multiply outvec_size must be invec1_size + invec2_size - 1. */
/* For monic polynomial multiply, the leading input coefficients of 1 are omitted as is the leading 1 output coefficient -- thus */
/* outvec_size must be invec1_size + invec2_size. */
void polymult (
	pmhandle *pmdata,		// Handle for polymult library
	gwnum	*invec1,		// First input poly
	int	invec1_size,		// Size of the first input polynomial (or fft size if POLYMULT_INVEC1_FFTED)
	gwnum	*invec2,		// Second input poly
	int	invec2_size,		// Size of the second input polynomial (or fft size if POLYMULT_INVEC2_FFTED)
	gwnum	*outvec,		// Output poly
	int	outvec_size,		// Size of the output polynomial (or fft_size if POLYMULT_FORWARD_FFT or POLYMULT_CIRCULAR)
	int	options)
{
	// Forward FFT not allowed on already FFTed polynomial
	ASSERTG (! ((options & POLYMULT_FORWARD_FFT) && (options & POLYMULT_INVEC1_FFTED)));
	// An already FFTed input cannot be multiplied with a monic polynomial
	ASSERTG (! ((options & POLYMULT_INVEC1_FFTED) && !(options & POLYMULT_INVEC2_MONIC)));
	ASSERTG (! ((options & POLYMULT_INVEC2_FFTED) && !(options & POLYMULT_INVEC1_MONIC)));
	// An already FFTed input must be at least as big as the output vector
	ASSERTG (! ((options & POLYMULT_INVEC1_FFTED) && invec1_size >= outvec_size));
	ASSERTG (! ((options & POLYMULT_INVEC2_FFTED) && invec2_size >= outvec_size));
	// Verify output vector size when we know the input vector sizes (no pre-FFTed inputs)
#define insz1	(invec1_size + !!(options & POLYMULT_INVEC1_MONIC))		// Adjust size for monic poly
#define insz2	(invec2_size + !!(options & POLYMULT_INVEC2_MONIC))
#define insz1a	((options & POLYMULT_INVEC1_RLP) ? 2*insz1-1 : insz1)		// Adjust size for RLP poly
#define insz2a	((options & POLYMULT_INVEC2_RLP) ? 2*insz2-1 : insz2)
#define outsz	(insz1a + insz2a - 1)						// Output poly size
#define outsz1	((options & POLYMULT_INVEC1_RLP && options & POLYMULT_INVEC2_RLP) ? (outsz + 1) / 2 : outsz)	// Output size adjusted for RLP
#define outsz2	(outsz1 - !!(options & POLYMULT_INVEC1_MONIC && options & POLYMULT_INVEC2_MONIC))		// Output size adjusted for monic
	ASSERTG ((options & (POLYMULT_FORWARD_FFT | POLYMULT_CIRCULAR | POLYMULT_INVEC1_FFTED | POLYMULT_INVEC2_FFTED) || outvec_size == outsz2));
#undef insz1
#undef insz2
#undef insz1a
#undef insz2a
#undef outsz
#undef outsz1
#undef outsz2

	// FFT inputs if needed
	if (!(options & POLYMULT_INVEC1_FFTED)) for (int j = 0; j < invec1_size; j++) if (invec1[j] != NULL) gwfft (pmdata->gwdata, invec1[j], invec1[j]);
	if (!(options & POLYMULT_INVEC2_FFTED)) for (int j = 0; j < invec2_size; j++) if (invec2[j] != NULL) gwfft (pmdata->gwdata, invec2[j], invec2[j]);

	// Prepare for polymult in parallel
	pmdata->invec1 = invec1;
	pmdata->invec1_size = invec1_size;
	pmdata->invec2 = invec2;
	pmdata->invec2_size = invec2_size;
	pmdata->outvec = outvec;
	pmdata->outvec_size = outvec_size;
	pmdata->options = options;
	pmdata->next_line = 0;

	// Fire up each helper thread
	pmdata->helpers_doing_polymult = TRUE;			// The helpers are doing polymult work, not user defined work
	polymult_launch_helpers (pmdata);

	// This thread works on the polymult too.  Select which optimized routines to use.
	if (pmdata->gwdata->cpu_flags & CPU_AVX512F)
		polymult_line_avx512 (pmdata, invec1, invec1_size, invec2, invec2_size, outvec, outvec_size, options);
	else if (pmdata->gwdata->cpu_flags & CPU_FMA3)
		polymult_line_fma (pmdata, invec1, invec1_size, invec2, invec2_size, outvec, outvec_size, options);
	else if (pmdata->gwdata->cpu_flags & CPU_AVX)
		polymult_line_avx (pmdata, invec1, invec1_size, invec2, invec2_size, outvec, outvec_size, options);
	else if (pmdata->gwdata->cpu_flags & CPU_SSE2)
		polymult_line_sse2 (pmdata, invec1, invec1_size, invec2, invec2_size, outvec, outvec_size, options);
#ifndef X86_64
	else
		polymult_line_dbl (pmdata, invec1, invec1_size, invec2, invec2_size, outvec, outvec_size, options);
#endif

	// Wait for all helper threads to complete
	polymult_wait_on_helpers (pmdata);
	pmdata->helpers_doing_polymult = FALSE;			// The helpers are not doing polymult work, user is free to have the threads do user work

	// Figure out which output entries (if any) need to get 1.0 added to the result because a monic poly is multiplied by a monic RLP poly.
	// We used to add the FFT(1) in monic_adjustments but this led to excessive round off error as FFT(1) is NOT a random signal.  It only
	// affects the least significant bits and is apt to get drowned out by the many multiplied signals.

	int	addin1 = -1;
	int	addin2 = -1;
	if ((options & POLYMULT_INVEC1_MONIC) && (options & POLYMULT_INVEC2_MONIC) && (options & POLYMULT_INVEC2_RLP)) addin1 = outvec_size - (2 * invec2_size);
	if ((options & POLYMULT_INVEC2_MONIC) && (options & POLYMULT_INVEC1_MONIC) && (options & POLYMULT_INVEC1_RLP)) addin2 = outvec_size - (2 * invec1_size);
	if (options & POLYMULT_INVEC2_MONIC_TIMES_MONIC_RLP_OK) addin2 = invec2_size;
			
	// Unfft all output coefficients
	if (!(options & POLYMULT_FORWARD_FFT)) {
		for (int j = 0; j < outvec_size; j++) {
			if (outvec[j] == NULL) continue;
			FFT_state (outvec[j]) = FULLY_FFTed;
			if ((options & POLYMULT_NO_UNFFT) && j != addin1 && j != addin2) continue;
			gwunfft2 (pmdata->gwdata, outvec[j], outvec[j], (j != addin1 && j != addin2 && (options & POLYMULT_STARTNEXTFFT)) ? GWMUL_STARTNEXTFFT : 0);
			if (j == addin1) gwaddsmall (pmdata->gwdata, outvec[j], 1);
			if (j == addin2) gwaddsmall (pmdata->gwdata, outvec[j], 1);
		}
	} else {
		// Mark all results poly FFTed
//		for (int j = 0; j < outvec_size; j++)
//GW - set unnorm count too? Don't need to do that if always unffting.
//GW - have a special FFT state for POLYMULT_FORWARD_FFT?
//		ASSERTG (outvec[j] != NULL);
//		FFT_state (outvec[j]) = FULLY_FFTed;
	}
}

#endif

