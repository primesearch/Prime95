/*----------------------------------------------------------------------
| gwnum.h
|
| This file contains the headers and definitions that are used in the
| multi-precision IBDWT arithmetic routines.  That is, all routines
| that deal with the gwnum data type.
|
| Gwnums are great for applications that do a lot of multiplies modulo
| a number.  Only Intel x86-platforms are supported.  Add and subtract
| are also pretty fast.
|
| Gwnums are not suited to applications that need to convert to and from
| binary frequently or need to change the modulus frequently.
| 
|  Copyright 2002-2005 Just For Fun Software, Inc.
|  All Rights Reserved.
+---------------------------------------------------------------------*/

#ifndef _GWNUM_H
#define _GWNUM_H

/* This is a C library.  If used in a C++ program, don't let the C++ */
/* compiler mangle names. */

#ifdef __cplusplus
extern "C" {
#endif

/* Include common definitions */

#include "common.h"

/* Handle the difference between the naming conventions in */
/* C compilers.  We need to do this for global variables that */
/* referenced by the assembly routines.  Most non-Windows systems */
/* should #define ADD_UNDERSCORES before including this file. */

#ifdef ADD_UNDERSCORES
#include "gwrename.h"
#endif

/* The gwnum data type.  A gwnum points to an array of doubles - the */
/* FFT data.  In practice, there is data stored before the doubles. */
/* See the internals section below if you really must know. */

typedef double *gwnum;

/*---------------------------------------------------------------------+
|                     SETUP AND TERMINATION ROUTINES                   |
+---------------------------------------------------------------------*/

/* This is the version number for the gwnum libraries. It changes whenever */
/* there is a change to the gwnum code and will match the corresponding */
/* prime95 version.  Thus, you may see some strange jumps in version */
/* numbers.  This version number is also embedded in the assembly code and */
/* gwsetup verifies that the version numbers match.  This prevents bugs */
/* from accidentally linking in the wrong gwnum library. */

#define GWNUM_VERSION		"24.14"
#define GWNUM_MAJOR_VERSION	24
#define GWNUM_MINOR_VERSION	14

/* Prior to calling one of the gwsetup routines, you can have the library */
/* "play it safe" by reducing the maximum allowable bits per FFT data word. */
/* For example, the code normally tests a maximum of 22477 bits in a 1024 */
/* SSE2 FFT, or 21.95 bits per double.  If you set the safety margin to 0.5 */
/* then the code will only allow 21.45 bits per double, or a maximum of */
/* 21965 bits in a 1024 length FFT. */

#define gwset_safety_margin(m)	(GWSAFETY_MARGIN = m)

/* The gwsetup routines needs to know the maximum value that will be used */
/* in a call to gwsetmulbyconst.  By default this value is assumed to be 3, */
/* which is what you would use in a base-3 Fermat PRP test.  Gwsetup must */
/* switch to a generic modular reduction if k * mulbyconst or c * mulbyconst */
/* is too large.  Call this routine prior to calling gwsetup. */

#define gwsetmaxmulbyconst(c) (GWMAXMULBYCONST = c)

/* Error codes returned by the three gwsetup routines */

#define GWERROR_VERSION		1001	/* GWNUM.H and FFT code version */
					/* numbers do not match. */
#define GWERROR_TOO_LARGE	1002	/* Number too large for the FFTs. */
#define GWERROR_K_TOO_SMALL	1003	/* k < 1 is not supported */
#define GWERROR_K_TOO_LARGE	1004	/* k > 53 bits is not supported */
#define GWERROR_MALLOC		1005	/* Insufficient memory available */

/* There are three different setup routines.  The first, gwsetup, is for */
/* gwnum's primary use - support for fast operations modulo K*B^N+C. */
/* Smaller K and C values result in smaller FFT sizes and faster operations. */
/* Right now, if B<>2 defaults to the slower gwsetup_general_mod case. */
/* Only choose a specific FFT size if you know what you are doing!! */

int gwsetup (
	double	k,		/* K in K*B^N+C. Must be a positive integer. */
	unsigned long b,	/* B in K*B^N+C. Must be two. */
	unsigned long n,	/* N in K*B^N+C. Exponent to test. */
	signed long c,		/* C in K*B^N+C. Must be rel. prime to K. */
	unsigned long fftlen);	/* Zero or specific FFT size to use. */

/* This setup routine is for operations modulo an arbitrary binary number. */
/* This is three times slower than the special forms above. */
/* Only choose a specific FFT size if you know what you are doing!! */

int gwsetup_general_mod (
	unsigned long *array,	/* The modulus as an array of longs */
	unsigned long arraylen,	/* The number of longs in the array */
	unsigned long fftlen);	/* Zero or specific FFT size to use. */

/* This setup routine is for operations without a modulo. In essence, */
/* you are using gwnums as a general-purpose FFT multiply library. */
/* Only choose a specific FFT size if you know what you are doing!! */

int gwsetup_without_mod (
	unsigned long n,	/* Maximum number of bits in OUTPUT numbers. */
	unsigned long fftlen);	/* Zero or specific FFT size to use. */

/* Free all memory allocated by gwnum routines since gwsetup was called. */

void gwdone (void);

/*---------------------------------------------------------------------+
|                     GWNUM MEMORY ALLOCATION ROUTINES                 |
+---------------------------------------------------------------------*/

/* Allocate memory for a gwnum */
gwnum gwalloc (void);

/* Free a previously allocated gwnum */
void gwfree (gwnum);

/* Free all previously allocated gwnums */
void gwfreeall (void);

/*---------------------------------------------------------------------+
|                        GWNUM CONVERSION ROUTINES                     |
+---------------------------------------------------------------------*/

/* Convert a double (must be an integer) to a gwnum */
void dbltogw (double, gwnum);

/* Convert a binary value to a gwnum */
void binarytogw (
	unsigned long *array,	/* Array containing the binary value */
	unsigned long arraylen,	/* Length of the array */
	gwnum	n);		/* Destination gwnum */

/* Convert a gwnum to a binary value.  Returns the number of longs */
/* written to the array.  The array is NOT zero-padded.  Returns a */
/* negative number if an error occurs during the conversion.  An error */
/* can happen if the FFT data contains a NaN or infinity value. */
long gwtobinary (
	gwnum	n,		/* Source gwnum */
	unsigned long *array,	/* Array to contain the binary value */
	unsigned long arraylen);/* Maximum size of the array */

/*---------------------------------------------------------------------+
|                          GWNUM MATH OPERATIONS                       |
+---------------------------------------------------------------------*/

/* Macros to interface with assembly code */
/* The assembly routines are designed to provide a flexible way of */
/* multiplying two numbers.  If you will use a value in several multiplies */
/* you can perform the forward transform just once.  Furthermore, the */
/* multiply routines are tuned to allow one unnormalized addition prior */
/* to a multiply without introducing too much convolution error.  Thus: */
/* Legal:	gwaddquick (t1, t2); gwmul (t2, x); */
/* Legal:	gwfft (t1, t1); gwfft (t2, t2); */
/*		gwfftadd (t1, t2); gwfftmul (t2, x); */
/* Not Legal:	gwaddquick (t1, t2); gwaddquick (y, x); gwmul (t2, x); */
/* Not Legal:	gwfft (t1, t1); gwfft (t2, t2); */
/*		gwfftadd (t1, t2); gwfftfftmul (t2, t2); */

/* A brief description of each of the "gw" routines: */
/* gwswap	Quickly swaps two gw numbers */
/* gwcopy(s,d)	Copies gwnum s to d */
/* gwadd	Adds two numbers and normalizes result if necessary */
/* gwsub	Subtracts first number from second number and normalizes
/*		result if necessary */
/* gwadd3quick	Adds two numbers WITHOUT normalizing */
/* gwsub3quick	Subtracts second number from first WITHOUT normalizing */
/* gwadd3	Adds two numbers and normalizes them if necessary */
/* gwsub3	Subtracts second number from first number and normalizes
/*		result if necessary */
/* gwaddsub	Adds and subtracts 2 numbers (first+second and first-second) */
/*		normalizes the results if necessary */
/* gwaddsub4	Like, gwaddsub but can store results in separate variables */
/* gwaddsub4quick Like, gwaddsub4 but will not do a normalize */
/* gwfft	Perform the forward Fourier transform on a number */
/* gwsquare	Multiplies a number by itself */
/* gwsquare_carefully  Like gwsquare but uses a slower method that will */
/*		have a low roundoff error even if input is non-random data */
/* gwmul(s,d)	Computes d=s*d.  NOTE: s is replaced by its FFT */
/* gwsafemul	Like gwmul but s is not replaced with its FFT */
/* gwfftmul(s,d) Computes d=s*d.  NOTE: s must have been previously FFTed */
/* gwfftfftmul(s1,s2,d) Computes d=s1*s2.  Both s1 and s2 must have */
/*		been previously FFTed */

/* The routines below operate on numbers that have already been FFTed. */

/* gwfftadd	Adds two FFTed numbers */
/* gwfftsub	Subtracts first FFTed number from second FFTed number */
/* gwfftadd3	Adds two FFTed numbers */
/* gwfftsub3	Subtracts second FFTed number from first FFTed number */
/* gwfftaddsub	Adds and subtracts 2 FFTed numbers */
/* gwfftaddsub4	Like, gwfftaddsub but stores results in separate variables */

#define gwswap(s,d)	{gwnum t; t = s; s = d; d = t;}
#define gwaddquick(s,d)	gwadd3quick (s,d,d)
#define gwsubquick(s,d)	gwsub3quick (d,s,d)
#define gwadd(s,d)	gwadd3 (s,d,d)
#define gwsub(s,d)	gwsub3 (d,s,d)
#define gwaddsub(a,b)	gwaddsub4 (a,b,a,b)
#define gwaddsubquick(a,b) gwaddsub4quick (a,b,a,b)
#define gwtouch(s)	gwcopy (s,s)
#define gwfftadd(s,d)	gwfftadd3 (s,d,d)
#define gwfftsub(s,d)	gwfftsub3 (d,s,d)
#define gwfftaddsub(a,b) gwfftaddsub4 (a,b,a,b)

/* Set the constant which the results of a multiplication should be */
/* multiplied by.  Use this macro in conjunction with the c argument of */
/* gwsetnormroutine. */

void gwsetmulbyconst (long s);

/* The multiplication code has two options that you can set using this */
/* macro.  The e argument tells the multiplication code whether or not */
/* it should perform round-off error checking - returning the maximum */
/* difference from an integer result in MAXERR.  The c argument tells the */
/* multiplication code whether or not it should multiply the result by */
/* a small constant. */

#define gwsetnormroutine(z,e,c) {NORMNUM=2*(c)+(e);}

/* If you know the result of a multiplication will be the input to another */
/* multiplication (but not gwsquare_carefully), then a small performance */
/* gain can be had in larger FFTs by doing some of the next forward FFT at */
/* the end of the multiplication.  Call this macro to tell the */
/* multiplication code whether or not it can start the forward FFT on */
/* the result. */

#define gwstartnextfft(f) {if (!GENERAL_MOD) POSTFFT=f;}

void gwcopy (			/* Copy a gwnum */
	gwnum	s,		/* Source */
	gwnum	d);		/* Dest */
void gwfft (			/* Forward FFT */
	gwnum	s,		/* Source number */
	gwnum	d);		/* Destination (can overlap source) */
void gwsquare (			/* Square a number */
	gwnum	s);		/* Source and destination */
void gwmul (			/* Multiply source with dest */
	gwnum	s,		/* Source number (changed to FFTed source!) */
	gwnum	d);		/* Source and destination */
void gwsafemul (		/* Multiply source with dest */
	gwnum	s,		/* Source number (not changed) */
	gwnum	d);		/* Source and destination */
void gwfftmul (			/* Multiply already FFTed source with dest */
	gwnum	s,		/* Already FFTed source number */
	gwnum	d);		/* Non-FFTed source. Also destination */
void gwfftfftmul (		/* Multiply two already FFTed sources */
	gwnum	s,		/* Already FFTed source number */
	gwnum	s2,		/* Already FFTed source number */
	gwnum	d);		/* Destination (can overlap sources) */
void gwadd3quick (		/* Add two numbers without normalizing */
	gwnum	s1,		/* Source #1 */
	gwnum	s2,		/* Source #2 */
	gwnum	d);		/* Destination */
void gwsub3quick (		/* Compute s1 - s2 without normalizing */
	gwnum	s1,		/* Source #1 */
	gwnum	s2,		/* Source #2 */
	gwnum	d);		/* Destination */
void gwaddsub4quick (		/* Add & sub two numbers without normalizing */
	gwnum	s1,		/* Source #1 */
	gwnum	s2,		/* Source #2 */
	gwnum	d1,		/* Destination #1 */
	gwnum	d2);		/* Destination #2 */
void gwadd3 (			/* Add two numbers normalizing if needed */
	gwnum	s1,		/* Source #1 */
	gwnum	s2,		/* Source #2 */
	gwnum	d);		/* Destination */
void gwsub3 (			/* Compute s1 - s2 normalizing if needed */
	gwnum	s1,		/* Source #1 */
	gwnum	s2,		/* Source #2 */
	gwnum	d);		/* Destination */
void gwaddsub4 (		/* Add & sub two nums normalizing if needed */
	gwnum	s1,		/* Source #1 */
	gwnum	s2,		/* Source #2 */
	gwnum	d1,		/* Destination #1 */
	gwnum	d2);		/* Destination #2 */
void gwfftadd3 (		/* Add two FFTed numbers */
	gwnum	s1,		/* Source #1 */
	gwnum	s2,		/* Source #2 */
	gwnum	d);		/* Destination */
void gwfftsub3 (		/* Compute FFTed s1 - FFTed s2 */
	gwnum	s1,		/* Source #1 */
	gwnum	s2,		/* Source #2 */
	gwnum	d);		/* Destination */
void gwfftaddsub4 (		/* Add & sub two FFTed numbers */
	gwnum	s1,		/* Source #1 */
	gwnum	s2,		/* Source #2 */
	gwnum	d1,		/* Destination #1 */
	gwnum	d2);		/* Destination #2 */

/* Square a number using a slower method that will have reduced */
/* round-off error on non-random input data.  Caller must make sure the */
/* input number has not been partially (via gwstartnextfft) or fully FFTed. */

void gwsquare_carefully (
	gwnum	s);		/* Source and destination */

/* These routines can be used to add a constant to the result of a */
/* multiplication.  Using these routines lets prime95 do the -2 operation */
/* in a Lucas-Lehmer test and use the gwstartnextfft macro for a small */
/* speedup.  NOTE:  There are some number formats that cannot use this */
/* routine.  If abs(c) in k*b^n+c is 1, then gwsetaddin can be used. */
/* To use gwsetaddinatbit, k must also be 1. */

void gwsetaddin (long);
void gwsetaddinatbit (long, unsigned long);

/* This routine adds a small value to a gwnum.  This lets us apply some */
/* optimizations that cannot be performed by gwadd */

void gwaddsmall (gwnum g, int addin);

/*---------------------------------------------------------------------+
|                      GWNUM ERROR-CHECKING ROUTINES                   |
+---------------------------------------------------------------------*/

#define gw_test_for_error()		GWERROR
#define gw_test_illegal_sumout()	(GWERROR & 1)
#define gw_test_mismatched_sums()	(GWERROR & 2)
#define gwsuminp(g)			((g)[-2])
#define gwsumout(g)			((g)[-3])

/* Return TRUE if we are operating near the limit of this FFT length */
/* Input argument is the percentage to consider as near the limit. */
/* For example, if percent is 1.0 and the FFT can handle 20 bits per FFT */
/* data word, then if there are more than 19.98 bits per FFT data word */
/* this function will return TRUE. */

int gwnear_fft_limit (double pct);

/*---------------------------------------------------------------------+
|                    GWNUM MISC. INFORMATION ROUTINES                  |
+---------------------------------------------------------------------*/

/* Return TRUE if this is a GPL'ed version of the GWNUM library. */
#define gwnum_is_gpl()		(0)

/* Generate a human-readable description of the FFT size chosen */
void gwfft_description (char *buf);

/* Generate a human-readable string for k*b^n+c */
void gw_as_string(char *buf, double k, unsigned long b, unsigned long n,
		  signed long c);

/* A human-readable string for the modulus currently in use */
#define gwmodulo_as_string()	(GWSTRING_REP)

/*---------------------------------------------------------------------+
|                 ALTERNATIVE INTERFACES USING GIANTS                  |
+---------------------------------------------------------------------*/

/* The giants library from Dr. Richard Crandall, Perfectly Scientific, */
/* is used internally for a few infrequent operations.  It can optionally */
/* be used in the interfaces to convert between gwnum data type and binary. */
/* I do not recommend this.  There are many other faster and more robust */
/* libraries available. */

#include "giants.h"

/* Same as gwsetup_general_mod but uses giants instead of array of longs */
int gwsetup_general_mod_giant (
	giant n,		/* The modulus */
	unsigned long fftlen);	/* Zero or specific FFT size to use. */

/* Convert a giant to a gwnum */
void gianttogw (giant, gwnum);

/* Convert a gwnum to a giant.  WARNING: Caller must allocate an array that */
/* is several words larger than the maximum result that can be returned. */
/* This is a gross kludge that lets gwtogiant use the giant for intermediate */
/* calculations.  Returns a negative number if an error occurs.  Returns */
/* zero on success. */
int gwtogiant (gwnum, giant);

/*---------------------------------------------------------------------+
|          MISC. CONSTANTS YOU PROBABLY SHOULDN'T CARE ABOUT           |
+---------------------------------------------------------------------*/

/* The maximum value k * mulbyconst can be in a zero pad FFT.  Larger */
/* values must use generic modular reduction. */

#define MAX_ZEROPAD_K	2251799813685247.0	/* 51-bit k's are OK. */

/* The maximum value c * mulbyconst can be in a zero pad FFT.  Larger */
/* values must use generic modular reduction. */

#define MAX_ZEROPAD_C	8388607			/* 23-bit c's seem to work. */

/*---------------------------------------------------------------------+
|          SPECIAL ECM ROUTINE FOR GMP-ECM USING GWNUM LIBRARY         |
+---------------------------------------------------------------------*/

/* Return codes */

#define ES1_SUCCESS		0	/* Success, but no factor */
#define ES1_FACTOR_FOUND	1	/* Success, factor found */
#define ES1_CANNOT_DO_IT	2	/* This k,b,n,c cannot be handled */
#define ES1_MEMORY		3	/* Out of memory */
#define ES1_INTERRUPT		4	/* Execution interrupted */
#define ES1_CANNOT_DO_QUICKLY	5	/* Requires 3-multiply reduction */
#define ES1_HARDWARE_ERROR	6	/* An error was detected, most */
					/* likely a hardware error. */

/* Option codes */

#define ES1_DO_SLOW_CASE	0x1	/* Set this if ecmStage1 should do */
					/* slow 3-multiply reduction cases. */

/* INPUTS:

Input number (3 possibilities):

1) k,b,n,c set and num_being_factored_array = NULL.  k*b^n+c is factored.
2) k,b,n,c zero and num_being_factored_array set.  num_being_factored is
   worked on using generic 3-multiply reduction
3) k,b,n,c set and num_being_factored_array set.  num_being_factored is
   worked on - it must be a factor of k*b^n+c.

A_array, B1 are required

B1_done is optional.  Use it to resume a stage 1 calculation.

x_array, z_array is the starting point.  If z_array is not given, then
z is assumed to be one.

stop_check_proc is optional

options are defined above


   OUTPUTS:

On success:

   if z_array is NULL, then x_array is set to normalized point
   else x_array, z_array is set to the unnormalized point

On factor found:

   x_array is set to the factor found

On interrupt:

   B1_done is set to the last prime below B1 that was processed.
   If z_array != NULL (preferred) then x_array and z_array are set to the
current point.  The (x,z) point is not normalized because it will
be slow for large numbers.  This is unacceptable during system shutdown.
Caller must allocate x and z arrays large enough to hold any k*b^n+c value.
   If z_array == NULL, then a normalized x is returned. Caller must allocate
x array large enough to hold any value less than num_being_factored.

*/

int gwnum_ecmStage1 (
	double	k,			/* K in K*B^N+C */
	unsigned long b,		/* B in K*B^N+C */
	unsigned long n,		/* N in K*B^N+C */
	signed long c,			/* C in K*B^N+C */
	unsigned long *num_being_factored_array, /* Number to factor */
	unsigned long num_being_factored_array_len,
	unsigned long B1,		/* Stage 1 bound */
	unsigned long *B1_done,		/* Stage 1 that is already done */
	unsigned long *A_array,		/* A - caller derives it from sigma */
	unsigned long A_array_len,
	unsigned long *x_array,		/* X value of point */
	unsigned long *x_array_len,
	unsigned long *z_array,		/* Z value of point */
	unsigned long *z_array_len,
	int	(*stop_check_proc)(void),/* Ptr to proc that returns TRUE */
					/* if user interrupts processing */
	unsigned long options);

/*---------------------------------------------------------------------+
|                             GWNUM INTERNALS                          |
+---------------------------------------------------------------------*/

/* A psuedo declaration for our big numbers.  The actual pointers to */
/* these big numbers are to the data array.  The 96 bytes prior to the */
/* data contain: */
/* data-4:  integer containing number of unnormalized adds that have been */
/*	    done.  After a certain number of unnormalized adds, the next add */
/*	    must be normalized to avoid overflow errors during a multiply. */
/* data-8:  integer containing number of bytes in data area. Used by gwcopy. */
/* data-16: double containing the product of the two sums of the input FFT */
/*	    values. */
/* data-24: double containing the sum of the output FFT values.  These two */
/*	    values can be used as a sanity check when multiplying numbers. */
/*	    The two values should be "reasonably close" to one another. */
/* data-28: Flag indicating gwnum value has been partially FFTed. */
/* data-32: Pointer returned by malloc - used to free memory when done. */
/* data-88: Seven doubles (input FFT values near the halfway point */
/*	    when doing a zero-padded FFT). */
/* data-96: Eight unused bytes */
/* typedef struct {
/*	char	pad[96];	   Used to track unnormalized add/sub */
/*				   and original address */
/*	double	data[512];	   The big number broken into chunks */
/*				   This array is variably sized. */
/* } *gwnum; */
#define GW_HEADER_SIZE	96	/* Number of data bytes before a gwnum ptr */

/* Some mis-named #defines that describe the maximum Mersenne number */
/* exponent that the gwnum routines can process. */

#define MAX_PRIME	79300000L	/* Maximum number of x87 bits */
#define MAX_PRIME_SSE2	596000000L	/* SSE2 bit limit */
#define MAX_FFTLEN	4194304L	/* 4M FFT max for x87 */
#define MAX_FFTLEN_SSE2	33554432L	/* 32M FFT max for SSE2 */

/* global variables */

extern double GWSAFETY_MARGIN;	/* Reduce maximum allowable bits per */
				/* FFT data word by this amount. */
extern long GWMAXMULBYCONST;	/* Gwsetup needs to know the maximum value */
				/* the caller will use in gwsetmulbyconst. */
				/* The default value is 3, commonly used */
				/* in a base-3 Fermat PRP test. */
extern double KARG;		/* K in K*B^N+C */
extern unsigned long BARG;	/* B in K*B^N+C */
extern unsigned long PARG;	/* N in K*B^N+C */
extern signed long CARG;	/* C in K*B^N+C */
extern unsigned long FFTLEN;	/* The FFT size we are using */
extern unsigned long RATIONAL_FFT;/* TRUE if bits per FFT word is integer */
extern unsigned long BITS_PER_WORD;/* Bits in a little word */
extern int GWERROR;		/* Set if an error is detected */
extern double MAXERR;		/* Convolution error in a multiplication */
extern double MAXDIFF;		/* Maximum allowable difference between */
				/* sum of inputs and outputs */
extern void (*GWPROCPTRS[24])(void);/* Ptrs to assembly routines */
extern unsigned int NORMNUM;	/* The post-multiply normalization routine */
extern void *SRCARG;		/* For assembly language arg passing */
extern void *SRC2ARG;		/* For assembly language arg passing */
extern void *DESTARG;		/* For assembly language arg passing */
extern void *DEST2ARG;		/* For assembly language arg passing */
extern long NUMARG;		/* For assembly language arg passing */
extern long NUM2ARG;		/* For assembly language arg passing */
extern double fft_count;	/* Count of forward and inverse FFTs */
extern char *GW_BIGBUF;		/* Optional buffer to allocate gwnums in */
extern unsigned long GW_BIGBUF_SIZE;	/* Size of the optional buffer */
extern gwnum *gwnum_alloc;		/* Array of allocated gwnums */
extern unsigned int gwnum_alloc_count;	/* Count of allocated gwnums */
extern unsigned int gwnum_alloc_array_size;/* Size of gwnum_alloc array */
extern gwnum *gwnum_free;		/* Array of available gwnums */
extern unsigned int gwnum_free_count;	/* Count of available gwnums */
extern unsigned long GW_ALIGNMENT;	/* How to align allocated gwnums */
extern int GENERAL_MOD;		/* True if doing general-purpose mod */
extern char GWSTRING_REP[40];	/* The gwsetup modulo number as a string. */
extern unsigned long ASM_TIMERS[32];/* Assembly language timers */
extern unsigned long POSTFFT;	/* True if assembly code can start the */
				/* FFT process on the result of a multiply */
extern unsigned long PASS1_CACHE_LINES; /* Cache lines grouped together in */
				/* first pass of an FFT. */
extern unsigned long PASS2_LEVELS;	/* FFT levels done in pass 2. */

double *addr (gwnum, unsigned long);
unsigned long gwnum_size (unsigned long);
int get_fft_value (gwnum, unsigned long, long *);
void set_fft_value (gwnum, unsigned long, long);
int is_big_word (unsigned long);
void bitaddr (unsigned long, unsigned long *, unsigned long *);
void specialmodg (giant);
#define gw_set_max_allocs(n)	if (gwnum_alloc==NULL) gwnum_alloc_array_size=n

unsigned long gwmap_to_fftlen (double, unsigned long, unsigned long, signed long);
double gwmap_to_timing (double, unsigned long, unsigned long, signed long, int);
unsigned long gwmap_to_memused (double, unsigned long, unsigned long, signed long);
unsigned long map_fftlen_to_max_exponent (unsigned long);

/* Speed of other processors compared to a Pentium II of same clock speed */

#define REL_486_SPEED	8.4	/* 486 is over 8 times slower than PII */
#define REL_K6_SPEED	3.0	/* K6 is 3 times slower than PII */
#define REL_PENT_SPEED	1.2	/* Pentium is 20% slower than PII */
#define REL_K7_SPEED	0.7	/* Assume K7 is faster than a PII */
#define REL_P4_SPEED	0.7	/* Assume P4 is faster than a PII factoring*/

/* Specialized routines that let the internal giants code share the free */
/* memory pool used by gwnums. */

void gwfree_temporarily (gwnum);
void gwrealloc_temporarily (gwnum);

/* Other routines used internally */

unsigned long addr_offset (unsigned long, unsigned long);
extern void eset_mul_const (void);

/* When debugging gwnum and giants, I sometimes write code that "cheats" */
/* by calling a routine that is part of prime95 rather than the gwnum and */
/* giants library.  Prime95 will set this routine pointer so that gwnum */
/* code can cheat while keeping the gwnum library interface clean. */

extern void (*OutputBothRoutine)(char *);

/* This routine lets me time many assembly language building blocks -- used */
/* when optimizing these building blocks. */

int timeit (void);

/* This global forces the FFT selection code to pick the n-th possible */
/* implementation instead of the best one.  The prime95 benchmarking code */
/* uses this to time every FFT implementation so that we can find the best */
/* one for a new CPU architecture.  DO NOT set this variable. */

extern	int bench_pick_nth_fft;

#ifdef __cplusplus
}
#endif

#endif
