/**************************************************************
 *
 *	giants.h
 *
 *	Header file for large-integer arithmetic library giants.c.
 *
 *	Updates:
 *          30 Apr 98  JF   USE_ASSEMBLER_MUL removed
 *          29 Apr 98  JF   Function prototypes cleaned up
 *	    20 Apr 97  RDW
 *
 *	c. 1997 Perfectly Scientific, Inc.
 *	c. 1998-2005 Just For Fun Software, Inc.
 *	All Rights Reserved.
 *
 **************************************************************/

#ifndef _GIANTS_H_
#define _GIANTS_H_

/* This is a C library.  If used in a C++ program, don't let the C++ */
/* compiler mangle names. */

#ifdef __cplusplus
extern "C" {
#endif

/* Include common definitions */

#include "common.h"

/**************************************************************
 *
 * Mul options
 *
 **************************************************************/

#define AUTO_MUL 0
#define GRAMMAR_MUL 1
#define FFT_MUL 2
#define KARAT_MUL 3

/**************************************************************
 *
 * Structure definitions
 *
 **************************************************************/

typedef struct
{
	 int	sign;
	 unsigned long *n;		/* ptr to array of longs */
#ifdef GDEBUG
	int	maxsize;
#endif
} giantstruct;

typedef giantstruct *giant;

#ifdef GDEBUG
#define setmaxsize(g,s)	(g)->maxsize = s
#else
#define setmaxsize(g,s)
#endif

/**************************************************************
 *
 * Function Prototypes
 *
 **************************************************************/

/* Cleanup giants allocations, etc. */
void	term_giants (void);

/* Creates a new giant */
giant 	newgiant(int numshorts);

/* Returns the bit-length n; e.g. n=7 returns 3. */
int 	bitlen(giant n);

/* Returns the value of the pos bit of g. */
#define bitval(g,pos)	((g)->n[(pos)>>5] & (1 << ((pos)&31)))

/* Returns whether g is zero. */
#define isZero(g)	((g)->sign == 0)

/* Returns whether g is one. */
#define isone(g)	(((g)->sign == 1) && ((g)->n[0] == 1))

/* Returns whether g is one. */
#define istwo(g)	(((g)->sign == 1) && ((g)->n[0] == 2))

/* Copies one giant to another. */
void 	gtog(giant src, giant dest);

/* Integer -> giant. */
#define setzero(g)	(g)->sign = 0
#define setone(g)	(g)->sign = (g)->n[0] = 1
void 	itog(int n, giant g);
void 	ultog(unsigned long n, giant g);

/* Double -> giant. */
void 	dbltog(double n, giant g);

/* Character string <-> giant. */
void 	ctog (char *s, giant g);
void 	gtoc (giant g, char *s, int sbufsize);

/* Returns the sign of g: -1, 0, 1. */
int	gsign(giant g);

/* Returns 1, 0, -1 as a>b, a=b, a<b. */
int	gcompg(giant a, giant b);

/* Set AUTO_MUL for automatic FFT crossover (this is the
 * default), set FFT_MUL for forced FFT multiply, set
 * GRAMMAR_MUL for forced grammar school multiply. */
void	setmulmode(int mode);

/**************************************************************
 *
 * Math Functions
 *
 **************************************************************/

/* g := -g. */
#define negg(g)		((g)->sign = -(g)->sign)

/* g := |g|. */
#define absg(g)		if ((g)->sign < 0) (g)->sign = -(g)->sign

/* g += i, with i an int (for compatability with old giants.h */
#define	iaddg(i,g)	sladdg ((long)(i), g)

/* g += i, where i is an unsigned long */
void 	uladdg (unsigned long i, giant g);

/* g += i, where i is a signed long */
void 	sladdg (long i, giant g);

/* b += a. */
void 	addg(giant a, giant b);

/* g -= i, where i is an unsigned long */
void 	ulsubg (unsigned long i, giant g);

/* b -= a. */
void 	subg(giant a, giant b);

/* g *= g. */
void	squareg(giant g);

/* b *= a. */
void	mulg(giant a, giant b);

/* g *= i, where i is an unsigned long */
void 	ulmulg (unsigned long i, giant g);
void 	imulg (long i, giant g);

/* g *= n, where n is a double */
void 	dblmulg (double n, giant g);

/* num := num % den, any positive den. */
void 	modg(giant den, giant num);

/* num := [num/den], any positive den. */
void 	divg(giant den, giant num);
void 	dbldivg(double den, giant num);

/* Mask rightmost bits in g, same as g = g % 2^bits. */
void 	gmaskbits(int bits, giant g);

/* Shift g left by bits, introducing zeros on the right. */
void 	gshiftleft(int bits, giant g);

/* Shift g right by bits, losing bits on the right. */
#define	gshiftright(n,g)	{if (n) gtogshiftright (n, g, g);}
void 	gtogshiftright (int bits, giant src, giant dest);

/* If 1/x exists (mod n), then x := 1/x.  If
 * inverse does not exist, then x := - GCD(n, x). */
void 	invg(giant n, giant x);
int 	invgi(giant n, giant x);	/* Interruptable version */

/* General GCD, x:= GCD(n, x). */
void 	gcdg(giant n, giant x);
int 	gcdgi(giant n, giant x);	/* Interruptable version */

/* x becomes x^n, NO mod performed. */
void power (giant x, int n);
void powerg (giant x, giant n);

/* num := [num/den], any positive den. */
void 	powermod(giant x, int n, giant z);

/* x := x^n (mod z). */
void 	powermodg(giant x, giant n, giant z);

/* r becomes the steady-state reciprocal 2^(2b)/d, where
 * b = bit-length of d-1. */
void	make_recip(giant d, giant r);

/* n := [n/d], d positive, using stored reciprocal directly. */
void	divg_via_recip(giant d, giant r, giant n);

/* stack handling functions */
giant	popg(int);	/* Number of longs in data area */
void	pushg(int);	/* Number of items to return to stack */



/* Handle the difference between the naming conventions in */
/* C compilers.  We need to do this for global variables that */
/* referenced by the assembly routines.  Most non-Windows systems */
/* should #define ADD_UNDERSCORES before including this file. */

#ifdef ADD_UNDERSCORES
#include "gwrename.h"
#endif

/* Low-level math routines the caller can use for multi-precision */
/* arithmetic */

extern unsigned long CARRYH;	/* For multi-precision asm routines */
extern unsigned long CARRYL;
extern unsigned long RES;

extern void eaddhlp (void);
extern void esubhlp (void);
extern void emuladdhlp (void);
extern void emuladd2hlp (void);
extern void emulsubhlp (void);
#define addhlp(a)	NUMARG=(long)(a), eaddhlp()
#define subhlp(a)	NUMARG=(long)(a), esubhlp()
#define muladdhlp(a,b)	{NUMARG=(long)(a); NUM2ARG=(long)(b); emuladdhlp();}
#define muladd2hlp(a,b)	{NUMARG=(long)(a); NUM2ARG=(long)(b); emuladd2hlp();}
#define mulsubhlp(a,b)	{NUMARG=(long)(a); NUM2ARG=(long)(b); emulsubhlp();}

/* External routine pointers. */

extern int (*StopCheckRoutine)(void);

#ifdef __cplusplus
}
#endif

#endif
