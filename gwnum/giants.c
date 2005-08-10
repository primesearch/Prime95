/**************************************************************
 *
 *  giants.c
 *
 *  Library for large-integer arithmetic.
 * 
 *  Massive rewrite by G. Woltman for 32-bit support
 *
 *  c. 1997,1998 Perfectly Scientific, Inc.
 *  c. 1998-2005 Just For Fun Software, Inc.
 *  All Rights Reserved.
 *
 **************************************************************/

/* Include Files */

#include <stdio.h>
#include <stdlib.h>
#include <ctype.h>
#include <string.h>
#include <math.h>
#include "giants.h"
#include "gwnum.h"
#include "fftsg.c"
#include "gwutil.h"

/**************************************************************
 *
 * Preprocessor definitions
 *
 **************************************************************/

#define TWOPI (double)(2*3.1415926535897932384626433)
#define SQRT2 (double)(1.414213562373095048801688724209)
#define SQRTHALF (double)(0.707106781186547524400844362104)
#define TWO_TO_MINUS_19 (double)(0.0000019073486328125)
#define ONE_MINUS_TWO_TO_MINUS_19 (double)(0.9999980926513671875)

/* Next, number of words at which Karatsuba becomes better than GRAMMAR. */
#define KARAT_BREAK_SQUARE	48
#define KARAT_BREAK_MULT	43

/* Next, mumber of words at which FFT becomes better than Karatsuba. */
#define FFT_BREAK_SQUARE	96
#define FFT_BREAK_MULT1		47
#define FFT_BREAK_MULT2		64
#define FFT_BREAK_MULT3		74

/* The limit (in 32-bit words) below which hgcd is too ponderous */
#define GCDLIMIT 150

/* Size by which to increment the stack used in pushg() and popg(). */
#define	STACK_GROW	100

/* Size at which a brute force hgcd is done */
#define CHGCD_BREAK	24

typedef struct gmatrixstruct {
	 giant	ul;			/* upper left */
	 giant	ur;			/* upper right */
	 giant	ll;			/* lower left */
	 giant	lr;			/* lower right */
} gmatrixstruct;
typedef gmatrixstruct *gmatrix;

/* Routines to help manage pushing the right number of temporaries */
/* back on the stack.  Especially handy for routines where user can */
/* hit escape in the middle of the computation. */

#define stackstart()	cur_stack_elem
#define pushall(x)	pushg(cur_stack_elem - x)

typedef struct gstacknodestruct {
	 gwnum	gw;			/* gwnum used */
	 unsigned long space_left;
	 long	offset;
	 long	next_offset;
} gstacknode;


/* Global variables. */

int	mulmode = AUTO_MUL;
int	cur_stack_size = 0;
int	cur_stack_elem = 0;
gstacknode *stack = NULL;
double	*ooura_sincos = NULL;
int	*ooura_ip = NULL;
int	ooura_fft_size = 0;
giant	cur_recip = NULL;
giant	cur_den = NULL;

/* Private function prototypes. */

void		normal_addg(giant, giant);
void		normal_subg(giant, giant);
void		reverse_subg(giant, giant, int);
void 		automulg(giant a, giant b);
void 		grammarmulg(giant a, giant b);
void		grammarsquareg(giant b);
void 		karatmulg(giant a, giant b);
void 		karatsquareg(giant b);

void		init_sincos(int);
int 		lpt(int);
void 		addsignal(giant, int, double *, int);
void 		FFTsquareg(giant x);
void 		FFTmulg(giant y, giant x);
void 		giant_to_double(giant x, int sizex, double *z, int L);

#define gswap(p,q)  {giant tgq; tgq = *(p); *(p) = *(q); *(q) = tgq;}
#define ulswap(p,q)  {unsigned long tgq; tgq = p; p = q; q = tgq;}

int		invg_common (giant, giant, int);
int		gcdg_common (giant, giant, int);
int	 	cextgcdg (giant *, giant *, gmatrix A, int);
int		ggcd (giant *, giant *, gmatrix, int);
void 		onestep (giant *, giant *, gmatrix);
void 		mulvM (gmatrix, giant, giant);
void 		mulmM (gmatrix, gmatrix);
void 		mulmMsp (gmatrix, gmatrix, int);
void		punch (giant, gmatrix);
int		hgcd (int, giant *, giant *, gmatrix, int);
int		rhgcd (giant *, giant *, gmatrix, int);

unsigned long EGCD_A=0;	/* Returned by assembly GCD helper */
unsigned long EGCD_B=0;
unsigned long EGCD_C=0;
unsigned long EGCD_D=0;
unsigned long EGCD_ODD=0;

unsigned long CARRYH=0;	/* For multi-precision asm routines */
unsigned long CARRYL=0;
unsigned long RES=0;

/* Internal routines */

int egcdhlp (void);

/* External routines */

int (*StopCheckRoutine)(void) = NULL;
#define stopCheck()	(StopCheckRoutine != NULL && (*StopCheckRoutine)())

/**************************************************************
 *
 *	Functions
 *
 **************************************************************/

void term_giants (void)
{
	pushg (cur_stack_elem);
	free (stack);
	stack = NULL;
	cur_stack_elem = 0;
	cur_stack_size = 0;
	aligned_free (ooura_sincos);
	ooura_sincos = NULL;
	free (ooura_ip);
	ooura_ip = NULL;
	ooura_fft_size = 0;
	free (cur_recip);
	free (cur_den);
	cur_recip = NULL;
	cur_den = NULL;
}

giant newgiant (		/* Create a new giant */
	int 	numshorts)
{
	int 	numlongs, size;
	giant 	thegiant;

	ASSERTG (numshorts > 0);

	numlongs = (numshorts + 1) / 2;
	size = sizeof (giantstruct) + numlongs * sizeof (unsigned long);
	thegiant = (giant) malloc (size);
	thegiant->sign = 0;
	thegiant->n = (unsigned long *)
		((char *) thegiant + sizeof (giantstruct));
	setmaxsize (thegiant, numlongs);
	return (thegiant);
}

void itog (		/* The giant g becomes set to the integer value i. */
	int	i,
	giant	g)
{
	if (i > 0) {
		g->sign = 1;
		g->n[0] = i;
	} else if (i == 0) {
		g->sign = 0;
		g->n[0] = 0;
	} else {
		g->sign = -1;
		g->n[0] = -i;
	}
}

void ultog (		/* The giant g becomes set to the integer value i. */
	unsigned long i,
	giant	g)
{
	if (i == 0) {
		g->sign = 0;
		g->n[0] = 0;
	} else {
		g->sign = 1;
		g->n[0] = i;
	}
}

void dbltog (		/* The giant g becomes set to the double i. */
	double	i,
	giant	g)
{
	if (i < 0.0) {
		dbltog (-i, g);
		negg (g);
	} else if (i < 4294967296.0) {
		ultog ((unsigned long) i, g);
	} else {
		g->sign = 2;
		g->n[1] = (unsigned long) (i / 4294967296.0);
		g->n[0] = (unsigned long) (i - g->n[1] * 4294967296.0);
	}
}

void ctog (		/* The giant g is set to string s. */
	char	*s,
	giant	g)
{
	for (g->sign = 0; isdigit (*s); s++) {
		ulmulg (10, g);
		uladdg (*s - '0', g);
	}
}

void gtoc (		/* The giant g is converted to string s. */
	giant	g,
	char	*s,
	int	sbufsize)
{
	giant	x, y, ten;
	int	i, len;
	char	c;

	ASSERTG (g->sign >= 0);

	x = popg (g->sign); gtog (g, x);
	y = popg (g->sign);
	ten = popg (1); itog (10, ten);
	sbufsize--;
	for (len = 0; len < sbufsize && x->sign; len++) {
		gtog (x, y);
		modg (ten, y);
		s[len] = (char) (y->n[0] + '0');
		divg (ten, x);
	}
	for (i = 0; i < len / 2; i++) {
		c = s[i];
		s[i] = s[len-1-i];
		s[len-1-i] = c;
	}
	s[len] = 0;
	pushg (3);
}

void gtog (			/* destgiant becomes equal to srcgiant. */
	giant	srcgiant,
	giant	destgiant)
{
	destgiant->sign = srcgiant->sign;
	memcpy ((char *) destgiant->n,
		(char *) srcgiant->n,
		abs (srcgiant->sign) * sizeof (unsigned long));
	ASSERTG (abs (destgiant->sign) <= destgiant->maxsize);
}

int gcompg (		/* Returns -1,0,1 if a<b, a=b, a>b, respectively. */
	giant	a,
	giant	b)
{
	int	sa = a->sign, j, sb = b->sign, sgn;
	unsigned long va, vb;

	if (sa > sb)
		return (1);
	if (sa < sb)
		return (-1);
	if (sa < 0) {
		sa = -sa; /* Take absolute value of sa. */
		sgn = -1;
	} else {
		sgn = 1;
	}
	for (j = sa-1; j >= 0; j--) {
		va = a->n[j];
		vb = b->n[j];
		if (va > vb) return (sgn);
		if (va < vb) return (-sgn);
	}
	return (0);
}

/* New add/subtract routines.
	The basic subtract "subg()" uses the following logic table:

     a      b          if (b > a)          if (a > b)
     
     +      +          b := b - a          b := -(a - b)
     -      +          b := b + (-a)       N.A.
     +      -          N.A.                b := -((-b) + a)
     -      -          b := (-a) - (-b)    b := -((-b) - (-a))

   The basic addition routine "addg()" uses:

     a      b          if(b > -a)          if(-a > b)
     
     +      +          b := b + a          N.A. 
     -      +          b := b - (-a)       b := -((-a) - b)
     +      -          b := a - (-b)       b := -((-b) - a)
     -      -          N.A.                b := -((-b) + (-a))

   In this way, internal routines "normal_addg," "normal_subg," 
	and "reverse_subg;" each of which assumes non-negative
	operands and a non-negative result, are now used for greater
	efficiency.
*/

void normal_addg (	/* b := a + b, both a,b assumed non-negative. */
	giant	a,
	giant	b)
{
	int 	j, asize = a->sign, bsize = b->sign;
	unsigned long *aptr = a->n, *bptr = b->n;

	ASSERTG (asize >= 0 && bsize >= 0);
	ASSERTG (a->sign == 0 || a->n[abs(a->sign)-1] != 0);
	ASSERTG (b->sign == 0 || b->n[abs(b->sign)-1] != 0);

	RES = CARRYL = 0;
	for (j = 0; j < asize || RES; j++) {
		if (j < asize) {
			addhlp (*aptr);
			aptr++;
		}
		if (j < bsize) addhlp (*bptr);
		*bptr++ = RES;
		RES = CARRYL;
		CARRYL = 0;
	}
	if (j > bsize) b->sign = j;

	ASSERTG (b->sign == 0 || b->n[abs(b->sign)-1] != 0);
}

void normal_subg (	/* b := b - a; assumes b, a positive and b >= a. */
	giant	a,
	giant	b)
{
	int 	j, asize = a->sign, bsize = b->sign;
	unsigned long *aptr = a->n, *bptr = b->n;

	ASSERTG (asize >= 0 && bsize >= asize);

	RES = CARRYL = CARRYH = 0;
	for (j = 0; j < asize; j++) {
		addhlp (*bptr);
		subhlp (*aptr);
		aptr++;
		*bptr++ = RES;
		RES = CARRYL;
		CARRYL = CARRYH;
	}
	for ( ; RES; j++)
	{
		if (j < bsize) addhlp (*bptr);
		*bptr++ = RES;
		RES = CARRYL;
		CARRYL = CARRYH;
	}
	if (j < bsize) return;

	while (j > 0 && b->n[j-1] == 0) j--;
	b->sign = j;

	ASSERTG (b->sign == 0 || b->n[abs(b->sign)-1] != 0);
}

void reverse_subg (	/* b := a - b; assumes b, a positive and a >= b. */
	giant	a,
	giant	b,
	int	azeros)
{
	int 	j, asize = a->sign + azeros, bsize = b->sign;
	unsigned long *aptr = a->n, *bptr = b->n;

	ASSERTG (bsize >= 0 && asize >= bsize);

	RES = CARRYL = CARRYH = 0;
	for (j = 0; j < asize; j++) {
		if (azeros)
			azeros--;
		else {
			addhlp (*aptr);
			aptr++;
		}
		if (j < bsize) subhlp (*bptr);
		*bptr++ = RES;
		RES = CARRYL;
		CARRYL = CARRYH;
	}

	while (j > 0 && b->n[j-1] == 0) j--;
	b->sign = j;

	ASSERTG (b->sign == 0 || b->n[abs(b->sign)-1] != 0);
}

void addg (		/* b := b + a, any signs any result. */
	giant	a,
	giant	b)
{
	int 	asgn = a->sign, bsgn = b->sign;

	ASSERTG (a->sign == 0 || a->n[abs(a->sign)-1] != 0);
	ASSERTG (b->sign == 0 || b->n[abs(b->sign)-1] != 0);

	if (asgn == 0) return;
	if (bsgn == 0) {
		gtog (a, b);
		return;
	}
	if ((asgn < 0) == (bsgn < 0)) {
		if (bsgn > 0) {
			normal_addg (a, b);
		} else {
			absg (b);
			if (a != b) absg (a);
			normal_addg (a, b);
			negg (b);
			if (a != b) negg (a);
		}
	} else if (bsgn > 0) {
		negg (a);
		if (gcompg (b, a) >= 0) {
			normal_subg (a, b);
			negg (a);
		} else {
			reverse_subg (a, b, 0);
			negg (a);
			negg (b);
		}
	} else {
		negg (b);
		if (gcompg(b,a) < 0) {
			reverse_subg (a, b, 0);
		} else {
			normal_subg (a, b);
			negg (b);
		}
	}

	ASSERTG (b->sign == 0 || b->n[abs(b->sign)-1] != 0);
	ASSERTG (abs (b->sign) <= b->maxsize);
}

void sladdg (			/* Giant g becomes g + i. */
	long	i,
	giant	g)
{
	giantstruct tmp;

	if (i == 0) return;

	if (i < 0) {i = -i; tmp.sign = -1;}
	else tmp.sign = 1;
	tmp.n = (unsigned long *) &i;
	setmaxsize (&tmp, 1);
	addg (&tmp, g);
}

void uladdg (			/* Giant g becomes g + i. */
	unsigned long i,
	giant	g)
{
	giantstruct tmp;

	if (i == 0) return;

	tmp.sign = 1;
	tmp.n = &i;
	setmaxsize (&tmp, 1);
	addg (&tmp, g);
}

void ulsubg (			/* Giant g becomes g - i. */
	unsigned long i,
	giant	g)
{
	giantstruct tmp;

	if (i == 0) return;

	tmp.sign = 1;
	tmp.n = &i;
	setmaxsize (&tmp, 1);
	subg (&tmp, g);
}

void subg (		/* b := b - a, any signs, any result. */
	giant	a,
	giant	b)
{
	int	asgn = a->sign, bsgn = b->sign;

	ASSERTG (a->sign == 0 || a->n[abs(a->sign)-1] != 0);
	ASSERTG (b->sign == 0 || b->n[abs(b->sign)-1] != 0);

	if (asgn == 0) return;
	if (bsgn == 0) {
		gtog (a, b);
		negg (b);
		return;
	}
	if ((asgn < 0) != (bsgn < 0)) {
		if (bsgn > 0) {
			negg (a);
			normal_addg (a, b);
			negg (a);
		} else {
			negg (b);
			normal_addg (a, b);
			negg (b);
		}
	} else if (bsgn > 0) {
		if (gcompg (b, a) >= 0) {
			normal_subg (a, b);
		} else {
			reverse_subg (a, b, 0);
			negg (b);
		}
	} else {
		negg (a);
		negg (b);
		if (gcompg (b, a) >= 0) {
			normal_subg (a, b);
			negg (a);
			negg (b);
		} else {
			reverse_subg (a, b, 0);
			negg (a);
		}
	}

	ASSERTG (b->sign == 0 || b->n[abs(b->sign)-1] != 0);
	ASSERTG (abs (b->sign) <= b->maxsize);
}

void setmulmode (
	int 	mode)
{
	mulmode = mode;
}

void squareg (			/* b becomes b*b */
	giant	b
)
/* Optimized general square, b becomes b*b. Modes are:
 * AUTO_MUL: switch according to empirical speed criteria.
 * GRAMMAR_MUL: force grammar-school algorithm.
 * KARAT_MUL: force Karatsuba divide-conquer method.
 * FFT_MUL: force floating point FFT method. */
{
	int 	bsize;

	ASSERTG (b->sign == 0 || b->n[abs(b->sign)-1] != 0);

	if (b->sign == 0) return;

	absg (b);
	switch (mulmode) {
	case GRAMMAR_MUL:
		grammarsquareg (b);
		break;
	case FFT_MUL:
		FFTsquareg (b);
		break;
	case KARAT_MUL:
		karatsquareg (b);
		break;
	case AUTO_MUL:
		bsize = b->sign;
		if (bsize >= FFT_BREAK_SQUARE)
			FFTsquareg (b);
		else if (bsize >= KARAT_BREAK_SQUARE)
			karatsquareg (b);
		else
			grammarsquareg (b);
		break;
	}

	ASSERTG (b->sign > 0 && b->n[b->sign-1] != 0);
}

void mulg (			/* b becomes a*b */
	giant	a,
	giant	b)
/* Optimized general multiply, b becomes a*b. Modes are:
 * AUTO_MUL: switch according to empirical speed criteria.
 * GRAMMAR_MUL: force grammar-school algorithm.
 * KARAT_MUL: force Karatsuba divide-conquer method.
 * FFT_MUL: force floating point FFT method. */
{
	int 	neg, asign;

	ASSERTG (a->sign == 0 || a->n[abs(a->sign)-1] != 0);
	ASSERTG (b->sign == 0 || b->n[abs(b->sign)-1] != 0);

	if (a == b) {
		squareg (b);
		return;
	}

	if (a->sign == 0 || b->sign == 0) {
		b->sign = 0;
		return;
	}

	neg = 0;
	asign = a->sign;
	if (a->sign < 0) { neg = 1; a->sign = -a->sign; }
	if (b->sign < 0) { neg = !neg; b->sign = -b->sign; }

	switch (mulmode) {
	case GRAMMAR_MUL:
		grammarmulg (a, b);
		break;
	case FFT_MUL:
		FFTmulg (a, b);
		break;
	case KARAT_MUL:
		karatmulg (a, b);
		break;
	case AUTO_MUL:
		automulg (a, b);
		break;
	}

	if (asign < 0) a->sign = -a->sign;	/* Restore a's sign */
	if (neg) b->sign = -b->sign;
	ASSERTG (b->sign != 0 && b->n[abs(b->sign)-1] != 0);
	ASSERTG (abs (b->sign) <= b->maxsize);
}

void ulmulg (			/* Giant g becomes g * i. */
	unsigned long i,
	giant	g)
{
	giantstruct tmp;

	if (i == 0) {g->sign = 0; return;}
	if (i == 1) return;

	tmp.sign = 1;
	tmp.n = &i;
	setmaxsize (&tmp, 1);
	grammarmulg (&tmp, g);
}

void imulg (			/* Giant g becomes g * i. */
	long	i,
	giant	g)
{
	giantstruct tmp;
	unsigned long data;

	if (i == 0) {g->sign = 0; return;}
	if (i == 1) return;

	tmp.n = &data;
	itog (i, &tmp);
	mulg (&tmp, g);
}

void dblmulg (			/* Giant g becomes g * i. */
	double	i,
	giant	g)
{
	giantstruct tmp;
	unsigned long tmparray[2];

	tmp.n = (unsigned long *) &tmparray;
	setmaxsize (&tmp, 2);
	dbltog (i, &tmp);
	grammarmulg (&tmp, g);
}

void modg (		/* n becomes n%d. n is arbitrary, but the
			 * denominator d must be positive! returned
			 * n will always be positive. */
	giant 	d,
	giant 	n)
{
	giant	tmp;

	ASSERTG (d->sign > 0);
	ASSERTG (d->n[d->sign-1] != 0);
	ASSERTG (n->sign == 0 || n->n[abs(n->sign)-1] != 0);

/* Use divg to compute the mod */

	tmp = popg (abs(n->sign));
	gtog (n, tmp);
	divg (d, tmp);
	mulg (d, tmp);
	subg (tmp, n);
	if (n->sign < 0) addg (d, n);
	pushg (1);

	ASSERTG (n->sign >= 0);
	ASSERTG (n->sign == 0 || n->n[n->sign-1] != 0);
}

void dbldivg (			/* Giant g becomes g / i. */
	double	i,
	giant	g)
{
	giantstruct tmp;
	unsigned long tmparray[2];

	tmp.n = (unsigned long *) &tmparray;
	setmaxsize (&tmp, 2);
	dbltog (i, &tmp);
	divg (&tmp, g);
}

void divg (		/* n becomes n/d. n is arbitrary, but the
			 * denominator d must be positive! */
	giant 	d,
	giant 	n)
{
	int	nsign, nsize, dsize;
	giant	r;

	ASSERTG (d->sign > 0);
	ASSERTG (d->n[d->sign-1] != 0);
	ASSERTG (n->sign == 0 || n->n[abs(n->sign)-1] != 0);

/* Handle simple cases */

	if (isone (d)) return;
	nsign = n->sign;
	nsize = abs (nsign);
	dsize = d->sign;
	if (nsize < dsize) {
		n->sign = 0;
		return;
	}

/* Handle cases where we can (probably) guess the quotient using the FPU */
/* I say probably because the floating point result is accurate to 53 bits. */
/* If, for example, qflt is exactly 3 then it could really be */
/* 2.999999999999999995 and we need to do a full divide to find out whether */
/* we should return 2 or 3. */
/* This speedup is especially important during GCDs which tend to generate */
/* very small quotients. */

	if (nsize - dsize <= 1) {
		double	dflt, nflt, qflt;
		unsigned long q;

		nflt = n->n[nsize-1];
		if (nsize > 1) nflt = nflt * 4294967296.0 + n->n[nsize-2];
		if (nsize > 2) nflt = nflt * 4294967296.0 + n->n[nsize-3];
		if (nsize > dsize && nsize > 3) nflt *= 4294967296.0;
		dflt = d->n[dsize-1];
		if (dsize > 1) dflt = dflt * 4294967296.0 + d->n[dsize-2];
		if (dsize > 2) dflt = dflt * 4294967296.0 + d->n[dsize-3];

		qflt = nflt / dflt;
		if (qflt < 4294967295.0) {
			q = (unsigned long) qflt;
			qflt = qflt - (double) q;
			if (nsize == 1 ||	/* We have an exact result */
			    (qflt >= TWO_TO_MINUS_19 &&
			     qflt <= ONE_MINUS_TWO_TO_MINUS_19)) {
				n->sign = (q == 0) ? 0 : 1;
				n->n[0] = q;
				if (nsign < 0) negg (n);
				return;
			}
		}
	}

/* The reciprocal code is very inefficient when n is large and d is small. */
/* Compensate for that here by using "schoolboy" division */

	if (nsize > dsize + dsize) {
		giant	tmp;
		int	i, chunks;

		r = popg (dsize + dsize);
		tmp = popg (dsize + dsize);

		chunks = nsize / dsize - 1;

/* Shift n right and do the first chunk */

		n->sign -= chunks * dsize;
		n->n += chunks * dsize;

		gtog (n, r);
		divg (d, n);

		gtog (n, tmp);
		mulg (d, tmp);
		subg (tmp, r);

/* For each remaining chunk, shift the remainder r left and copy in words */
/* from the input number n. */
		
		for (i = 0; i < chunks; i++) {
			gshiftleft (dsize << 5, r);
			n->n -= dsize;
			n->sign += dsize;
			memcpy ((char *) r->n,
				(char *) n->n,
				dsize * sizeof (unsigned long));
			if (r->sign == 0) {
				r->sign = dsize;
				while (r->sign && r->n[r->sign-1] == 0) r->sign--;
			}

/* Now compute our next chunk of result and copy it (zero padded) */
/* to our output giant n */

			gtog (r, tmp);
			divg (d, tmp);
			memset ((char *) n->n, 0,
				dsize * sizeof (unsigned long));
			memcpy ((char *) n->n,
				(char *) tmp->n,
				tmp->sign * sizeof (unsigned long));

/* Now compute the remainder for the next iteration */

			mulg (d, tmp);
			subg (tmp, r);
		}

		pushg (2);
	}

/* Divide using reciprocals */

	else if (d->sign < 30) {
		r = popg ((d->sign << 1) + 1);
		make_recip (d, r);
		divg_via_recip (d, r, n);
		pushg (1);
	}
	else {
		if (cur_recip == NULL || gcompg (d, cur_den)) {
			free (cur_recip);
			free (cur_den);
			cur_recip = newgiant ((d->sign << 1) + 1);
			cur_den = newgiant (d->sign << 1);
			gtog (d, cur_den);
			make_recip (d, cur_recip);
		}
		divg_via_recip (d, cur_recip, n);
	}

	ASSERTG (n->sign == 0 || n->n[abs(n->sign)-1] != 0);
}

void powerg (			/* x becomes x^n, NO mod performed. */
	giant	x,
	giant	n)
{
	int 	len;
	giant	scratch;

	ASSERTG (x->sign > 0);

	scratch = popg (x->sign << 1);
	gtog (x, scratch);
	for (len = bitlen (n) - 1; len; len--) {
		squareg (x);
		if (bitval (n, len-1)) mulg (scratch, x);
	}
	pushg (1);
}

void power (			/* x becomes x^n. */
	giant	x,
	int 	n)
{
	giant ng;

	if (istwo (x)) {
		itog (1, x);
		gshiftleft (n, x);
	} else {
		ng = popg (1);
		itog (n, ng);
		powerg (x, ng);
		pushg (1);
	}
}

void powermodg (		/* x becomes x^n (mod g). */
	giant	x,
	giant	n,
	giant	g)
{
	int 	len;
	giant	scratch;

	ASSERTG (x->sign > 0);
	ASSERTG (g->sign > 0);

	scratch = popg (g->sign << 1);
	gtog (x, scratch);
	for (len = bitlen (n) - 1; len; len--) {
		squareg (scratch);
		modg (g, scratch);
		if (bitval (n, len-1)) {
			mulg (x, scratch);
			modg (g, scratch);
		}
	}
	gtog (scratch, x);
	pushg (1);
}

void powermod (			/* x becomes x^n (mod g). */
	giant	x,
	int 	n,
	giant 	g)
{
	giant ng;

	ng = popg (1);
	itog (n, ng);
	powermodg (x, ng, g);
	pushg (1);
}

void automulg (			/* b becomes a*b */
	giant	a,
	giant	b)
{
	int 	asize, bsize;

	ASSERTG (a->sign >= 0);
	ASSERTG (b->sign >= 0);
	ASSERTG (a->sign == 0 || a->n[a->sign-1] != 0);
	ASSERTG (b->sign == 0 || b->n[b->sign-1] != 0);

/* Use grammar multiply for small arguments */

	asize = a->sign;
	bsize = b->sign;
	if (asize < KARAT_BREAK_MULT || bsize < KARAT_BREAK_MULT) {
		grammarmulg (a, b);
	}

/* This special case when one argument is more than roughly twice as big */
/* as the other (which happens frequently in ggcd) is not much faster than */
/* one big multiply, but it does use less memory which is important. */

	else if (asize + asize <= bsize + 2) {
		giant	d;

		d = popg (bsize);	/* d is upper half of b */
		gtogshiftright (asize << 5, b, d);
		b->sign = asize;	/* b is lower half of b */
		while (b->sign && b->n[b->sign-1] == 0) b->sign--;

		automulg (a, d);	/* Compute a * upper part of b */
		if (b->sign)
			automulg (a, b);/* Compute a * lower part of b */
		else {
			memset (b->n, 0, asize * sizeof (unsigned long));
			b->sign = asize;
		}

		b->sign -= asize;	/* Trick to add shifted upper */
		b->n += asize;
		normal_addg (d, b);
		b->sign += asize;	/* Undo the trick */
		b->n -= asize;

		pushg (1);
		return;
	}

/* Do a Karatsuba or FFT multiply */

	else if ((asize < FFT_BREAK_MULT1 && bsize < FFT_BREAK_MULT1) ||
	         (asize >= FFT_BREAK_MULT2 && bsize >= FFT_BREAK_MULT2 &&
	          asize < FFT_BREAK_MULT3 && bsize < FFT_BREAK_MULT3))
		karatmulg (a, b);
	else
		FFTmulg (a, b);

	ASSERTG (b->sign != 0 && b->n[abs(b->sign)-1] != 0);
}

void grammarsquareg (		/* a := a^2. */
	giant	a)
{
	int	i, asize = a->sign, max = asize * 2 - 1;
	unsigned long *ptr, *ptr1, *ptr2;
	giant	scratch;
	
	ASSERTG (a->sign >= 0);
	ASSERTG (a->sign == 0 || a->n[a->sign-1] != 0);

	if (asize == 0) return;

	scratch = popg (asize);
	gtog (a, scratch);
	ptr = scratch->n;

	asize--;
	RES = CARRYH = CARRYL = 0;
	for (i = 0; i < max; i++) {
		if (i <= asize) {
			ptr1 = ptr;
			ptr2 = ptr + i;
		} else {
			ptr1 = ptr + i - asize;
			ptr2 = ptr + asize;
		}
		while (ptr1 < ptr2) {
			muladd2hlp (*ptr1, *ptr2);
			ptr1++;
			ptr2--;
		}
		if (ptr1 == ptr2) {
			muladdhlp (*ptr1, *ptr1);
		}
		a->n[i] = RES;
		RES = CARRYL;
		CARRYL = CARRYH;
		CARRYH = 0;
	}
	if (RES) {
		a->n[i] = RES;
		a->sign = i+1;
	} else
		a->sign = i;

	pushg (1);
	ASSERTG (a->sign == 0 || a->n[a->sign-1] != 0);
}

void grammarmulg (		/* b becomes a*b. */
	giant	a,
	giant	b)
{
	int 	i, max;
	int 	asize = a->sign, bsize = b->sign;
	unsigned long *aptr, *bptr, *destptr;
	giant	scratch;

	ASSERTG (a->sign >= 0);
	ASSERTG (a->sign == 0 || a->n[a->sign-1] != 0);
	ASSERTG (b->sign >= 0);
	ASSERTG (b->sign == 0 || b->n[b->sign-1] != 0);

	if (bsize == 0) return;
	if (asize == 0) {
		b->sign = 0;
		return;
	}
	if (asize == 1 && a->n[0] == 1)
		return;

	scratch = popg (bsize);
	gtog (b, scratch);

	destptr = b->n;
	max = asize + bsize - 1;
	bsize--;
	RES = CARRYH = CARRYL = 0;
	for (i = 0; i < max; i++) {
		if (i <= bsize) {
			aptr = a->n;
			bptr = scratch->n + i;
		} else {
			aptr = a->n + i - bsize;
			bptr = scratch->n + bsize;
		}
		while (aptr < a->n + asize && bptr >= scratch->n) {
			muladdhlp (*aptr, *bptr);
			aptr++;
			bptr--;
		}
		*destptr++ = RES;
		RES = CARRYL;
		CARRYL = CARRYH;
		CARRYH = 0;
	}
	if (RES) {
		*destptr++ = RES;
		b->sign = i+1;
	} else
		b->sign = i;

	pushg (1);
	ASSERTG (b->sign == 0 || b->n[b->sign-1] != 0);
}

void karatsquareg (		/* x becomes x^2. */
	giant	x)
{
	int	s = x->sign, w;
	giantstruct a, b;
	giant	c;

	ASSERTG (x->sign >= 0);
	ASSERTG (x->sign == 0 || x->n[x->sign-1] != 0);

				/* IDEA: save a memcpy by having b point */
				/* to x's memory and writing its square to */
				/* the upper half of x */
	w = (s + 1) / 2;

	a.sign = w;		/* a is the lower half of the number */
	a.n = x->n;
	setmaxsize (&a, w);
	while (a.sign && a.n[a.sign-1] == 0) a.sign--;

	b.sign = s - w;		/* b is the upper half of the number */
	b.n = x->n + w + w;
	setmaxsize (&b, b.sign);
	memmove ((char *) b.n,
		 (char *) (a.n + w),
		 b.sign * sizeof (unsigned long));

	c = popg (w + w + 3);	/* c is the upper and lower half added */
	gtog (&a, c);
	normal_addg (&b, c);

	if (a.sign < KARAT_BREAK_SQUARE)
		grammarsquareg (&a);	/* recurse */
	else
		karatsquareg (&a);	/* recurse */
	if (b.sign < KARAT_BREAK_SQUARE)
		grammarsquareg (&b);
	else
		karatsquareg (&b);
	if (b.sign < KARAT_BREAK_SQUARE)
		grammarsquareg (c);
	else
		karatsquareg (c);

	subg (&a, c);		/* Isolate 2 * upper * lower */
	subg (&b, c);

	while (a.sign < w + w)	/* zero pad squared lower half */
		a.n[a.sign++] = 0;

	x->sign = b.sign + w;	/* Trick to add in a shifted value of c */
	x->n = a.n + w;
	normal_addg (c, x);
	x->sign += w;		/* Undo the trick */
	x->n = a.n;

	pushg (1);
	ASSERTG (x->sign == 0 || x->n[x->sign-1] != 0);
}

/* Next, improved Karatsuba routines from A. Powell. */

void karatmulg (		/* y becomes x*y. */
	giant	x,
	giant	y)
{
	int	s = x->sign, t = y->sign, w;
	giantstruct a, b, c;
	giant	d, e, f;

	ASSERTG (x->sign >= 0);
	ASSERTG (x->sign == 0 || x->n[x->sign-1] != 0);
	ASSERTG (y->sign >= 0);
	ASSERTG (y->sign == 0 || y->n[y->sign-1] != 0);

	w = (s + t + 2) / 4;

	a.sign = (s < w) ? s : w;	/* a is the lower half of x */
	a.n = x->n;
	setmaxsize (&a, a.sign);
	while (a.sign && a.n[a.sign-1] == 0) a.sign--;

	b.sign = (s > w) ? s - w : 0;	/* b is the upper half of x */
	b.n = a.n + w;
	setmaxsize (&b, b.sign);

	c.sign = (t < w) ? t : w;	/* c is the lower half of y */
	c.n = y->n;
	setmaxsize (&c, y->maxsize);
	while (c.sign && c.n[c.sign-1] == 0) c.sign--;

	d = popg (s + t);		/* d is the upper half of y */
	d->sign = (t > w) ? t - w : 0;
	memcpy ((char *) d->n,
		(char *) (c.n + w),
		d->sign * sizeof (unsigned long));

	e = popg (s + t);		/* e is the x halves added */
	gtog (&a, e); normal_addg (&b, e);

	f = popg (s + t + 1);		/* f is the y halves added */
	gtog (&c, f); normal_addg (d, f);

	if (a.sign < KARAT_BREAK_MULT || c.sign < KARAT_BREAK_MULT)
		grammarmulg (&a, &c);		/* Recurse: mul lowers */
	else
		karatmulg (&a, &c);		/* Recurse: mul lowers */
	if (b.sign < KARAT_BREAK_MULT || d->sign < KARAT_BREAK_MULT)
		grammarmulg (&b, d);		/* mul uppers */
	else
		karatmulg (&b, d);		/* mul uppers */
	if (e->sign < KARAT_BREAK_MULT || f->sign < KARAT_BREAK_MULT)
		grammarmulg (e, f);		/* mul sums */
	else
		karatmulg (e, f);		/* mul sums */

	normal_subg (&c, f);		/* Isolate cross product */
	normal_subg (d, f);

	if (d->sign) {
		memcpy (c.n + w + w,	/* Copy muled uppers to result */
			d->n,
			d->sign * sizeof (unsigned long));
		while (c.sign < w + w)	/* zero pad mul'ed lowers */
			c.n[c.sign++] = 0;
		y->sign = d->sign + w;	/* Trick to add shifted f */
	} else {
		while (c.sign < w)	/* zero pad mul'ed lowers */
			c.n[c.sign++] = 0;
		y->sign = c.sign - w;	/* Trick to add shifted f */
	}
	y->n = c.n + w;
	normal_addg (f, y);
	y->sign += w;			/* Undo the trick */
	y->n = c.n;

	pushg (3);
	ASSERTG (y->sign == 0 || y->n[y->sign-1] != 0);
	ASSERTG (y->sign <= y->maxsize);
}

void gmaskbits (	/* keep rightmost bits. Equivalent to g = g%2^bits. */
	int	bits,
	giant	g)
{
	int 	rem = bits&31, words = bits>>5;
	int 	size = abs (g->sign);

	ASSERTG (bits >= 0);
	ASSERTG (g->sign >= 0 || g->n[abs(g->sign)-1] != 0);

	if (size <= words) return;

	if (rem) {
		g->n[words] &= ((1 << rem) - 1);
		words++;
	}
	g->sign = (g->sign < 0) ? -words : words;
	while (g->sign && g->n[g->sign-1] == 0) g->sign--;
}

void gshiftleft (	/* shift g left bits. Equivalent to g = g*2^bits. */
	int	bits,
	giant	g)
{
	int 	rem = bits&31, crem = 32-rem, words = bits>>5;
	int 	size = abs (g->sign), j;
	unsigned long carry;

	ASSERTG (bits >= 0);
	ASSERTG (g->sign == 0 || g->n[abs(g->sign)-1] != 0);

	if (bits == 0) return;
	if (size == 0) return;

	if (rem == 0) {
		memmove (g->n + words, g->n, size * sizeof (unsigned long));
		memset (g->n, 0, words * sizeof (unsigned long));
		g->sign += (g->sign < 0) ? -words : words;
		return;
	}
	carry = g->n[size-1] >> crem;
	for (j = size-1; j > 0; j--) {
		g->n[j+words] = (g->n[j] << rem) | (g->n[j-1] >> crem);
	}
	g->n[words] = g->n[0] << rem;
	memset (g->n, 0, words * sizeof (unsigned long));
	size = size + words;
	if (carry) g->n[size++] = carry;
	g->sign = (g->sign < 0) ? -size : size;

	ASSERTG (g->sign == 0 || g->n[abs(g->sign)-1] != 0);
}

void gtogshiftright (	/* shift src right. Equivalent to dest = src/2^bits. */
	int	bits,
	giant	src,
	giant	dest)
{
	register int j, size = abs (src->sign);
	register unsigned long carry, *sptr, *dptr;
	int 	words = bits >> 5;
	int 	remain = bits & 31, cremain = (32-remain);

	ASSERTG (bits >= 0);
	ASSERTG (bits || src != dest);
	ASSERTG (src->sign == 0 || src->n[abs(src->sign)-1] != 0);

	if (words >= size) {
		dest->sign = 0;
		return;
	}
	if (remain == 0) {
		memcpy (dest->n,
			src->n + words,
			(size - words) * sizeof (unsigned long));
		dest->sign = src->sign + (src->sign < 0 ? words : -words);
		return;
	}

	size -= words;
	sptr = src->n + words;
	dptr = dest->n;
	for (j = 0; j < size-1; j++, sptr++, dptr++) {
		carry = sptr[1] << cremain;
		*dptr = (sptr[0] >> remain) | carry;
	}
	*dptr = *sptr >> remain;
	if (*dptr == 0) size--;
	dest->sign = (src->sign > 0) ? size : -size;

	ASSERTG (dest->sign == 0 || dest->n[abs(dest->sign)-1] != 0);
}

void invg (		/* Computes 1/y, that is the number n such that */
			/* n * y mod x = 1.  If x and y are not */
			/* relatively prime, y is the -1 * GCD (x, y). */
	giant 	xx,
	giant 	yy)
{
	invg_common (xx, yy, 0);
}

int invgi (		/* Interruptable version of invg */
	giant 	xx,
	giant 	yy)
{
	return (invg_common (xx, yy, 1));
}

int invg_common (	/* Common invg code */
	giant 	xx,
	giant 	yy,
	int	interruptable)
{
	giant	x, y;
	gmatrixstruct A;
	giantstruct ul, ll;
	int	ss;

	ASSERTG (xx->sign > 0 && xx->n[xx->sign-1] != 0);
	ASSERTG (yy->sign != 0 && yy->n[abs(yy->sign)-1] != 0);

/* Copy the first argument, we can trash the second */

	ss = stackstart ();
	x = popg (xx->sign); gtog (xx, x);
	y = yy;

/* Make y positive and less than x */

	while (y->sign < 0) addg (x, y);
	modg (x, y);

/* Compute the GCD and inverse the fastest way possible */
/* The inverse is computed in the matrix A, and only the */
/* right side of the matrix is needed.  However, the recursive */
/* ggcd code needs the left side of the array allocated. */

	A.ur = popg (x->sign); setzero (A.ur);
	A.lr = popg (x->sign); setone (A.lr);
	if (y->sign <= GCDLIMIT) {
		A.ul = &ul; setzero (A.ul); setmaxsize (A.ul, 1);
		A.ll = &ll; setzero (A.ll); setmaxsize (A.ll, 1);
		if (!cextgcdg (&x, &y, &A, interruptable)) goto esc;
	} else {
		A.ul = popg (x->sign); setone (A.ul);
		A.ll = popg (x->sign); setzero (A.ll);
		if (!ggcd (&x, &y, &A, interruptable)) goto esc;
		pushg (2);
	}

/* If a GCD was found, return that in yy (times -1) */

	if (! isone (x)) {
		if (x != yy) gtog (x, yy);
		yy->sign = -yy->sign;
		ASSERTG (yy->sign < 0 && yy->n[-yy->sign-1] != 0);
	}

/* Otherwise, return the inverse of yy in yy */

	else {
		gtog (A.ur, yy);
		if (A.ur->sign < 0) addg (xx, yy);
		ASSERTG (yy->sign > 0 && yy->n[yy->sign-1] != 0);
	}

/* Cleanup and return */

	pushg (3);
	return (TRUE);
esc:	pushall (ss);
	return (FALSE);
}

/* A wrapper routine for the assembly language egcdhlp routine */
/* Egcdhlp returns A,B,C,D as defined in Knuth vol. 2's description */
/* of extended GCD for large numbers. */

int egcdhlp_wrapper (
	giant	u,
	giant	v,
	gmatrix	A)
{
	int	i;
	unsigned long r, ch, cl;

	SRCARG = u->n;
	SRC2ARG = v->n;
	NUMARG = u->sign;
	NUM2ARG = v->sign;
	if (!egcdhlp ()) return (FALSE);
	if (A != NULL) {
		gmatrixstruct B;
		giantstruct ul, ur, ll, lr;
		setmaxsize (&ul, 1); setmaxsize (&ur, 1);
		setmaxsize (&ll, 1); setmaxsize (&lr, 1);
		B.ul = &ul; B.ur = &ur;
		B.ll = &ll; B.lr = &lr;
		ul.sign = EGCD_A ? 1 : 0; ul.n = &EGCD_A;
		ur.sign = EGCD_B ? 1 : 0; ur.n = &EGCD_B;
		ll.sign = EGCD_C ? 1 : 0; ll.n = &EGCD_C;
		lr.sign = EGCD_D ? 1 : 0; lr.n = &EGCD_D;
		if (EGCD_ODD) {
			ul.sign = -ul.sign;
			lr.sign = -lr.sign;
		} else {
			ur.sign = -ur.sign;
			ll.sign = -ll.sign;
		}
		mulmM (&B, A);
	}

/* Now do a mulvM (&B, u, v) using the special known properties */
/* of matrix B.  Not only is this code faster than mulvM, but it */
/* uses less memory and does not require u and v to be copied */
/* prior to calling mulvM. */

	RES = CARRYH = CARRYL = 0;
	r = ch = cl = 0;
	for (i = 0; i < v->sign; i++) {
		if (EGCD_ODD) {
			mulsubhlp (u->n[i], EGCD_A);
			muladdhlp (v->n[i], EGCD_B);
		} else {
			muladdhlp (u->n[i], EGCD_A);
			mulsubhlp (v->n[i], EGCD_B);
		}
		ulswap (RES, r); ulswap (CARRYH, ch); ulswap (CARRYL, cl);
		if (EGCD_ODD) {
			muladdhlp (u->n[i], EGCD_C);
			mulsubhlp (v->n[i], EGCD_D);
		} else {
			mulsubhlp (u->n[i], EGCD_C);
			muladdhlp (v->n[i], EGCD_D);
		}
		ulswap (RES, r); ulswap (CARRYH, ch); ulswap (CARRYL, cl);
		u->n[i] = RES;
		RES = CARRYL;
		CARRYL = CARRYH;
		CARRYH = (int) CARRYH >> 4;
		v->n[i] = r;
		r = cl;
		cl = ch;
		ch = (int) ch >> 4;
	}
	u->sign = v->sign;
	while (u->sign > 0 && u->n[u->sign-1] == 0) u->sign--;
	while (v->sign > 0 && v->n[v->sign-1] == 0) v->sign--;
	ASSERTG (u->sign > v->sign ||
		 (u->sign == v->sign && u->n[u->sign-1] >= v->n[v->sign-1]));
	return (TRUE);
}


void gcdg (	/* Computes the GCD of x and y and returns the GCD in y */
		/* The x argument is not destroyed */
	giant	xx,
	giant	yy)
{
	gcdg_common (xx, yy, 0);
}

int gcdgi (	/* Interruptable version of the aove */
	giant	xx,
	giant	yy)
{
	return (gcdg_common (xx, yy, 1));
}

int gcdg_common (	/* Common code for above */
	giant	xx,
	giant	yy,
	int	interruptable)
{
	giant	x, y;
	int	ss;

	ASSERTG (xx->sign > 0 && xx->n[xx->sign-1] != 0);
	ASSERTG (yy->sign > 0 && yy->n[yy->sign-1] != 0);

/* Copy the first argument, we can trash the second */

	ss = stackstart ();
	x = popg (xx->sign + 1); gtog (xx, x);
	y = yy;

/* Make y less than x */

	if (gcompg (x, y) < 0) gswap (&x, &y);

/* Compute the GCD the fastest way possible */

	if (abs (y->sign) <= GCDLIMIT) {
		if (!cextgcdg (&x, &y, NULL, interruptable)) goto esc;
	} else {
		if (!ggcd (&x, &y, NULL, interruptable)) goto esc;
	}

/* If the routines we called happened to return the result in yy, then great */
/* Otherwise, copy the result to yy */

	if (x != yy) gtog (x, yy);

/* Cleanup and return */

	ASSERTG (yy->sign > 0 && yy->n[yy->sign-1] != 0);
	pushg (1);
	return (TRUE);
esc:	pushall (ss);
	return (FALSE);
}

int cextgcdg (		/* Classical Euclid GCD. a becomes gcd(a, b). */
	giant 	*a,
	giant 	*b,
	gmatrix	A,
	int	interruptable)
{
	giant	u, v;

	ASSERTG ((*a)->sign > 0);
	ASSERTG ((*b)->sign >= 0);

/* Dereference arguments for faster access */

	u = *a;
	v = *b;

/* Do the extended GCD */

	while (v->sign > 0) {
		unsigned int diff_length;

/* Check for an interrupt */

		if (interruptable && stopCheck ()) return (FALSE);

/* If the quotient will fit in one word, use the GCD helper */
/* function which can do several GCD steps in single precision, */
/* postponing multi-precision operations as long as possible. */

/* There is a remote chance that the GCD helper routine can not */
/* correctly determine a single quotient.  In that case egcdhlp */
/* returns FALSE and we use the slower code that can compute */
/* the proper quotient. */

		diff_length = u->sign - v->sign;
		if (diff_length > 1 ||
		    (diff_length == 1 && u->n[u->sign-1] >= v->n[v->sign-1]) ||
		    ! egcdhlp_wrapper (u, v, A))
			onestep (&u, &v, A);
	}

/* Dereference arguments for faster access */

	*a = u;
	*b = v;

/* Cleanup and return */

	return (TRUE);
}


/**************************************************************
 *
 * Initialization and utility functions
 *
 **************************************************************/

#ifdef GDEBUG
unsigned long totmemsize = 0;
unsigned long maxtotmemsize = 0;
#endif

giant popg (
	int	size)
{
	giant	g;
	gstacknode *s, *sprev;
static	unsigned long gwnumsize = 0;
	unsigned long memsize;

	ASSERTG (size >= 0);

/* Initialize the stack if we're just starting out. */

	if (cur_stack_size == 0) {
		cur_stack_size = STACK_GROW;
		stack = (gstacknode *)
			malloc (cur_stack_size * sizeof (gstacknode));
		gwnumsize = 0;
	}

/* Expand the stack if we need to. */

	else if (cur_stack_elem >= cur_stack_size) {
		cur_stack_size += STACK_GROW;
		stack = (gstacknode *)
			realloc (stack, cur_stack_size * sizeof (gstacknode));
	}

/* Set gwnum_size if FFTs have been initialized.  This lets giants and gwnum */
/* share the same memory space.  Pretty handy in conserving allocated memory */

	if (gwnumsize == 0 && FFTLEN)
		gwnumsize = gwnum_size (FFTLEN);

/* Malloc our giant */

	s = &stack[cur_stack_elem];
	memsize = sizeof (giantstruct) + size * sizeof (unsigned long);
#ifdef GDEBUG
	totmemsize += memsize;
	if (totmemsize > maxtotmemsize)
		maxtotmemsize = totmemsize;
#endif
	if (memsize > gwnumsize) {
		s->gw = (gwnum) malloc (memsize);
		s->offset = -1;
		s->next_offset = 0;
		s->space_left = 0;
		g = (giant) s->gw;
	} else if (cur_stack_elem == 0 ||
		   (sprev = &stack[cur_stack_elem-1])->space_left < memsize) {
				/* IDEA:  Look for holes in previous gwnums */
		s->gw = gwalloc ();
		s->offset = 0;
		s->next_offset = memsize;
		s->space_left = gwnumsize - memsize;
		g = (giant) s->gw;
	} else {
		s->gw = sprev->gw;
		s->offset = sprev->next_offset;
		s->next_offset = s->offset + memsize;
		s->space_left = sprev->space_left - memsize;
		g = (giant) ((char *) s->gw + s->offset);
	}
	cur_stack_elem++;

	g->n = (unsigned long *) ((char *) g + sizeof (giantstruct));
	setmaxsize (g, size);

	return (g);
}

void pushg (
	int	a)
{
	gstacknode *s;

	for ( ; a; a--) {
		s = &stack[--cur_stack_elem];
#ifdef GDEBUG
		totmemsize -= s->next_offset - s->offset;
#endif
		if (s->offset == -1) free (s->gw);
		if (s->offset == 0) gwfree (s->gw);
	}
}

int bitlen (
	giant	g)
{
	int 	b = 32, c = 1<<31, w, size;

	ASSERTG (g->sign == 0 || g->n[abs(g->sign)-1] != 0);

	size = abs (g->sign);
	if (size == 0) return (0);
	w = g->n[size - 1];
	while ((w&c) == 0) {
		b--;
		c >>= 1;
	}
	return (32 * (size-1) + b);
}

int gsign (		/* Returns the sign of g. */
	giant 	g)
{
	if (isZero(g)) return (0);
	if (g->sign > 0) return (1);
	return (-1);
}


/**************************************************************
 *
 * Private Math Functions
 *
 **************************************************************/


void make_recip (	/* r becomes the steady-state reciprocal
			 * 2^(2b)/d, where b = bit-length of d-1. */
	giant 	d, 
	giant 	r)
{
	int	b;
	giant 	tmp, tmp2;

	ASSERTG (d->sign > 0);

	tmp = popg ((d->sign << 1) + 1);
	tmp2 = popg ((d->sign << 1) + 1);
	setone (r);
	normal_subg (r, d);
	b = bitlen (d);
	normal_addg (r, d);
	gshiftleft (b, r);
	gtog (r, tmp2);
	while (1) {
		gtog (r, tmp);
		squareg (tmp);
		gshiftright (b, tmp);
		mulg (d, tmp);
		gshiftright (b, tmp);
		addg (r, r); 
		subg (tmp, r);
		if (gcompg(r, tmp2) <= 0) 
			break;
		gtog (r, tmp2);
	}
	setone (tmp);
	gshiftleft (2*b, tmp);
	gtog (r, tmp2); 
	mulg (d, tmp2);
	subg (tmp2, tmp);
	setone (tmp2);
	while (tmp->sign < 0) {
		subg (tmp2, r);
		addg (d, tmp);
	}
	pushg(2);
}

void divg_via_recip (		/* n := n/d, where r is the precalculated
				 * steady-state reciprocal of d. */
	giant 	d, 
	giant 	r, 
	giant 	n)
{
	int 	s = 2*(bitlen(r)-1), sign = gsign(n);
	giant 	tmp, tmp2;

	ASSERTG (d->sign > 0);
	
	tmp = popg (abs (n->sign) + r->sign);
	tmp2 = popg (abs (n->sign) + r->sign);
	
	n->sign = abs (n->sign);
	setzero (tmp2);
	while (1) {
		gtog (n, tmp);	
		mulg (r, tmp);
		gshiftright (s, tmp);
		addg (tmp, tmp2);
		mulg (d, tmp);
		subg (tmp, n);
		if (gcompg(n,d) >= 0) {
			subg (d,n);
			iaddg (1, tmp2);
		}
		if (gcompg (n,d) < 0) break;
	}
	gtog (tmp2, n);
	n->sign *= sign;
	pushg (2);
}


/**************************************************************
 *
 * FFT multiply Functions
 *
 **************************************************************/

int lpt (		/* Returns least power of two greater than n. */
	int	n)
{
	register int	i = 1;

	while (i < n) {
		i <<= 1;
	}
	return(i);
}

void makewt (int nw, int *ip, double *w);
void makect (int nc, int *ip, double *c);

void init_sincos (
	int 	n)
{
	if (n <= ooura_fft_size) return;
	aligned_free (ooura_sincos);
	ooura_sincos = (double *) aligned_malloc ((n/2) * sizeof (double), sizeof (double));
	if (ooura_ip) free (ooura_ip);
	ooura_ip = (int *) malloc (((int) sqrt((double)(n/2)) + 2) * sizeof (int));
	ooura_ip[0] = 0;
	ooura_fft_size = n;
        makewt (n >> 2, ooura_ip, ooura_sincos);
        makect (n >> 2, ooura_ip, ooura_sincos + (n >> 2));
}

void mp_squ_cmul(int nfft, double dinout[])
{
	int j;
	double xr, xi;
    
	dinout[0] *= dinout[0];
	dinout[1] *= dinout[1];
	for (j = 2; j < nfft; j += 2) {
		xr = dinout[j];
		xi = dinout[j + 1];
		dinout[j] = xr * xr - xi * xi;
		dinout[j + 1] = xr * xi * 2.0;
	}
}

void mp_mul_cmul(int nfft, double din[], double dinout[])
{
	int j;
	double xr, xi, yr, yi;
    
	dinout[0] *= din[0];
	dinout[1] *= din[1];
	for (j = 2; j < nfft; j += 2) {
		xr = din[j];
		xi = din[j + 1];
		yr = dinout[j];
		yi = dinout[j + 1];
		dinout[j] = xr * yr - xi * yi;
		dinout[j + 1] = xr * yi + xi * yr;
	}
}

void addsignal (
	giant	x,
	int	size,
	double 	*z,
	int 	n)
{
#define TWO24		((double)(16777216.0))
#define TWO24PLUSHALF	((double)(16777216.5))
#define TWOM24		(double)(0.000000059604644775390625)
	register int	j, m1, m2, value, carry;
	register double scale, f, fltcarry;

/* Scale each element down */

	scale = 2.0 / (double) n;

/* Convert each double to an integer.  Extract lower 16 bits and leave the */
/* rest as a carry.  Each pair of conversions results in one output value.*/
/* This is tricky as we can't convert the 53-bit double to a 31-bit int. */
/* Therefore we convert in two stages - the upper 29 bits and the */
/* lower 24-bits. */

	carry = 0;
	fltcarry = 0.0;
	for (j = 0; ; ) {
		f = z[j+j] * scale + fltcarry * 256.0;
		fltcarry = (double) ((int) (f * TWOM24));
		f = f - fltcarry * TWO24;
		m1 = (int) (f + TWO24PLUSHALF) + carry;
		carry = (m1 >> 16) - 256;

		f = z[j+j+1] * scale + fltcarry * 256.0;
		fltcarry = (double) ((int) (f * TWOM24));
		f = f - fltcarry * TWO24;
		m2 = (int) (f + TWO24PLUSHALF) + carry;
		carry = (m2 >> 16) - 256;

		value = ((m2 & 0xFFFF) << 16) + (m1 & 0xFFFF);
		if (++j == size) break;
		x->n[j-1] = value;
	}
	if (value == 0) j--;
	else x->n[j-1] = value;
	x->sign = j;
	ASSERTG (x->n[x->sign-1] != 0);
}

void FFTsquareg (
	giant	x)
{
	int	size = x->sign;
	double	*z1;
	register int 	L;

	ASSERTG (x->sign >= 4 && x->n[x->sign-1] != 0);

	L = lpt (size+size) << 1;
	init_sincos (L);
	z1 = (double *) aligned_malloc (L * sizeof (double), sizeof (double));

	giant_to_double (x, size, z1, L);
	rdft (L, 1, z1, ooura_ip, ooura_sincos);
	mp_squ_cmul (L, z1);
	rdft (L, -1, z1, ooura_ip, ooura_sincos);
	addsignal (x, size+size, z1, L);

	aligned_free (z1);
}

void FFTmulg (			/* x becomes y*x. */
	giant	y,
	giant	x)
{
	int	sizex = x->sign, sizey = y->sign;
	double	*z1, *z2;
	register int	L;

	ASSERTG (y->sign >= 4 && y->n[y->sign-1] != 0);
	ASSERTG (x->sign >= 4 && x->n[x->sign-1] != 0);

/* Do the FFT multiply.  Make sure the arrays of doubles are aligned on */
/* an eight byte boundaries. */

	L = lpt (sizex+sizey) << 1;
	init_sincos (L);
	z1 = (double *) aligned_malloc (L * sizeof (double), sizeof (double));
	z2 = (double *) aligned_malloc (L * sizeof (double), sizeof (double));

	giant_to_double (x, sizex, z1, L);
	giant_to_double (y, sizey, z2, L);
	rdft (L, 1, z1, ooura_ip, ooura_sincos);
	rdft (L, 1, z2, ooura_ip, ooura_sincos);
	mp_mul_cmul (L, z2, z1);
	rdft (L, -1, z1, ooura_ip, ooura_sincos);
	addsignal (x, sizex+sizey, z1, L);

	aligned_free (z1);
	aligned_free (z2);
	ASSERTG (y->sign > 0 && y->n[y->sign-1] != 0);
}

void giant_to_double (
	giant 	x,
	int 	sizex,
	double 	*z,
	int 	L)
{
	register int j, w, carry;

	carry = 0;
	for (j = 0; j < sizex; j++) {
		w = (x->n[j] & 0xFFFF) + carry;
		if (w < 32768) {
			z[j+j] = w;
			carry = 0;
		} else {
			z[j+j] = w - 65536;
			carry = 1;
		}
		w = (x->n[j] >> 16) + carry;
		if (w < 32768) {
			z[j+j+1] = w;
			carry = 0;
		} else {
			z[j+j+1] = w - 65536;
			carry = 1;
		}
	}
	if (carry) z[j+j-1] += 65536;
	for (j = sizex + sizex; j < L; j++) {
		z[j] = 0.0;
	}
}

void gsetlength (	/* Set the length of g to n bits (g = g mod 2^n) */
	int	n,
	giant	g)
{
	int	size;

	ASSERTG (g->sign >= (n >> 5));

	size = (n + 31) >> 5;
	if (n & 31) g->n[size-1] &= (1 << (n & 31)) - 1;
	while (size && g->n[size-1] == 0) size--;
	g->sign = size;
}

void addshiftedg (	/* Shift x left n words then add to g */
			/* This lets us allocate smaller temporaries in hgcd */
	int	n,
	giant	x,	/* x must be positive */
	giant	g)
{
	ASSERTG (n >= 0 && x->sign >= 0);

	if (x->sign == 0) return;

	if (g->sign >= 0) {
		while (g->sign < n) g->n[g->sign++] = 0;
		g->sign -= n;
		g->n += n;
		normal_addg (x, g);
		g->sign += n;
		g->n -= n;
		while (g->sign && g->n[g->sign-1] == 0) g->sign--;
		return;
	}

	negg (g);
	if (g->sign < n) {
		reverse_subg (x, g, n);
		return;
	} 

	g->sign -= n;
	g->n += n;
	if (gcompg (g, x) >= 0) {
		normal_subg (x, g);
		g->sign += n;
		g->n -= n;
		while (g->sign && g->n[g->sign-1] == 0) g->sign--;
		negg (g);
	} else {
		g->sign += n;
		g->n -= n;
		reverse_subg (x, g, n);
	}
}

void onestep (		/* Do one step of the euclidean algorithm and modify
			 * the matrix A accordingly. */
	giant 	*x,
	giant 	*y,
	gmatrix A)
{
	giant q = popg ((*x)->sign);

	ASSERTG ((*x)->sign >= (*y)->sign);

	gtog (*x, q);		/* Set q = x / y */ 
	divg (*y, q);
	if (A != NULL) punch (q, A);
	mulg (*y, q);		/* Now set x = x - q * y */
	subg (q, *x);

	gswap (x, y);
	pushg(1);
}

void mulvM (		/* Multiply vector by Matrix; changes x,y. */
			/* Caller must make sure x and y variables */
			/* can hold larger intermediate results */
	gmatrix A,
	giant 	x,
	giant 	y)
{
	giant s0 = popg (abs(A->ll->sign) + x->sign);
	giant s1 = popg (abs(A->ur->sign) + y->sign);

	gtog (x, s0);
	gtog (y, s1);
	mulg (A->ul, x); mulg (A->ur, s1); addg (s1, x);
	mulg (A->lr, y); mulg (A->ll, s0); addg (s0, y);

	pushg (2);
}

void mulmM (		/* Multiply matrix by Matrix; changes second matrix. */
	gmatrix A,
	gmatrix B)
{
	giant s0 = popg (abs(A->ll->sign) + abs(B->ur->sign));
	giant s1 = popg (abs(A->ur->sign) + abs(B->lr->sign));

	gtog (B->ul, s0);
	gtog (B->ll, s1);
	mulg (A->ul, B->ul); mulg (A->ur, s1); addg (s1, B->ul);
	mulg (A->lr, B->ll); mulg (A->ll, s0); addg (s0, B->ll);
	gtog (B->ur, s0);
	gtog (B->lr, s1);
	mulg (A->ul, B->ur); mulg (A->ur, s1); addg (s1, B->ur);
	mulg (A->lr, B->lr); mulg (A->ll, s0); addg (s0, B->lr);

	pushg (2);
}

void mulmMsp (		/* Like mulmM except that the data areas of A */
	gmatrix A,	/* are in the upper half B (see hgcd) */
	gmatrix B,
	int	maxsize)
{
	giant tmp0 = popg (maxsize);
	giant tmp1 = popg (maxsize);
	giant tmp2 = popg (maxsize);
	giant tmp3 = popg (maxsize);

	gtog (A->ur, tmp0);	/* Copy A->ur before the mul destroys it */
	gtog (A->lr, tmp2);	/* Copy A->lr before the mul destroys it */
	gtog (B->ur, tmp1);	/* Copy B->ur before the mul destroys it */
	mulg (A->ul, B->ur);	/* A->ul * B->ur (destroys A->ur & A->lr) */
	gtog (tmp0, tmp3);
	mulg (B->lr, tmp3);	/* A->ur * B->lr */
	addg (tmp3, B->ur);	/* B->ur = A->ul * B->ur + A->ur * B->lr */

	mulg (A->ll, tmp1);	/* A->ll * B->ur */
	mulg (tmp2, B->lr);	/* A->lr * B->lr */
	addg (tmp1, B->lr);	/* B->lr = A->ll * B->ur + A->lr * B->lr */

	mulg (B->ll, tmp0);	/* A->ur * B->ll */
	gtog (B->ul, tmp1);	/* Copy B->ul before the mul destroys it */
	gtog (A->ll, tmp3);	/* Copy A->ll before the mul destroys it */
	mulg (A->ul, B->ul);	/* A->ul * B->ul (destroys A->ul & A->ll) */
	addg (tmp0, B->ul);	/* B->ul = A->ul * B->ul + A->ur * B->ll */

	mulg (tmp3, tmp1);	/* A->ll * B->ul */
	mulg (tmp2, B->ll);	/* A->lr * B->ll */
	addg (tmp1, B->ll);	/* B->ll = A->ll * B->ul + A->lr * B->ll */

	pushg (4);
}

void punch (		/* Multiply the matrix A on the left by [0,1,1,-q]. */
	giant 	q,
	gmatrix A)
{
	giant tmp = popg (abs(A->lr->sign) + q->sign);

	gtog (q, tmp);
	mulg (A->ll, tmp);
	subg (tmp, A->ul);
	gswap (&A->ul, &A->ll);

	gtog (q, tmp);
	mulg (A->lr, tmp);
	subg (tmp, A->ur);
	gswap (&A->ur, &A->lr);

	pushg (1);
}

int ggcd (		/* A giant gcd.  Modifies its arguments. */
	giant 	*x,
	giant 	*y,
	gmatrix	R,
	int	interruptable)
{
	gmatrixstruct A;

/* To avoid continually expanding the sincos array, figure out (roughly) */
/* the maximum size table we will need and allocate it now. */

	if ((*x)->sign / 2 > FFT_BREAK_MULT1)
		init_sincos (lpt ((*x)->sign / 2) << 1);

/* If R is not NULL then we are doing an extended GCD.  Recursively */
/* do half GCDs and then multiply the matrices in reverse order for */
/* optimum efficiency. */

	if (R != NULL) {
		if (! hgcd (0, x, y, R, interruptable)) return (FALSE);
		return (rhgcd (x, y, R, interruptable));
	}

/* Call half GCDs until the numbers get pretty small.  The half GCD code */
/* is most efficient when computing 1/3 of the GCD result size.  But more */
/* importantly, the FFT code is most efficient when operating on a power */
/* of two words.  Thus, if x->size is near a power of two we make the */
/* traditional two hgcd calls each computing 1/4 of the result size. */
/* Otherwise, we make one hgcd call computing 1/3 of the result size. */

/* If FFTmulg ever implemented non-power-of-2 FFTs, then we likely would */
/* always do a 1/3 hgcd call */

	while ((*y)->sign > GCDLIMIT) {
		int 	quarter_size, third_size, a_size;

		quarter_size = ((*y)->sign - 1) >> 2;
		third_size = ((*y)->sign - 1) / 3;
		if (lpt (quarter_size) == lpt (third_size))
			a_size = third_size;
		else
			a_size = quarter_size;
		A.ul = popg (a_size); setone (A.ul);
		A.ll = popg (a_size); setzero (A.ll);
		A.ur = popg (a_size); setzero (A.ur);
		A.lr = popg (a_size); setone (A.lr);

/* Do the first recursion */

		if (! hgcd ((*y)->sign - (a_size + a_size + 1),
			    x, y, &A, interruptable))
			return (FALSE);

/* Do a single step if the hgcd call didn't make any progress */

		if (isone (A.lr)) {
			pushg (4);		/* Free matrix A */
			onestep (x, y, NULL);
			continue;
		}

/* If we did a 1/3 bits hgcd call, then we're done.  If we did a 1/4 bits */
/* hgcd call, then the next iteration will do the matching 1/4 bits hgcd call */
/* (since result x is now 3/4 size and 1/3 of 3/4 is 1/4! */

		pushg (4);		/* Free matrix A */
	}

/* Do the last few words in a brute force way */

	return (cextgcdg (x, y, NULL, interruptable));
}

int rhgcd (	/* recursive hgcd calls accumulating extended GCD info */
		/* when done. */
	giant	*x,
	giant 	*y,
	gmatrix R,
	int	interruptable)
{
	gmatrixstruct A;

	A.ul = popg ((*x)->sign); setone (A.ul);
	A.ll = popg ((*x)->sign); setzero (A.ll);
	A.ur = popg ((*x)->sign); setzero (A.ur);
	A.lr = popg ((*x)->sign); setone (A.lr);
	if ((*y)->sign <= GCDLIMIT) {
		if (! cextgcdg (x, y, &A, interruptable)) return (FALSE);
	} else {
		if (! hgcd (0, x, y, &A, interruptable)) return (FALSE);
		if (isone (A.lr)) onestep (x, y, &A);
		if (! rhgcd (x, y, &A, interruptable)) return (FALSE);
	}
	mulmM (&A, R);
	pushg (4);
	return (TRUE);
}

int hgcd (	/* hgcd(n,x,y,A) chops n words off x and y and computes the
		 * 2 by 2 matrix A such that A[x y] is the pair of terms
		 * in the remainder sequence starting with x,y that is
		 * half the size of x. */
	int 	n,
	giant	*xx,
	giant 	*yy,
	gmatrix A,
	int	interruptable)
{
	giant 	x, y;
	int	half_size;

	ASSERTG (n >= 0);

/* Don't do anything if y isn't more than n words long */

	if ((*yy)->sign <= n) return (TRUE);

/* Finagle x and y giant structures to emulate a shift right n words */

	x = *xx;
	x->sign -= n;
	x->n += n;
	setmaxsize (x, x->maxsize - n);
	y = *yy;
	y->sign -= n;
	y->n += n;
	setmaxsize (y, y->maxsize - n);

/* If the numbers are small or the size of the first quotient won't let */
/* split the numbers further, then do the final level of recursion */
/* using a brute force algorithm */

	half_size = (y->sign - 1) >> 1;
	if (half_size <= CHGCD_BREAK || x->sign - y->sign > half_size >> 1) {
		for ( ; ; ) {
			int	quot_length, a_length;

/* Determine if applying quotient to the matrix will put us over half_size. */
/* Be wary of the final matrix multiply of A in egcdhlp and onestep.  It */
/* adds two multiplied products together - so that a 32-bit quotient applied */
/* to a 20 word A->lr value can produce a 22 word result if the top word */
/* of A->lr is also 32 bits long. */

			quot_length = x->sign - y->sign;
			if (x->n[x->sign-1] >= y->n[y->sign-1]) quot_length++;
			a_length = abs (A->lr->sign);
			if (A->lr->n[a_length-1] & 0x80000000) a_length++;
			if (a_length + quot_length > half_size) break;
		
/* If the quotient will fit in one word, use the GCD helper */
/* function which can do several GCD steps in single precision, */
/* postponing multi-precision operations as long as possible. */

			if (quot_length > 1 || !egcdhlp_wrapper (x, y, A))
				onestep (&x, &y, A);
		}
	}

/* Otherwise do two recursions to compute the half-gcd */

	else {
		gmatrixstruct B;
		giantstruct ul, ur, ll, lr;
		int	a_size, b_size;

/* Check for an interrupt */

		if (interruptable && stopCheck ()) return (FALSE);

/* Do the first recursion */

		a_size = half_size >> 1;
		if (! hgcd (y->sign - (a_size + a_size + 1),
			    &x, &y, A, interruptable))
			return (FALSE);
		a_size = abs (A->lr->sign);

/* Do the second recursion.  Note how we use the upper half of the gmatrix A */
/* in order to save a lot of memory.  Be wary of the matrix multiply of B */
/* and A in mulmMsp.  It adds two multiplied products together - so that */
/* if B->lr's top word is 32 bits A->lr's top word is 32 bits, then a */
/* carry overflow will occur.  Also be careful that we don't ask hgcd to */
/* compute more than y->sign / 2 words. */

		b_size = half_size - a_size;
		if (A->lr->n[a_size-1] & 0x80000000) b_size--;
		ul.n = A->ul->n + a_size; setmaxsize (&ul, b_size);
		ur.n = A->ur->n + a_size; setmaxsize (&ur, b_size);
		ll.n = A->ll->n + a_size; setmaxsize (&ll, b_size);
		lr.n = A->lr->n + a_size; setmaxsize (&lr, b_size);
		B.ul = &ul; setone (&ul);
		B.ur = &ur; setzero (&ur);
		B.ll = &ll; setzero (&ll);
		B.lr = &lr; setone (&lr);
		if (b_size > ((y->sign - 1) >> 1)) b_size = (y->sign - 1) >> 1;
		if (! hgcd (y->sign - (b_size + b_size + 1),
			    &x, &y, &B, interruptable))
			return (FALSE);
		mulmMsp (&B, A, half_size);
	}

/* Copy the x and y values, then undo the changes we made to the input */
/* giants to emulate a shift right n words and instead point to just */
/* the lower n words. */

	if (n) {
		giant 	xinp, yinp, tmp;

		tmp = x; x = popg (x->sign); gtog (tmp, x);
		tmp = y; y = popg (y->sign); gtog (tmp, y);
		xinp = *xx;
		xinp->sign = n;
		xinp->n -= n;
		while (xinp->sign && xinp->n[xinp->sign-1] == 0) xinp->sign--;
		setmaxsize (xinp, xinp->maxsize + n);
		yinp = *yy;
		yinp->sign = n;
		yinp->n -= n;
		while (yinp->sign && yinp->n[yinp->sign-1] == 0) yinp->sign--;
		setmaxsize (yinp, yinp->maxsize + n);

/* Now apply the matrix A to the bits of xx and yy that were shifted off */
/* This lets us compute the final x and y values for the caller */

		mulvM (A, xinp, yinp);
		addshiftedg (n, x, xinp);
		addshiftedg (n, y, yinp);
		if (xinp->sign < 0) {
			negg (xinp); negg (A->ul); negg (A->ur);
		}
		if (yinp->sign < 0) {
			negg (yinp); negg (A->ll); negg (A->lr);
		}
		if (gcompg (xinp, yinp) < 0) {
			gswap (xx, yy);
			gswap (&A->ul, &A->ll);
			gswap (&A->ur, &A->lr);
		}
		pushg (2);
	}

/* If we didn't do any shifting, make sure *xx is the larger value */

	else {
		if (*xx != x) gswap (xx, yy);
	}

/* All done */

	return (TRUE);
}
