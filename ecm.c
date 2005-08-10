/*
   is sieve tuned properly?  could be even faster.  use assembly?
	do sieve once and store in a file?
   handle B1 and/or B2 above 2^32
   Implement an FFT stage 2
*/

/**************************************************************
 *
 *	ecm.c
 *
 *	ECM and P-1 factoring program
 *
 *	Original author:  Richard Crandall - www.perfsci.com
 *	Adapted to Mersenne numbers and optimized by George Woltman
 *	Further optimizations from Paul Zimmerman's GMP-ECM program
 *	Other important ideas courtesy of Peter Montgomery.
 *
 *	c. 1997 Perfectly Scientific, Inc.
 *	c. 1998-2005 Just For Fun Software, Inc.
 *	All Rights Reserved.
 *
 *************************************************************/

/* Global variables */

unsigned long D;		/* Stage 2 loop size */
unsigned long E;		/* Suyama's power in stage 2 */
giant	N = NULL;		/* Number being factored */
giant	FAC = NULL;		/* Found factor */

gwnum	Ad4 = NULL;
gwnum	*nQx = NULL;		/* Array of data used in stage 2 */
gwnum	*eQx = NULL;		/* Array of data used in stage 2 of P-1 */
char	*pairings = NULL;	/* Bits used in determining if primes pair */

#define POOL_3MULT	2	/* Use an algorithm that takes 3 multiplies */
#define POOL_N_SQUARED	4	/* Use an O(N^2) multiplies algorithm */
int	pool_type;		/* Algorithm type to use */
unsigned int pool_count = 0;	/* Count of pooled normalizes */
int	pool_ffted = 0;		/* TRUE if pooled values were pre-ffted */
gwnum	pool_modinv_value = NULL;/* Value we will eventually do a modinv on */
gwnum	*pool_values = NULL;	/* Array of values to normalize */
gwnum	*poolz_values = NULL;	/* Array of z values we are normalize */

gwnum Qprevmx, Qprevmz, Qmx, Qmz, Q2Dxplus1, Q2Dxminus1, *mQx;
unsigned int mQx_count;

unsigned long modinv_count = 0;
int	TWO_FFT_STAGE2 = 0;	/* Type of ECM stage 2 to execute */

int	QA_IN_PROGRESS = FALSE;
int	QA_TYPE = 0;
int	QA_SAVE_FILES = 0;
#define QA_SAVE_TEST 33

/* Bit manipulation macros */

#define bitset(a,i)	{ a[(i) >> 3] |= (1 << ((i) & 7)); }
#define bitclr(a,i)	{ a[(i) >> 3] &= ~(1 << ((i) & 7)); }
#define bittst(a,i)	(a[(i) >> 3] & (1 << ((i) & 7)))

/* Perform setup functions.  This includes decding how big an FFT to */
/* use, allocating memory, calling the FFT setup code, etc. */

int ecm_setup1 (
	double	k,
	unsigned long b,
	unsigned long n,
	signed long c)
{
	int	res;

/* Setup the assembly code */

	res = gwsetup (k, b, n, c, 0);
	if (res) return (res);
	pool_count = 0;

/* A kludge so that the error checking code is not as strict. */

	MAXDIFF *= IniGetInt (INI_FILE, "MaxDiffMultiplier", 1);
	return (0);
}

/* Perform setup functions, part 2. */

void ecm_setup2 (void)
{
	unsigned long i, max;

/* Allocate more memory.  D/3 is enough for the nQx values and */
/* we need an additional D/3 or E*2 values for pooling in case we */
/* are using the POOL_3MULT algorithm */

	max = (D/3 > E*2) ? D/3 : E*2;
	gw_set_max_allocs (D/3 + max + 20);
	nQx = (gwnum *) malloc ((D/2) * sizeof (gwnum));
	for (i = 0; i < D/2; i++) nQx[i] = NULL;
	pool_values = (gwnum *) malloc (max * sizeof (gwnum));
	poolz_values = (gwnum *) malloc (max * sizeof (gwnum));
	mQx = (gwnum *) malloc (E * sizeof (gwnum));
	pairings = (char *) malloc ((D + 15) >> 4);
}

/* Perform cleanup functions. */

void ecm_cleanup (void)
{
	free (N);
	free (nQx);
	free (pool_values);
	free (poolz_values);
	free (mQx);
	free (pairings);
	gwdone ();
}


/* Use a simple sieve to find prime numbers */

#define MAX_PRIMES	6600
static	unsigned int *primes = NULL;
static	struct sieve_info {
	unsigned long first_number;
	unsigned int bit_number;
	unsigned int num_primes;
	unsigned long start;
	char	array[4096];
} si = {0};

/* Fill up the sieve array */

void fill_sieve (void)
{
	unsigned int i, fmax;

/* Determine the first bit to clear */

	fmax = (unsigned int)
		sqrt ((double) (si.first_number + sizeof (si.array) * 8 * 2));
	for (i = si.num_primes; i < MAX_PRIMES * 2; i += 2) {
		unsigned long f, r, bit;
		f = primes[i];
		if (f > fmax) break;
		if (si.first_number == 3) {
			bit = (f * f - 3) >> 1;
		} else {
			r = (unsigned long) (si.first_number % f);
			if (r == 0) bit = 0;
			else if (r & 1) bit = (f - r) / 2;
			else bit = (f + f - r) / 2;
			if (f == si.first_number + 2 * bit) bit += f;
		}
		primes[i+1] = bit;
	}
	si.num_primes = i;

/* Fill the sieve wth ones, then zero out the composites */

	memset (si.array, 0xFF, sizeof (si.array));
	for (i = 0; i < si.num_primes; i += 2) {
		unsigned int f, bit;
		f = primes[i];
		for (bit = primes[i+1]; bit < sizeof (si.array) * 8; bit += f)
			bitclr (si.array, bit);
		primes[i+1] = bit - sizeof (si.array) * 8;
	}
	si.bit_number = 0;
}

/* Start sieve by allocate a sieve info structure */

void start_sieve (
	unsigned long start)
{
	unsigned int i;

/* Remember starting point (in case its 2) and make real start odd */

	if (start < 2) start = 2;
	si.start = start;
	start |= 1;

/* See if we can just reuse the existing sieve */

	if (si.first_number &&
	    start >= si.first_number &&
	    start < si.first_number + sizeof (si.array) * 8 * 2) {
		si.bit_number = (start - si.first_number) / 2;
		return;
	}

/* Initialize sieve */

	if (primes == NULL) {
		unsigned int f;
		primes = (unsigned int *)
			malloc (MAX_PRIMES * 2 * sizeof (unsigned int));
		for (i = 0, f = 3; i < MAX_PRIMES * 2; f += 2)
			if (isPrime (f)) primes[i] = f, i += 2;
	}

	si.first_number = start;
	si.num_primes = 0;
	fill_sieve ();
}

/* Return next prime from the sieve */

unsigned long sieve (void)
{
	if (si.start == 2) {
		si.start = 3;
		return (2);
	}
	for ( ; ; ) {
		unsigned int bit;
		if (si.bit_number == sizeof (si.array) * 8) {
			si.first_number += 2 * sizeof (si.array) * 8;
			fill_sieve ();
		}
		bit = si.bit_number++;
		if (bittst (si.array, bit))
			return (si.first_number + 2 * bit);
	}
}

/* Simple routine to determine if two numbers are relatively prime */

int relatively_prime (
	unsigned long i,
	unsigned long D)
{
	unsigned long f;
	for (f = 3; f * f <= i; f += 2) {
		if (i % f != 0) continue;
		if (D % f == 0) return (FALSE);
		do {
			i = i / f;
		} while (i % f == 0);
	}
	return (i == 1 || D % i != 0);
}

/**************************************************************
 *
 *	Functions
 *
 **************************************************************/

/* computes 2P=(x2:z2) from P=(x1:z1), uses the global variables Ad4 */

void ell_dbl (
	gwnum	x1,
	gwnum	z1,
	gwnum	x2,
	gwnum	z2)
{					/* 10 FFTs */
	gwnum	t1, t3;
	t1 = gwalloc ();
	t3 = gwalloc ();
	gwaddsub4 (x1, z1, t1, x2);
	gwsquare (t1);			/* t1 = (x1 + z1)^2 */
	gwsquare (x2);			/* t2 = (x1 - z1)^2 (store in x2) */
	gwsub3 (t1, x2, t3);		/* t3 = t1 - t2 = 4 * x1 * z1 */
	gwfft (t3, t3);
	gwfft (x2, x2);
	gwfftadd3 (t3, x2, t1);		/* Compute the fft of t1! */
	gwfftfftmul (Ad4, x2, x2);	/* x2 = t2 * Ad4 */
	gwfft (x2, x2);
	gwfftadd3 (x2, t3, z2);		/* z2 = (t2 * Ad4 + t3) */
	gwfftfftmul (t3, z2, z2);	/* z2 = z2 * t3 */
	gwfftfftmul (t1, x2, x2);	/* x2 = x2 * t1 */
	gwfree (t1);
	gwfree (t3);
}

/* adds Q=(x2:z2) and R=(x1:z1) and puts the result in (x3:z3),
   Assumes that Q-R=P or R-Q=P where P=(xdiff:zdiff). */

#ifdef ELL_ADD_USED
void ell_add (
	gwnum 	x1,
	gwnum 	z1,
	gwnum 	x2,
	gwnum 	z2,
	gwnum	xdiff,
	gwnum	zdiff,
	gwnum	x3,
	gwnum	z3)
{					/* 16 FFTs */
	gwnum	t1, t2, t3;
	t1 = gwalloc ();
	t2 = gwalloc ();
	t3 = gwalloc ();
	gwaddsub4 (x1, z1, t1, t2);	/* t1 = (x1 + z1)(x2 - z2) */
					/* t2 = (x1 - z1)(x2 + z2) */
	gwsub3 (x2, z2, t3);
	gwmul (t3, t1);
	gwadd3 (x2, z2, t3);
	gwmul (t3, t2);
	gwaddsub (t2, t1);		/* x3 = (t2 + t1)^2 * zdiff */
	gwsquare (t2);
	gwmul (zdiff, t2);
	gwsquare (t1);			/* z3 = (t2 - t1)^2 * xdiff */
	gwmul (xdiff, t1);
	gwcopy (t2, x3);
	gwcopy (t1, z3);
	gwfree (t1);
	gwfree (t2);
	gwfree (t3);
}
#endif

/* Like ell_add except that x1, z1, xdiff, and zdiff have been FFTed */
/* NOTE: x2 and z2 represent the FFTs of (x2+z2) and (x2-z2) respectively. */

void ell_add_special (
	gwnum 	x1,
	gwnum 	z1,
	gwnum 	x2,
	gwnum 	z2,
	gwnum	xdiff,
	gwnum	zdiff,
	gwnum	x3,
	gwnum	z3)
{				/* 10 FFTs */
	gwnum	t1, t2;
	t1 = gwalloc ();
	t2 = gwalloc ();
	gwfftaddsub4 (x1, z1, t1, t2);	/* t1 = (x1 + z1)(x2 - z2) */
					/* t2 = (x1 - z1)(x2 + z2) */
	gwfftfftmul (z2, t1, t1);
	gwfftfftmul (x2, t2, t2);
	gwaddsub (t2, t1);		/* x3 = (t2 + t1)^2 * zdiff */
	gwsquare (t2);
	gwfftmul (zdiff, t2);
	gwsquare (t1);			/* z3 = (t2 - t1)^2 * xdiff */
	gwfftmul (xdiff, t1);
	gwcopy (t2, x3);
	gwcopy (t1, z3);
	gwfree (t1);
	gwfree (t2);
}

/* This routine is called prior to a series of many ell_add_fft and */
/* ell_dbl_fft calls.  The sequence ends by calling ell_add_fft_last. */
/* Note: We used to simply just FFT x1 and z1.  However, convolution error */
/* in computing (x1+z1)^2 and the like was too great.  Instead, we now */
/* save the FFTs of (x1+z1) and (x1-z1).  The multiplication by xdiff */
/* and zdiff is now more complicated, but convolution errors are reduced */
/* since only one argument of any multiply will involve a value that is */
/* the sum of two FFTs rather than computing a properly normalized sum */
/* and then taking the FFT. */

void ell_begin_fft (
	gwnum	x1,
	gwnum	z1,
	gwnum	x2,
	gwnum	z2)
{
	gwaddsub4 (x1, z1, x2, z2);	/* x2 = x1 + z1, z2 = x1 - z1 */
	gwfft (x2, x2);
	gwfft (z2, z2);
}

/* Like ell_dbl, but the input arguments are FFTs of x1=x1+z1, z1=x1-z1 */
/* The output arguments are also FFTs of x2=x2+z2, z2=x2-z2 */

void ell_dbl_fft (
	gwnum	x1,
	gwnum	z1,
	gwnum	x2,
	gwnum	z2)
{					/* 10 FFTs, 4 adds */
	gwnum	t1, t3;
	t1 = gwalloc ();
	t3 = gwalloc ();
	gwfftfftmul (x1, x1, t1);	/* t1 = (x1 + z1)^2 */
	gwfftfftmul (z1, z1, x2);	/* t2 = (x1 - z1)^2 (store in x2) */
	gwsub3 (t1, x2, t3);		/* t3 = t1 - t2 = 4 * x1 * z1 */
	gwfft (t3, t3);
	gwfft (x2, x2);
	gwfftadd3 (t3, x2, t1);		/* Compute fft of t1! */
	gwfftfftmul (Ad4, x2, x2);	/* x2 = t2 * Ad4 */
	gwfft (x2, x2);
	gwfftadd3 (x2, t3, z2);		/* z2 = (t2 * Ad4 + t3) * t3 */
	gwfftfftmul (t3, z2, z2);
	gwfftfftmul (t1, x2, x2);	/* x2 = x2 * t1 */
	gwaddsub (x2, z2);		/* x2 = x2 + z2, z2 = x2 - z2 */
	gwfft (x2, x2);
	gwfft (z2, z2);
	gwfree (t1);
	gwfree (t3);
}

/* Like ell_add but input arguments are FFTs of x1=x1+z1, z1=x1-z1, */
/* x2=x2+z2, z2=x2-z2, xdiff=xdiff+zdiff, zdiff=xdiff-zdiff. */
/* The output arguments are also FFTs of x3=x3+z3, z3=x3-z3 */

void ell_add_fft (
	gwnum 	x1,
	gwnum 	z1,
	gwnum 	x2,
	gwnum 	z2,
	gwnum 	xdiff,
	gwnum 	zdiff,
	gwnum 	x3,
	gwnum 	z3)
{				/* 12 FFTs, 6 adds */
	gwnum	t1, t2;
	t1 = gwalloc ();
	t2 = gwalloc ();
	gwfftfftmul (x1, z2, t1);/* t1 = (x1 + z1)(x2 - z2) */
	gwfftfftmul (x2, z1, t2);/* t2 = (x1 - z1)(x2 + z2) */
	gwaddsub (t2, t1);
	gwsquare (t2);		/* t2 = (t2 + t1)^2 (will become x3) */
	gwsquare (t1);		/* t1 = (t2 - t1)^2 (will become z3) */
	gwfftaddsub4 (xdiff, zdiff, x3, z3);
				/* x3 = xdiff = (xdiff + zdiff) */
				/* z3 = zdiff = (xdiff - zdiff) */
	gwfftmul (z3, t2);	/* t2 = t2 * zdiff (new x3) */
	gwfftmul (x3, t1);	/* t1 = t1 * xdiff (new z3) */
	gwaddsub (t2, t1);	/* t2 = x3 + z3, t1 = x3 - z3 */
	gwfft (t2, x3);
	gwfft (t1, z3);
	gwfree (t1);
	gwfree (t2);
}

/* Like ell_add_fft but output arguments are not FFTed. */

void ell_add_fft_last (
	gwnum 	x1,
	gwnum 	z1,
	gwnum 	x2,
	gwnum 	z2,
	gwnum 	xdiff,
	gwnum 	zdiff,
	gwnum 	x3,
	gwnum 	z3)
{				/* 10 FFTs, 6 adds */
	gwnum	t1, t2;
	t1 = gwalloc ();
	t2 = gwalloc ();
	gwfftfftmul (x1, z2, t1);/* t1 = (x1 + z1)(x2 - z2) */
	gwfftfftmul (x2, z1, t2);/* t2 = (x1 - z1)(x2 + z2) */
	if (xdiff != x3) {
		gwaddsub4 (t2, t1, x3, z3);
		gwsquare (x3);		/* x3 = (t2 + t1)^2 */
		gwsquare (z3);		/* z3 = (t2 - t1)^2 */
		gwfftaddsub4 (xdiff, zdiff, t1, t2);
				/* t1 = xdiff = (xdiff + zdiff) */
				/* t2 = zdiff = (xdiff - zdiff) */
		gwfftmul (t2, x3);	/* x3 = x3 * zdiff */
		gwfftmul (t1, z3);	/* z3 = z3 * xdiff */
	} else {
		gwaddsub (t2, t1);
		gwsquare (t2); gwfft (t2, t2);
		gwsquare (t1); gwfft (t1, t1);
		gwfftaddsub4 (xdiff, zdiff, z3, x3);
		gwfftfftmul (t2, x3, x3);
		gwfftfftmul (t1, z3, z3);
	}
	gwfree (t1);
	gwfree (t2);
}

/* Perform an elliptic multiply using an algorithm developed by */
/* Peter Montgomery.  Basically, we try to find a near optimal */
/* Lucas chain of additions that generates the number we are */
/* multiplying by.  This minimizes the number of calls to ell_dbl */
/* and ell_add. */

/* The costing function assigns an ell_dbl call a cost of 12 and */
/* an ell_add call a cost of 12.  This cost estimates the number */
/* of forward and inverse transforms performed. */

#define swap(a,b)	{t=a;a=b;b=t;}

unsigned long lucas_cost (
	unsigned long n,
	double	v)
{
	unsigned long c, d, e, t, dmod3, emod3;

	c = 0;
	while (n != 1) {
	    d = (unsigned long) (n/v+0.5); e = n - d;
	    d = d - e;

	    c += 12;

	    while (d != e) {
		if (d < e) {
			swap (d,e);
		}
		if (d <= e + (e >> 2)) {
			if ((dmod3 = d%3) == 3 - (emod3 = e%3)) {
				t = d;
				d = (d+d-e)/3;
				e = (e+e-t)/3;
				c += 36;
				continue;
			}
			if (dmod3 == emod3 && (d&1) == (e&1)) {
				d = (d-e) >> 1;
				c += 22;
				continue;
			}
		}
		if (d <= (e << 2)) {
			d = d-e;
			c += 12;
		} else if ((d&1) == (e&1)) {
			d = (d-e) >> 1;
			c += 22;
		} else if ((d&1) == 0) {
			d = d >> 1;
			c += 22;
		} else if ((dmod3 = d%3) == 0) {
			d = d/3-e;
			c += 46;
		} else if (dmod3 == 3 - (emod3 = e%3)) {
			d = (d-e-e)/3;
			c += 46;
		} else if (dmod3 == emod3) {
			d = (d-e)/3;
			c += 46;
		} else {
			e = e >> 1;
			c += 22;
		}
	    }
	    c += 12;
	    n = d;
	}

	return (c);
}

void lucas_mul (
	gwnum	xx,
	gwnum	zz,
	unsigned long n,
	double	v)
{
	unsigned long d, e, t, dmod3, emod3;
	gwnum	xA, zA, xB, zB, xC, zC, xs, zs, xt, zt;

	xA = gwalloc ();
	zA = gwalloc ();
	xB = gwalloc ();
	zB = gwalloc ();
	xC = gwalloc ();
	zC = gwalloc ();
	xs = xx;
	zs = zz;
	xt = gwalloc ();
	zt = gwalloc ();

	while (n != 1) {
	    ell_begin_fft (xx, zz, xA, zA);			/* A */
	    ell_dbl_fft (xA, zA, xB, zB);			/* B = 2*A */
	    gwcopy (xA, xC); gwcopy (zA, zC);			/* C = A */

	    d = (unsigned long) (n/v+0.5); e = n - d;
	    d = d - e;

	    while (d != e) {
		if (d < e) {
			swap (d, e);
			gwswap (xA, xB); gwswap (zA, zB);
		}
		if (d <= e + (e >> 2)) {
			if ((dmod3 = d%3) == 3 - (emod3 = e%3)) {
				ell_add_fft (xA, zA, xB, zB, xC, zC, xs, zs);/* S = A+B */
				ell_add_fft (xA, zA, xs, zs, xB, zB, xt, zt);/* T = A+S */
				ell_add_fft (xs, zs, xB, zB, xA, zA, xB, zB);/* B = B+S */
				gwswap (xt, xA); gwswap (zt, zA);/* A = T */
				t = d;
				d = (d+d-e)/3;
				e = (e+e-t)/3;
				continue;
			}
			if (dmod3 == emod3 && (d&1) == (e&1)) {
				ell_add_fft (xA, zA, xB, zB, xC, zC, xB, zB);/* B = A+B */
				ell_dbl_fft (xA, zA, xA, zA);	/* A = 2*A */
				d = (d-e) >> 1;
				continue;
			}
		}
		if (d <= (e << 2)) {
			ell_add_fft (xA, zA, xB, zB, xC, zC, xC, zC);/* B = A+B */
			gwswap (xB, xC); gwswap (zB, zC);	/* C = B */
			d = d-e;
		} else if ((d&1) == (e&1)) {
			ell_add_fft (xA, zA, xB, zB, xC, zC, xB, zB);/* B = A+B */
			ell_dbl_fft (xA, zA, xA, zA);		/* A = 2*A */
			d = (d-e) >> 1;
		} else if ((d&1) == 0) {
			ell_add_fft (xA, zA, xC, zC, xB, zB, xC, zC);/* C = A+C */
			ell_dbl_fft (xA, zA, xA, zA);		/* A = 2*A */
			d = d >> 1;
		} else if ((dmod3 = d%3) == 0) {
			ell_dbl_fft (xA, zA, xs, zs);		/* S = 2*A */
			ell_add_fft (xA, zA, xB, zB, xC, zC, xt, zt);/* T = A+B */
			ell_add_fft (xs, zs, xA, zA, xA, zA, xA, zA);/* A = S+A */
			ell_add_fft (xs, zs, xt, zt, xC, zC, xC, zC);/* B = S+T */
			gwswap (xB, xC); gwswap (zB, zC);	/* C = B */
			d = d/3-e;
		} else if (dmod3 == 3 - (emod3 = e%3)) {
			ell_add_fft (xA, zA, xB, zB, xC, zC, xs, zs);/* S = A+B */
			ell_add_fft (xA, zA, xs, zs, xB, zB, xB, zB);/* B = A+S */
			ell_dbl_fft (xA, zA, xs, zs);		/* S = 2*A */
			ell_add_fft (xs, zs, xA, zA, xA, zA, xA, zA);/* A = S+A */
			d = (d-e-e)/3;
		} else if (dmod3 == emod3) {
			ell_add_fft (xA, zA, xB, zB, xC, zC, xt, zt);/* T = A+B */
			ell_add_fft (xA, zA, xC, zC, xB, zB, xC, zC);/* C = A+C */
			gwswap (xt, xB); gwswap (zt, zB);	/* B = T */
			ell_dbl_fft (xA, zA, xs, zs);		/* S = 2*A */
			ell_add_fft (xs, zs, xA, zA, xA, zA, xA, zA);/* A = S+A */
			d = (d-e)/3;
		} else {
			ell_add_fft (xB, zB, xC, zC, xA, zA, xC, zC);/* C = C-B */
			ell_dbl_fft (xB, zB, xB, zB);		/* B = 2*B */
			e = e >> 1;
		}
	    }

	    ell_add_fft_last (xB, zB, xA, zA, xC, zC, xx, zz);	/* A = A+B */

	    n = d;
	}
	gwfree (xA);
	gwfree (zA);
	gwfree (xB);
	gwfree (zB);
	gwfree (xC);
	gwfree (zC);
	gwfree (xt);
	gwfree (zt);
}

/* Multiplies the point (xx,zz) by n using a combination */
/* of ell_dbl and ell_add calls */

void bin_ell_mul (
	gwnum	xx,
	gwnum	zz,
	unsigned long n)
{
	unsigned long c, zeros;
	gwnum	xorg, zorg, xs, zs;

	xorg = gwalloc ();
	zorg = gwalloc ();
	xs = gwalloc ();
	zs = gwalloc ();

	for (zeros = 0; (n & 1) == 0; zeros++) n >>= 1;

	if (n > 1) {
		ell_begin_fft (xx, zz, xorg, zorg);

		c = (unsigned long)(1<<31);
		while ((c&n) == 0) c >>= 1;
		c >>= 1;

		/* If the second bit is zero, we can save one ell_dbl call */

		if (c&n) {
			gwcopy (xorg, xx); gwcopy (zorg, zz);
			ell_dbl_fft (xx, zz, xs, zs);
		} else {
			ell_dbl_fft (xorg, zorg, xx, zz);
			ell_add_fft (xorg, zorg, xx, zz, xorg, zorg, xs, zs);
			c >>= 1;
		}

		/* Do the rest of the bits */

		do {
			if (c&n) {
				if (c == 1) {
					ell_add_fft_last (xs, zs, xx, zz,
							  xorg, zorg, xx, zz);
				} else {
					ell_add_fft (xs, zs, xx, zz,
						     xorg, zorg, xx, zz);
					ell_dbl_fft (xs, zs, xs, zs);
				}
			} else {
				ell_add_fft (xx, zz, xs, zs,
					     xorg, zorg, xs, zs);
				ell_dbl_fft (xx, zz, xx, zz);
			}
			c >>= 1;
		} while (c);
	}

	gwfree (xorg); 
	gwfree (zorg); 
	gwfree (xs); 
	gwfree (zs); 

	while (zeros--) ell_dbl (xx, zz, xx, zz);
}

/* Try a series of Lucas chains to find the cheapest. */
/* First try v = (1+sqrt(5))/2, then (2+v)/(1+v), then (3+2*v)/(2+v), */
/* then (5+3*v)/(3+2*v), etc.  Finally, execute the cheapest. */
/* This is much faster than bin_ell_mul, but uses more memory. */

void ell_mul (
	gwnum	xx,
	gwnum	zz,
	unsigned long n)
{
	unsigned long zeros;

	for (zeros = 0; (n & 1) == 0; zeros++) n >>= 1;

	if (n > 1) {
		unsigned long c, min;
		double	minv;

		min = lucas_cost (n, minv = 1.61803398875);/*v=(1+sqrt(5))/2*/

		c = lucas_cost (n, 1.38196601125);	/*(2+v)/(1+v)*/
		if (c < min) min = c, minv = 1.38196601125;

		c = lucas_cost (n, 1.72360679775);	/*(3+2*v)/(2+v)*/
		if (c < min) min = c, minv = 1.72360679775;

		c = lucas_cost (n, 1.580178728295);	/*(5+3*v)/(3+2*v)*/
		if (c < min) min = c, minv = 1.580178728295;

		c = lucas_cost (n, 1.632839806089);	/*(8+5*v)/(5+3*v)*/
		if (c < min) min = c, minv = 1.632839806089;

		c = lucas_cost (n, 1.612429949509);	/*(13+8*v)/(8+5*v)*/
		if (c < min) min = c, minv = 1.612429949509;

		c = lucas_cost (n, 1.620181980807);	/*(21+13*v)/(13+8*v)*/
		if (c < min) min = c, minv = 1.620181980807;

		c = lucas_cost (n, 1.617214616534);	/*(34+21*v)/(21+13*v)*/
		if (c < min) min = c, minv = 1.617214616534;

		c = lucas_cost (n, 1.618347119656);	/*(55+34*v)/(34+21*v)*/
		if (c < min) min = c, minv = 1.618347119656;

		c = lucas_cost (n, 1.617914406529);	/*(89+55*v)/(55+34*v)*/
		if (c < min) min = c, minv = 1.617914406529;

		lucas_mul (xx, zz, n, minv);
	}
	while (zeros--) ell_dbl (xx, zz, xx, zz);
}

/* Test if factor divides N, return TRUE if it does */

int testFactor (
	double	k,		/* K in K*B^N+C */
	unsigned long b,	/* B in K*B^N+C */
	unsigned long n,	/* N in K*B^N+C */
	signed long c,		/* C in K*B^N+C */
	giant	f)
{
	giant	tmp;
	int	divides_ok;

	tmp = popg (f->sign + 5);	/* Allow room for mul by KARG */
	itog (b, tmp);
	powermod (tmp, n, f);
	dblmulg (k, tmp);
	iaddg (c, tmp);
	modg (f, tmp);
	divides_ok = isZero (tmp);
	pushg (1);
	return (divides_ok);
}

/* Set N, the number we are trying to factor */

int setN (
	double	k,		/* K in K*B^N+C */
	unsigned long b,	/* B in K*B^N+C */
	unsigned long n,	/* N in K*B^N+C */
	signed long c)		/* C in K*B^N+C */
{
	unsigned long bits, p;
	FILE	*fd;
	char	buf[2500];

/* Create the binary representation of the number we are factoring */
/* Allocate 10 extra words to handle any possible k value. */

	bits = (unsigned long) (n * log ((double) b) / log ((double) 2.0));
	N = newgiant ((bits >> 4) + 10);
	ultog (b, N);
	power (N, n);
	dblmulg (k, N);
	iaddg (c, N);

/* Ignore file of known factors when QAing */

	if (QA_IN_PROGRESS) return (TRUE);

/* Open file of known factors */

	if (k != 1.0 || b != 2 || abs(c) != 1) return (TRUE);
	fd = fopen (c == 1 ? "lowp.txt" : "lowm.txt", "r");
	if (fd == NULL) return (TRUE);

/* Loop until the entire file is processed */
/* We are looking for lines of the form: "M( 2843 )C: 142151" */

	while (fscanf (fd, "%s", buf) != EOF) {
		giant	f;

		if (buf[0] != 'M' && buf[0] != 'P') continue;
		fscanf (fd, "%ld", &p);
		if (p > n) break;
		if (p < n) continue;
		fscanf (fd, "%s", buf);
		if (buf[1] != 'C') continue;

/* Get the factor */

		fscanf (fd, "%s", buf);
		f = newgiant (strlen (buf));
		ctog (buf, f);

/* Divide N by factor - but first verify the factor */

		if (!testFactor (k, b, n, c, f)) {
			free (f);
			fclose (fd);
			OutputBoth (c == 1 ?
				"Factor in lowp.txt does not divide 2^P+1\n" :
				"Factor in lowm.txt does not divide 2^P-1\n");
			return (FALSE);
		}
		divg (f, N);
		free (f);
	}

/* Close file and return */

	fclose (fd);
	return (TRUE);
}

/* Do a GCD of the input value and N to see if a factor was found. */
/* The GCD is returned in FAC iff a factor is found. */
/* Returns TRUE if GCD completed, FALSE if it was interrupted */

int gcd (
	gwnum	gg)
{
	giant	v, save;
	int	retval;

/* Convert input number to binary */

	v = popg ((PARG >> 5) + 10);
	save = popg ((PARG >> 5) + 10);
	gwtogiant (gg, v);
	gtog (v, save);

/* Do the GCD and let the gcdg code use gwnum gg's memory. */

	gwfree_temporarily (gg);
	retval = gcdgi (N, v);
	gwrealloc_temporarily (gg);

/* Restore the input argument */

	gianttogw (save, gg);

/* If a factor was found, save it in FAC */

	if (retval && ! isone (v) && gcompg (N, v)) {
		FAC = newgiant ((bitlen (v) >> 4) + 1);
		gtog (v, FAC);
	}

/* Cleanup and return */

	pushg (2);
	return (retval);
}

/* Computes the modular inverse of a number */
/* This is done using the extended GCD algorithm */
/* The GCD is returned in FAC.  Function returns FALSE */
/* if it was interrupted by an escape. */

int modinv (
	gwnum b)
{
	giant	v;
	int	retval;

/* Convert input number to binary */

	v = popg ((PARG >> 5) + 10);
	gwtogiant (b, v);

/* Let the invg code use gwnum b's memory. */
/* Compute 1/v mod N */

	gwfree_temporarily (b);
	retval = invgi (N, v);
	gwrealloc_temporarily (b);
	if (!retval) return (FALSE);

/* If a factor was found, save it in FAC */

	if (v->sign < 0) {
		negg (v);
		FAC = newgiant ((bitlen (v) >> 4) + 1);
		gtog (v, FAC);
	}

/* Otherwise, convert the inverse to FFT-ready form */

	else {
		gianttogw (v, b);
	}
	pushg (1);

/* Increment count and return */

	modinv_count++;
	return (TRUE);
}

/* Computes the modular inverse of an array of numbers */
/* Uses extra multiplications to make only one real modinv call */
/* Uses the simple formula 1/a = b * 1/ab, 1/b = a * 1/ab */
/* If we accidentally find a factor, it is returned in U. */
/* Return FALSE if there is a user interrupt. */

int grouped_modinv (
	gwnum	*b,
	unsigned int size,
	gwnum	*tmp)
{
	unsigned int i;
	gwnum	*orig_tmp;

/* Handle group of 1 as a special case */

	if (size == 1) return (modinv (*b));

/* Handle an odd size */

	orig_tmp = tmp;
	if (size & 1) {
		gwswap (b[0], *tmp);
		tmp++;
	}

/* Multiply each pair of numbers */

	for (i = (size & 1); i < size; i += 2) {
		gwfft (b[i], b[i]);
		gwfft (b[i+1], b[i+1]);
		gwfftfftmul (b[i], b[i+1], *tmp);
		tmp++;
	}

/* Recurse */

	if (!grouped_modinv (orig_tmp, (size+1) / 2, tmp)) return (FALSE);
	if (FAC != NULL) return (TRUE);

/* Handle an odd size */

	if (size & 1) {
		gwswap (b[0], *orig_tmp);
		orig_tmp++;
	}

/* Now perform multiplications on each pair to get the modular inverse */

	for (i = (size & 1); i < size; i += 2) {
		gwfft (*orig_tmp, *orig_tmp);
		gwfftfftmul (*orig_tmp, b[i], b[i]);
		gwfftfftmul (*orig_tmp, b[i+1], b[i+1]);
		gwswap (b[i], b[i+1]);
		orig_tmp++;
	}

/* All done, return TRUE */

	return (TRUE);
}

/* Takes a point (a,b) and multiplies it by a value such that b will be one */
/* If we accidentally find a factor it is returned in FAC.  Function returns */
/* FALSE if it was interrupted by an escape. */

int normalize (
	gwnum	a,
	gwnum	b)
{

/* Compute the modular inverse and scale up the first input value */

	if (!modinv (b)) return (FALSE);
	if (FAC != NULL) return (TRUE);
	gwmul (b, a);
	return (TRUE);
}

/* Adds a point (a,b) to the list of numbers that need normalizing. */
/* This is done in such a way as to minimize the amount of memory used. */

/* This is an interesting bit of code with a variety of algorithms */
/* available.  Assuming there are N pairs to normalize, then you can: */
/* 1) Use 3*N memory and use as few as 2 multiplies per pair. */
/* 2) Use 2*N memory and use 3 multiplies per pair. */
/* 3) Use N+log N memory and use O(log N) multiplies */
/* 4) Use N memory and use O(N^2) multiplies. */

void add_to_normalize_pool (
	gwnum	a,
	gwnum	b,
	int	ffted)		/* TRUE if input arguments have been FFTed */
{

/* Switch off the type of pooling we are going to do */

	switch (pool_type) {

/* Implement algorithm 2 above */

	case POOL_3MULT:

/* If this is the first call allocate memory for the gwnum we use */

		if (pool_count == 0) {
			pool_modinv_value = gwalloc ();
			gwcopy (b, pool_modinv_value);
			pool_ffted = ffted;
		}

/* Otherwise, multiply a by the accumulated b values */

		else if (ffted) {
			if (pool_count != 1)
				gwfft (pool_modinv_value, pool_modinv_value);
			gwfftfftmul (pool_modinv_value, a, a);
			poolz_values[pool_count] = gwalloc ();
			gwcopy (b, poolz_values[pool_count]);
			gwfftfftmul (poolz_values[pool_count],
				     pool_modinv_value,
				     pool_modinv_value);
		} else {
			gwmul (pool_modinv_value, a);
			poolz_values[pool_count] = gwalloc ();
			gwfft (b, poolz_values[pool_count]);
			gwfftfftmul (poolz_values[pool_count],
				     pool_modinv_value,
				     pool_modinv_value);
		}

/* Add a to array of values to normalize */

		pool_values[pool_count++] = a;
		break;

/* Implement algorithm 4 above */

	case POOL_N_SQUARED:

/* If this is the first call allocate memory for the gwnum we use */

		if (pool_count == 0) {
			pool_modinv_value = gwalloc ();
			gwcopy (b, pool_modinv_value);
			pool_ffted = ffted;
		}

/* Otherwise, multiply a by the accumulated b values */
/* and multiply all previous a's by this b */

		else if (ffted) {
			unsigned int i;
			if (pool_count != 1)
				gwfft (pool_modinv_value, pool_modinv_value);
			gwfftfftmul (pool_modinv_value, a, a);
			gwfftfftmul (b, pool_modinv_value, pool_modinv_value);
			for (i = 0; i < pool_count; i++)
				if (i == 0 && pool_ffted) {
					gwfftfftmul (b, pool_values[i], pool_values[i]);
					pool_ffted = FALSE;
				} else
					gwfftmul (b, pool_values[i]);
		} else {
			unsigned int i;
			gwnum	tmp;
			gwmul (pool_modinv_value, a);
			tmp = gwalloc ();
			gwfft (b, tmp);
			gwfftfftmul (tmp, pool_modinv_value,pool_modinv_value);
			for (i = 0; i < pool_count; i++)
				if (i == 0 && pool_ffted) {
					gwfftfftmul (tmp, pool_values[i], pool_values[i]);
					pool_ffted = FALSE;
				} else
					gwfftmul (tmp, pool_values[i]);
			gwfree (tmp);
		}

/* Add a to array of values to normalize */

		pool_values[pool_count++] = a;
		break;
	}
}

/* Takes each point from add_to_normalize_pool and normalizes it. */
/* If we accidentally find a factor, it is returned in U. */
/* Return FALSE if there is a user interrupt. */

int normalize_pool (void)
{
	unsigned int i;

/* Compute the modular inverse */

	if (!modinv (pool_modinv_value)) return (FALSE);
	if (FAC != NULL) goto exit;

/* Now invert each value */
/* Switch off the type of pooling we are going to do */

	switch (pool_type) {

/* Implement algorithm 2 above */

	case POOL_3MULT:
		for (i = pool_count-1; ; i--) {
			if (i == 0 && pool_ffted) {
				gwfft (pool_modinv_value, pool_modinv_value);
				gwfftfftmul (pool_modinv_value,
					     pool_values[i],
					     pool_values[i]);
			} else
				gwmul (pool_modinv_value, pool_values[i]);
			if (i == 0) break;
			gwfftfftmul (poolz_values[i],
				     pool_modinv_value, pool_modinv_value);
			gwfree (poolz_values[i]);
		}
		break;

/* Implement algorithm 4 above */

	case POOL_N_SQUARED:
		gwfft (pool_modinv_value, pool_modinv_value);
		for (i = 0; i < pool_count; i++)
			if (i == 0 && pool_ffted) {
				gwfftfftmul (pool_modinv_value,
					     pool_values[i],
					     pool_values[i]);
			} else
				gwfftmul (pool_modinv_value, pool_values[i]);
		break;
	}

/* Cleanup and reinitialize */

exit:	pool_count = 0;
	gwfree (pool_modinv_value);
	return (TRUE);
}


/* Test if N is a probable prime */
/* Compute i^(N-1) mod N for i = 3,5,7 */

int isProbablePrime (void)
{
	int	i, j, len, retval;
	gwnum	t1, t2;
	giant	x;

	if (isone (N)) return (TRUE);

	retval = TRUE;		/* Assume it is a probable prime */
	t1 = gwalloc ();
	len = bitlen (N);
	for (i = 3; retval && i <= 7; i += 2) {
		t2 = gwalloc ();
		dbltogw ((double) 1.0, t1);
		dbltogw ((double) i, t2);
		gwfft (t2, t2);
		for (j = 1; j <= len; j++) {
			gwsquare (t1);
			if (bitval (N, len-j)) gwfftmul (t2, t1);
		}
		gwfree (t2);
		x = popg ((PARG >> 5) + 10);
		gwtogiant (t1, x);
		modg (N, x);
		iaddg (-i, x);
		if (!isZero (x)) retval = FALSE;	/* Not a prime */
		pushg (1);
	}
	gwfree (t1);
	return (retval);
}

/* Print the factor we just found */

int printFactor (
	double	k,		/* K in K*B^N+C */
	unsigned long b,	/* B in K*B^N+C */
	unsigned long n,	/* N in K*B^N+C */
	signed long c)		/* C in K*B^N+C */
{
	int	msglen;
	char	*msg;

	if (!testFactor (k, b, n, c, FAC)) {
		OutputBoth ("ERROR: Factor doesn't divide N!\n");
		ultog (1, FAC);
		return (FALSE);
	}

	msglen = FAC->sign * 10 + 80;
	msg = (char *) malloc (msglen);
	sprintf (msg, "%s has a factor: ", gwmodulo_as_string ());
	gtoc (FAC, msg+strlen(msg), msglen);
	strcat (msg, "\n");
	OutputBoth (msg);
	spoolMessage (PRIMENET_RESULT_MESSAGE, msg);
	free (msg);

	if (QA_IN_PROGRESS) return (TRUE);

	if (n < 10000) {
		divg (FAC, N);
		if (isProbablePrime ()) {
			OutputBoth ("Cofactor is a probable prime!\n");
			updateWorkToDo (k, b, n, c, WORK_FACTOR, 0);
			return (TRUE);
		}
		mulg (FAC, N);
	}
	return (FALSE);
}

/* From R. P. Brent, priv. comm. 1996:
Let s > 5 be a pseudo-random seed (called $\sigma$ in the Tech. Report),

	u/v = (s^2 - 5)/(4s)

Then starting point is (x_1, y_1) where

	x_1 = (u/v)^3
and
	a = (v-u)^3(3u+v)/(4u^3 v) - 2
*/
void choose12 (
	double	k,			/* K in K*B^N+C */
	unsigned long b,		/* B in K*B^N+C */
	unsigned long n,		/* N in K*B^N+C */
	signed long c,			/* C in K*B^N+C */
	gwnum 	x,
	gwnum 	z,
	double 	curve)
{
	gwnum	xs, zs, t1, t2, t3;

	xs = gwalloc ();
	zs = gwalloc ();
	t1 = gwalloc ();
	t2 = gwalloc ();
	t3 = gwalloc ();

again:	dbltogw (curve, zs);
	gwcopy (zs, xs);
	gwsquare (xs);			/* s^2 */
	dbltogw (5.0, t1);
	gwsub (t1, xs);			/* u = s^2 - 5 */
	dbltogw (4.0, t1);
	gwmul (t1, zs);			/* v = 4*s */
	gwcopy (xs, x);
	gwsquare (x);
	gwsafemul (xs, x);		/* x = u^3 */
	gwcopy (zs, z);
	gwsquare (z);
	gwsafemul (zs, z);		/* z = v^3 */

	/* Now for A. */
	gwcopy (zs, t2);
	gwsub (xs, t2);
	gwcopy (t2, t3);
	gwsquare (t2);
	gwmul (t3, t2);		/* (v-u)^3 */
	gwcopy (xs, t3);
	gwadd (t3, t3);
	gwadd (xs, t3);
	gwadd (zs, t3);
	gwmul (t3, t2);		/* An = (v-u)^3 (3u+v) */
	gwcopy (zs, t3);
	gwsafemul (xs, t3);
	gwsquare (xs);
	gwsafemul (xs, t3);
	Ad4 = gwalloc ();
	dbltogw (4.0, Ad4);
	gwmul (t3, Ad4);	/* An/Ad is now A + 2 */
	normalize (Ad4, t2);	/* Normalize so that An is one */
	if (FAC != NULL) {	/* If a factor was found, then normalize */
				/* failed to find a modular inverse. */
		gwfree (Ad4);
		printFactor (k, b, n, c);
		divg (FAC, N);
		free (FAC);
		FAC = NULL;
		goto again;
	}
	dbltogw (4.0, t1);	/* For extra speed, precompute Ad * 4 */
	gwmul (t1, Ad4);
	gwfft (Ad4, Ad4);	/* Even more speed, save FFT of Ad4 */

	gwfree (xs);
	gwfree (zs);
	gwfree (t1);
	gwfree (t2);
	gwfree (t3);
}

/* Print message announcing the start of this curve */

void curve_start_msg (
	unsigned long curve,
	double	 sigma,
	unsigned long B,
	unsigned long C)
{
	char	buf[120];

	sprintf (buf, "curve #%ld", curve);
	title (buf);

	sprintf (buf,
		 "ECM on %s: curve #%ld with s=%.0f, B1=%lu, B2=%lu\n",
		 gwmodulo_as_string (), curve, sigma, B, C);
	OutputStr (buf);
}

/* These routines manage the computing of Q^m in stage 2 */

void mQ_init (gwnum x, unsigned long m, gwnum Q2Dx)
{
	Qprevmx = gwalloc (); Qprevmz = gwalloc ();
	Qmx = gwalloc (); Qmz = gwalloc ();
	gwcopy (x, Qprevmx); dbltogw (1.0, Qprevmz);
	bin_ell_mul (Qprevmx, Qprevmz, m - 4*D);
	gwfft (Qprevmx, Qprevmx); gwfft (Qprevmz, Qprevmz);
	gwcopy (x, Qmx); dbltogw (1.0, Qmz);
	bin_ell_mul (Qmx, Qmz, m - 2*D);
	gwfft (Qmx, Qmx); gwfft (Qmz, Qmz);

	/* There will be no more ell_dbl calls */
	gwfree (Ad4);

	/* Precompute the FFTs of Q2Dx+1 and Q2Dx-1 */
	Q2Dxplus1 = Q2Dx;
	Q2Dxminus1 = gwalloc ();
	gwaddsmall (Q2Dx, -1);
	gwfft (Q2Dx, Q2Dxminus1);
	gwaddsmall (Q2Dx, 2);
	gwfft (Q2Dx, Q2Dxplus1);

	/* Init the arrays used in pooled normalizes of mQx values */
	if (TWO_FFT_STAGE2) {
		unsigned long i;
		for (i = 0; i < E; i++) mQx[i] = gwalloc ();
		mQx_count = 0;
	}
}
int mQ_next (gwnum *retx, gwnum *retz)
{

/* The non-normalized case - simple multiply the last Q^m value */
/* by Q^2D to get the next Q^m value */

	if (!TWO_FFT_STAGE2) {
		ell_add_special (Qmx, Qmz, Q2Dxplus1, Q2Dxminus1,
				 Qprevmx, Qprevmz, Qprevmx, Qprevmz);
		gwswap (Qmx, Qprevmx); gwswap (Qmz, Qprevmz);
		gwfft (Qmx, Qmx); gwfft (Qmz, Qmz);
		*retx = Qmx;
		*retz = Qmz;
		return (TRUE);
	}

/* The normalized case - batch up a bunch of Q^m values and normalize */
/* them.  Then return them one at a time.  Obviously retz need not be */
/* returned since it is always one. */

	if (mQx_count == 0) {
		for ( ; mQx_count < E; mQx_count++) {
			ell_add_special (Qmx, Qmz, Q2Dxplus1, Q2Dxminus1,
					 Qprevmx, Qprevmz, Qprevmx, Qprevmz);
			gwswap (Qmx, Qprevmx); gwswap (Qmz, Qprevmz);
			gwfft (Qmx, Qmx); gwfft (Qmz, Qmz);
			gwcopy (Qmx, mQx[mQx_count]);
			add_to_normalize_pool (mQx[mQx_count], Qmz, 1);
		}
		if (!normalize_pool ()) return (FALSE);
		if (FAC != NULL) return (TRUE);
	}
	*retx = mQx[E-mQx_count];
	gwfft (*retx, *retx);
	mQx_count--;
	return (TRUE);
}
void mQ_term (void)
{
	gwfree (Qprevmx);
	gwfree (Qprevmz);
	gwfree (Qmx);
	gwfree (Qmz);
	gwfree (Q2Dxminus1);
}

/* Choose 4 FFT stage 2 of the 2 FFT stage 2.  Also choose a good */
/* value for D and a good algorithm for normalize_pool. */
/* We try to choose the above such that the number of multiplications */
/* are minimized, yet too much memory isn't used. */

void choose_stage2_plan (
	double	k,
	unsigned long b,
	unsigned long n,
	signed long c,
	unsigned long B,		/* Stage 1 bound */
	unsigned long C,		/* Stage 2 bound */
	int	memory)			/* MB of memory we can use */
{
	unsigned long numvals, d, e;
	unsigned long relprime, beste;
	double	numprimes, numpairings, numsections, numgcdsections;
	double	cost, bestcost, density, gcd_cost;

/* Define constants for the number of transforms for various operations */
/* The GCD cost is based on our timings and an Excel spreadsheet */

#define ELL_ADD_COST		12
#define N_SQUARED_POOL_COST	2
#define MULT3_POOL_COST		7
	gcd_cost = 861.0 * log ((double) n) - 7775.0;
	if (gcd_cost < 100.0) gcd_cost = 100.0;

/* Will there even be a stage 2?  If not, set D and E appropriately */

	if (C <= B) {
		D = 0;
		E = 0;
		return;
	}

/* Figure out how many gwnum values fit in our MB limit */

	numvals = (unsigned long)
			(((double) memory * 1000000.0 -
			  (double) gwmap_to_memused (k, b, n, c)) /
			 (double) gwnum_size (FFTLEN));

/* If memory is really tight, then the 4 FFT - O(n^2) pooling is the */
/* most memory efficient ECM implementation.  Note: D=30 (8 nQx values) */
/* requires 20 gwnums.  The next D value (60) requires 28 gwnums. */

	if (numvals < 28) {
		D = 30;
		E = 0;
		TWO_FFT_STAGE2 = FALSE;
		pool_type = POOL_N_SQUARED;
	}

/* Numprimes below C approximately equals C / (ln(C)-1) */
/* Compute numprimes between B and C */

	numprimes = ceil ((C / (log ((double) C) - 1)) - (B / (log ((double) B) - 1)));

/* Figure out the best value for E when using the O(N^2) pool method */

	beste = (unsigned long) sqrt (gcd_cost / N_SQUARED_POOL_COST) + 1;

/* Loop through various D values choosing the most cost effective one */

	bestcost = 1.0E99;
	d = ((unsigned long) sqrt ((double) (C-B)) / 2310 + 3) * 2310;
	for ( ; ; ) {
		if (d >= 2310) {
			relprime = d / 2310 * 480;
			density = 480.0 / 2310.0;
		} else if (d >= 210) {
			relprime = d / 210 * 48;
			density = 48.0 / 210.0;
		} else {
			relprime = d / 30 * 8;
			density = 8.0 / 30.0;
		}

/* Half the primes are eligible for pairing (numprimes / 2). */
/* The chance that a pairing occurs is numprimes / area.  Area would */
/* normally be C-B.  However, the relprime algorithm makes */
/* our primes much denser than that. */

		numpairings = ceil (
			(numprimes / 2.0 * numprimes / ((C-B) * density)));

/* There will be (C-B)/2D sections */

		numsections = ceil ((double) (C-B) / (double) (d+d));

/* Cost out the 4FFT stage 2 using this D			*/
/* The cost will be:						*/
/*	D/2 ell_add_ffts  + pool_cost (relprime) +		*/
/*	(C-B)/2D ell_add_specials + #primes*4			*/
/* The memory consumed will be:					*/
/*	13 + relprime gwnums if N^2 pooling			*/
/* or	13 + 2*relprime gwnums if 3N pooling			*/
/* Note that MQ_init requires B is at least 4 times D		*/

		if (B >= 4*d && 13 + relprime <= numvals &&
		    (QA_TYPE == 0 || QA_TYPE == 1)) {
			cost = d/2 * ELL_ADD_COST +
			       relprime * relprime * N_SQUARED_POOL_COST +
			       numsections * ELL_ADD_COST +
			       (numprimes - numpairings) * 4;
			if (cost < bestcost) {
				TWO_FFT_STAGE2 = FALSE;
				pool_type = POOL_N_SQUARED;
				D = d;
				E = 0;
				bestcost = cost;
			}
		}
		if (B >= 4*d && 13 + relprime*2 <= numvals &&
		    (QA_TYPE == 0 || QA_TYPE == 2)) {
			cost = d/2 * ELL_ADD_COST +
			       relprime * MULT3_POOL_COST +
			       numsections * ELL_ADD_COST +
			       (numprimes - numpairings) * 4;
			if (cost < bestcost) {
				TWO_FFT_STAGE2 = FALSE;
				pool_type = POOL_3MULT;
				D = d;
				E = 0;
				bestcost = cost;
			}
		}

/* Cost out the 2FFT stage 2 using this D			*/
/* The cost will be:						*/
/*	D/2 ell_add_ffts  + pool_cost (relprime) +		*/
/*	(C-B)/2D ell_add_specials + #primes*2 +			*/
/*	(C-B)/2D/E * pool_cost (e)				*/
/*	(C-B)/2D/E * gcd_cost					*/
/* The memory consumed will be:					*/
/*	13 + relprime gwnums if N^2 pooling			*/
/* or	13 + 2*relprime gwnums if 3N pooling			*/

		if (B >= 4*d && 13 + relprime <= numvals &&
		    (QA_TYPE == 0 || QA_TYPE == 3)) {
			e = numvals - relprime - 13;
			if (e == 0) e = 1;
			if (e > beste) e = beste;
			numgcdsections = ceil (numsections / e);
			cost = d/2 * ELL_ADD_COST +
			       relprime * relprime * N_SQUARED_POOL_COST +
			       numsections * ELL_ADD_COST +
			       (numprimes - numpairings) * 2 +
			       numgcdsections * e * e * N_SQUARED_POOL_COST +
			       numgcdsections * gcd_cost;
			if (cost < bestcost) {
				TWO_FFT_STAGE2 = TRUE;
				pool_type = POOL_N_SQUARED;
				D = d;
				E = e;
				bestcost = cost;
			}
		}
		if (B >= 4*d && 13 + relprime*2 <= numvals &&
		    (QA_TYPE == 0 || QA_TYPE == 4)) {
			e = (numvals - relprime - 13) / 2;
			if (e == 0) e = 1;
			numgcdsections = ceil (numsections / e);
			e = (unsigned long) ceil (numsections / numgcdsections);
			cost = d/2 * ELL_ADD_COST +
			       relprime * MULT3_POOL_COST +
			       numsections * (ELL_ADD_COST + 1) +
			       (numprimes - numpairings) * 2.0 +
			       numgcdsections * e * MULT3_POOL_COST +
			       numgcdsections * gcd_cost;
			if (cost < bestcost) {
				TWO_FFT_STAGE2 = TRUE;
				pool_type = POOL_3MULT;
				D = d;
				E = e;
				bestcost = cost;
			}
		}

/* Cost out the next possible value of D */

		if (d > 2310) d = d - 2310;
		else if (d > 210) d = d - 210;
		else if (d > 30) d = d - 30;
		else break;
	}

/* Print out our selection */

	if (QA_IN_PROGRESS) {
		char	buf[200];
		if (pool_type == POOL_N_SQUARED)
			sprintf (buf, "num temps = %d, D = %d, E = %d\n",
				 numvals, D, E);
		else
			sprintf (buf, "num temps = %d, D = %d\n",
				 numvals, D);
		OutputBoth (buf);
		sprintf (buf, "Stage 2 FFTs = %d, Pool type = %s\n",
			 TWO_FFT_STAGE2 ? 2 : 4,
			 pool_type == POOL_N_SQUARED ? "N squared" :
						       "3 multiplies");
		OutputBoth (buf);
	}
}

/* Routines to read and write a byte array from and to a save file */

int read_array (
	int	fd,
	char	*buf,
	unsigned long len,
	long	*sum)
{
	unsigned long i;
	signed char *sbuf;

	if (_read (fd, buf, len) != len) return (FALSE);
	sbuf = (signed char *) buf;  /* Watcom default for char is unsigned! */
	for (i = 0; i < len; i++) *sum += sbuf[i];
	return (TRUE);
}

int write_array (
	int	fd,
	char	*buf,
	unsigned long len,
	long	*sum)
{
	unsigned long i;
	signed char *sbuf;

	if (len == 0) return (TRUE);
	if (_write (fd, buf, len) != len) return (FALSE);
	sbuf = (signed char *) buf;  /* Watcom default for char is unsigned! */
	for (i = 0; i < len; i++) *sum += sbuf[i];
	return (TRUE);
}

/* Routines to read and write a gwnum from and to a save file */

int read_gwnum (
	int	fd,
	gwnum	g,
	long	*sum)
{
	giant	tmp;
	long	i, len, bytes;

	if (_read (fd, &len, sizeof (long)) != sizeof (long)) return (FALSE);
	if (len == 0) return (FALSE);
	tmp = popg ((PARG >> 5) + 10);
	bytes = len * sizeof (long);
	if (_read (fd, tmp->n, bytes) != bytes) return (FALSE);
	tmp->sign = len;
	*sum += len;
	for (i = 0; i < len; i++) *sum += tmp->n[i];
	gianttogw (tmp, g);
	pushg (1);
	return (TRUE);
}

int write_gwnum (
	int	fd,
	gwnum	g,
	long	*sum)
{
	giant	tmp;
	long	i, len, bytes;

	tmp = popg ((PARG >> 5) + 10);
	gwtogiant (g, tmp);
	len = tmp->sign;
	if (len == 0) return (FALSE);
	if (_write (fd, &len, sizeof (long)) != sizeof (long)) return (FALSE);
	bytes = len * sizeof (long);
	if (_write (fd, tmp->n, bytes) != bytes) return (FALSE);
	*sum += len;
	for (i = 0; i < len; i++) *sum += tmp->n[i];
	pushg (1);
	return (TRUE);
}

/* Routines to read and write longs from and to a save file */

int read_long (
	int	fd,
	unsigned long *val,
	long	*sum)
{
	if (_read (fd, val, sizeof (long)) != sizeof (long)) return (FALSE);
	*sum += *val;
	return (TRUE);
}

int write_long (
	int	fd,
	unsigned long val,
	long	*sum)
{
	if (_write (fd, &val, sizeof (long)) != sizeof (long)) return (FALSE);
	*sum += val;
	return (TRUE);
}

/* Routines to create and read save files for an ECM factoring job */

#define ECM_STAGE1	0
#define ECM_STAGE2	1
void ecm_save (
	char	*filename,
	int	stage,
	unsigned long curve,
	double	sigma,
	unsigned long B,
	unsigned long B_processed,
	unsigned long C_processed,
	gwnum	x,
	gwnum	gg)
{
	char	newfilename[16];
	int	fd;
	unsigned long magicnum, version;
	long	sum = 0, i;

/* If we are allowed to create multiple intermediate files, then */
/* write to a file called yNNNNNNN. */

	strcpy (newfilename, filename);
	if (TWO_BACKUP_FILES) newfilename[0] = 'y';

/* Create the intermediate file */

	fd = _open (newfilename, _O_BINARY|_O_WRONLY|_O_TRUNC|_O_CREAT, CREATE_FILE_ACCESS);
	if (fd < 0) return;

/* Write the file header. */

	magicnum = 0x1a2b3cd4;
	if (_write (fd, &magicnum, sizeof (long)) != sizeof (long))
		goto writeerr;
	version = 1;
	if (_write (fd, &version, sizeof (long)) != sizeof (long))
		goto writeerr;

/* Write the file data */

	i = stage;
	if (! write_long (fd, i, &sum)) goto writeerr;
	if (! write_long (fd, curve, &sum)) goto writeerr;

	if (_write (fd, &sigma, sizeof (double)) != sizeof (double))
		goto writeerr;

	if (! write_long (fd, B, &sum)) goto writeerr;
	if (! write_long (fd, B_processed, &sum)) goto writeerr;
	if (! write_long (fd, C_processed, &sum)) goto writeerr;

/* Write the data values */

	if (! write_gwnum (fd, x, &sum)) goto writeerr;
	if (! write_gwnum (fd, gg, &sum)) goto writeerr;

/* Write the checksum */

	if (_write (fd, &sum, sizeof (long)) != sizeof (long)) goto writeerr;
	_commit (fd);
	_close (fd);

/* Now rename the intermediate files */

	if (TWO_BACKUP_FILES) {
		_unlink (filename);
		rename (newfilename, filename);
	}
	return;

/* An error occured.  Close and delete the current file. */

writeerr:
	_close (fd);
	_unlink (newfilename);
}

/* Read a save file */

int ecm_restore (
	char	*filename,
	int	*stage,
	unsigned long *curve,
	double	*sigma,
	unsigned long *B,
	unsigned long *B_processed,
	unsigned long *C_processed,
	gwnum	x,
	gwnum	gg)
{
	int	fd;
	unsigned long magicnum, version;
	unsigned long tmp;
	long	sum = 0, i;

/* Open the intermediate file */

	fd = _open (filename, _O_BINARY | _O_RDONLY);
	if (fd < 0) goto error;

/* Read the file header */

	if (_read (fd, &magicnum, sizeof (long)) != sizeof (long))
		goto readerr;
	if (magicnum != 0x1a2b3cd4) goto readerr;

	if (_read (fd, &version, sizeof (long)) != sizeof (long)) goto readerr;
	if (version != 1) goto readerr;

/* Read the file data */

	if (! read_long (fd, &tmp, &sum)) goto readerr;
	*stage = (int) tmp;
	if (! read_long (fd, curve, &sum)) goto readerr;

	if (_read (fd, sigma, sizeof (double)) != sizeof (double))
		goto readerr;

	if (! read_long (fd, B, &sum)) goto readerr;
	if (! read_long (fd, B_processed, &sum)) goto readerr;
	if (! read_long (fd, C_processed, &sum)) goto readerr;

/* Read the values */

	if (! read_gwnum (fd, x, &sum)) goto readerr;
	if (! read_gwnum (fd, gg, &sum)) goto readerr;

/* Read and compare the checksum */

	if (_read (fd, &i, sizeof (long)) != sizeof (long)) goto readerr;
	if (i != sum) goto readerr;
	_close (fd);
	return (TRUE);

/* An error occured.  Delete the current intermediate file. */
/* Set stage to -1 to indicate an error. */

readerr:
	OutputStr ("Error reading ECM save file.\n");
	_close (fd);
error:
	_unlink (filename);
	return (FALSE);
}


/**************************************************************
 *
 *	Main ECM Function
 *
 **************************************************************/

int ecm (
	double	k,			/* K in K*B^N+C */
	unsigned long b,		/* B in K*B^N+C */
	unsigned long n,		/* N in K*B^N+C */
	signed long c,			/* C in K*B^N+C */
	unsigned long B,		/* Stage 1 bound */
	unsigned long C_start,		/* Stage 2 starting point (usually equals B) */
	unsigned long C,		/* Stage 2 ending point */
	unsigned long curves_to_do,	/* Number of curves to test */
	double	specific_sigma)		/* Debug only - sigma to test */
{
	unsigned int memory;
	unsigned long sieve_start, SQRT_B, orig_B;
	double	sigma, last_output;
	unsigned long i, j, m, curve, prime;
	char	filename[16], buf[100], fft_desc[100];
	int	res, retval, stage, escaped;
	long	write_time = DISK_WRITE_TIME * 60;
	time_t	start_time, current_time;
	gwnum	x, z, t1, t2, gg;
	gwnum	Q2x, Q2z, Qiminus2x, Qiminus2z, Qdiffx, Qdiffz;

/* Unless a save file indicates otherwise, we are testing our first curve */

	curve = 1;

/* Clear all timers */

restart:
	clear_timers ();

/*#define TIMING1*/
#ifdef TIMING1
if (n == 598) {
	gwnum	x, y;
	ecm_setup1 (k, b, n, c);
	gwsetnormroutine (0, ERRCHK, 0);
	x = gwalloc ();
	y = gwalloc ();
	dbltogw (100.0, x);
	dbltogw (10.0, y);
	gwaddsub (x, y);
	OutputStr ("\n");
}
if (n == 599) {
	giant x, y;
	x = newgiant (1050000);
	y = newgiant (1050000);

	itog (1, x);
	gshiftleft ((1L<<16), x);
	iaddg (-1, x);
	gtog (x, y);
//setzero(y);
	iaddg (-257, y);
	start_timer (0);
//invg (x,y);
	gcdg (x, y);
	end_timer (0);
	OutputStr ("GCD is ");
    	gtoc (y, buf, 100);
	OutputStr (buf);
	OutputStr (", ");
	print_timer (0, TIMER_CLR);
	OutputStr ("\n");
}
if (n == 600) {
int i, j;
giant	x, y, z, a, m;
#define TESTSIZE	200
RDTSC_TIMING = 12;
x = newgiant(2*TESTSIZE); y = newgiant (4*TESTSIZE);
z = newgiant (4*TESTSIZE), a = newgiant (4*TESTSIZE);
m = newgiant (2*TESTSIZE);
srand ((unsigned) time (NULL));
for (i = 0; i < TESTSIZE; i++) {
	x->n[i] = (rand () << 17) + rand ();
	m->n[i] = (rand () << 17) + rand ();
}
x->n[TESTSIZE-1] &= 0x00FFFFFF;
m->n[TESTSIZE-1] &= 0x00FFFFFF;
for (i = TESTSIZE; i >= 40; i--) {
	x->sign = i;
	m->sign = i;
	setmulmode (GRAMMAR_MUL);
	for (j = 0; j < 10; j++) {
		gtog (x, y);
		start_timer (0);
		if (B&1) mulg (m, y);
		else squareg (y);
		end_timer (0);
		if (timers[1] == 0 || timers[1] > timers[0]) timers[1] = timers[0];
		timers[0] = 0;
	}
	setmulmode (KARAT_MUL);
	for (j = 0; j < 10; j++) {
		gtog (x, z);
		start_timer (0);
		if (B&1) mulg (m, z);
		else squareg (z);
		end_timer (0);
		if (timers[2] == 0 || timers[2] > timers[0]) timers[2] = timers[0];
		timers[0] = 0;
	}
	setmulmode (FFT_MUL);
	for (j = 0; j < 10; j++) {
		gtog (x, a);
		start_timer (0);
		if (B&1) mulg (m, a);
		else squareg (a);
		end_timer (0);
		if (timers[3] == 0 || timers[3] > timers[0]) timers[3] = timers[0];
		timers[0] = 0;
	}
	sprintf (buf, "Size: %ld  ", i);
	OutputStr (buf);
	OutputStr ("G: ");
	print_timer (1, TIMER_MS | TIMER_CLR);
	OutputStr (", K: ");
	print_timer (2, TIMER_MS | TIMER_CLR);
	OutputStr (", F: ");
	print_timer (3, TIMER_MS | TIMER_NL | TIMER_CLR);
	if (gcompg (y, z) != 0)
		i--;
	if (gcompg (y, a) != 0)
		i--;
	Sleep (100);
}
return 0;
}
if (n == 601) {
int i, j;
giant	x, a, m;
#define TESTSIZE	260000
RDTSC_TIMING = 12;
x = newgiant(2*TESTSIZE);
a = newgiant (4*TESTSIZE);
m = newgiant (2*TESTSIZE);
srand ((unsigned) time (NULL));
for (i = 0; i < TESTSIZE; i++) {
	x->n[i] = (rand () << 17) + rand ();
	m->n[i] = (rand () << 17) + rand ();
}
x->n[TESTSIZE-1] &= 0x00FFFFFF;
m->n[TESTSIZE-1] &= 0x00FFFFFF;
for (i = 30; i < TESTSIZE/2; i<<=1) {
	x->sign = i;
	m->sign = i;
	setmulmode (FFT_MUL);
	for (j = 0; j < 10; j++) {
		gtog (x, a);
		start_timer (0);
		if (B&1) mulg (m, a);
		else squareg (a);
		end_timer (0);
		if (timers[3] == 0 || timers[3] > timers[0]) timers[3] = timers[0];
		timers[0] = 0;
	}
	sprintf (buf, "Size: %ld  ", i);
	OutputStr (buf);
	OutputStr (", F: ");
	print_timer (3, TIMER_NL | TIMER_CLR | TIMER_MS);
	Sleep (100);
}
return 0;
}
#endif

/* Include timing code when building the debug version of prime95 */
#ifdef GDEBUG
if (n == 600) {
int j, min_test, max_test, test, cnt, NUM_X87_TESTS, NUM_SSE2_TESTS;
gwsetup (1.0, 2, 10000000, -1, 0);
SRCARG = (void *) aligned_malloc (40000000, 4096);
memset (SRCARG, 0, 40000000);
RDTSC_TIMING = 12;
min_test = IniGetInt (INI_FILE, "MinTest", 0);
max_test = IniGetInt (INI_FILE, "MaxTest", min_test);
NUMARG = (long) -1; NUM_X87_TESTS = timeit ();
NUMARG = (long) -2; NUM_SSE2_TESTS = timeit ();
//SetThreadPriority (CURRENT_THREAD, THREAD_PRIORITY_TIME_CRITICAL);
for (j = 0; j < NUM_X87_TESTS + NUM_SSE2_TESTS; j++) {
	cnt = 0;
	test = (j < NUM_X87_TESTS ? j : 1000 + j - NUM_X87_TESTS);
	if (min_test && (test < min_test || test > max_test)) continue;
	if (! (CPU_FLAGS & CPU_SSE2) && test >= 1000) break;
for (i = 1; i <= 50; i++) {
	start_timer (0);
	NUMARG = (long) test;
	timeit ();
	end_timer (0);
	if (timers[1] == 0 || timers[1] > timers[0]) timers[1] = timers[0];
	if (i > 1 && timers[0] < 3.0 * timers[1]) {
		if (timers[0] > 1.5 * timers[1])
			i++;
		timers[2] += timers[0];
		cnt++;
	}
	timers[0] = 0;
}
sprintf (buf, "Test %d: ", test);
OutputBoth (buf);
print_timer (1, TIMER_OUT_BOTH | TIMER_CLR);
timers[2] /= cnt;
OutputBoth (", avg: ");
print_timer (2, TIMER_OUT_BOTH | TIMER_NL | TIMER_CLR);
if (min_test) exit (0);
}
aligned_free (SRCARG);
gwdone ();
return 0;
}
#endif

/* Init filename */

	tempFileName (filename, n);
	strcat (filename, EXTENSION);
	filename[0] = 'e';
	if (c == 1) filename[0]--;

/* Get the current time */

	time (&start_time);

/* Init the random number generator */

	srand ((unsigned) time (NULL));

/* Choose a default value for the second bound if none was specified */
/* MQ_init also requires that B is at least 120 (4 times the minimum D */

	orig_B = B;
	if (C == 0)
		C = (B >= 42900000) ? 4290000000UL : B * 100;
	if (C <= B) C = B;
	else if (B < 120) {
		OutputStr ("Using minimum bound #1 of 120\n");
		B = 120;
	}

/* Set other constants */

	SQRT_B = (unsigned long) sqrt ((double) B);

/* Perform setup functions */

	res = ecm_setup1 (k, b, n, c);
	if (res) {
		sprintf (buf, "Cannot initialize FFT code, errcode=%d\n", res);
		OutputBoth (buf);
		return (FALSE);
	}
	gwsetnormroutine (0, ERRCHK, 0);
	last_output = fft_count = modinv_count = 0;

/* Default is to use up to 24MB of memory */

	memory = avail_mem ();
	if (memory == 0) memory = 24;

/* Choose a good value for D.  One that reduces the number of */
/* multiplications, yet doesn't use too much memory. */

	choose_stage2_plan (k, b, n, c, B, C, memory);

/* Perform the rest of the setup */

	ecm_setup2 ();

/* Compute the number we are factoring */

	if (!setN (k, b, n, c)) {
		ecm_cleanup ();
		return (FALSE);
	}

#ifdef CHECK_LOWM
prime_check ();
#endif

/*#define TIMING*/
#ifdef TIMING
{
unsigned long i, j, limit;
gwnum	n1, n2, n3;
n1 = gwalloc ();
n2 = gwalloc ();
n3 = gwalloc ();
dbltogw (283457283657.0, n2);
for (i = 1; i <= 50; i++) gwsquare (n2); /* gen big random number */
gwcopy (n2, n3);
gwcopy (n2, n1);
if (n < 20000) limit = 100; else limit = 10;
for (i = 1; i <= limit; i++) {
	start_timer (0); gwsquare (n2); end_timer (0);
	start_timer (1); gwmul (n2, n3); end_timer (1);
	start_timer (2); gwfftmul (n2, n3); end_timer (2);
	start_timer (3); normalize (n1, n3); end_timer (3);
	start_timer (4); gwfftfftmul (n2, n2, n2); end_timer (4);
	start_timer (5); gwfft (n2, n2); end_timer (5);
	start_timer (6); gwadd (n2, n2); end_timer (6);
	gwcopy (n1, n2);
}
OutputStr ("100 squares: ");  print_timer (0, TIMER_NL | TIMER_CLR);
OutputStr ("100 muls: ");  print_timer (1, TIMER_NL | TIMER_CLR);
OutputStr ("100 fftmuls: ");  print_timer (2, TIMER_NL | TIMER_CLR);
OutputStr ("100 ffts: ");  print_timer (5, TIMER_NL | TIMER_CLR);
OutputStr ("100 fftfftmuls: ");  print_timer (4, TIMER_NL | TIMER_CLR);
OutputStr ("100 normalizes: ");  print_timer (3, TIMER_NL | TIMER_CLR);
OutputStr ("100 adds: ");  print_timer (6, TIMER_NL | TIMER_CLR);
	start_timer (7);
	start_sieve (2);
	for (i = 0; sieve () < 0xFFFFFFFF; i++);
	end_timer (7);
sprintf (buf, "Sieve: %ld primes found.  ", i);
OutputStr (buf);  print_timer (7, TIMER_NL | TIMER_CLR);
}		
#endif

/*#define TESTING*/
#ifdef TESTING	
{ unsigned long i; long j;
for (i = 0; i < FFTLEN; i++)
	* addr (t1, i) = (i == 0) ? 10000.0 : 0.0;
for (i = 0; i < FFTLEN; i++)
	* addr (t2, i) = (i < 4) ? 1.0 : 0.0;
gwcopy (t1,t3);
//gwsafemul (t1,t2);
gwmul (t1, t2);
gwfftmul (t1, t2);
gwadd (t2, t2);
gwcopy (t1, t2);
gwfftfftmul (t1, t2, t2);
gwsquare (t3);

for (i = 0; i < FFTLEN; i++)
	set_fft_value (t1, i, (i >= 0) ? 1 : 0);
for (i = 0; i < FFTLEN; i++)
	set_fft_value (t2, i, (i >= 0) ? 1 : 0);
for (i = 0; i < FFTLEN; i++)
	set_fft_value (t3, i, (i >= 0) ? 1 : 0);
gwfft(t2, t2);
//gwfft (t3, t3);
gwfft (t1, t1);
gwfftfftmul (t1, t2, t2);
//gwfftfftmul (t3, t3, t3);
gwsquare (t3);
gwadd (t3, t3);
	for (i = 0; i < FFTLEN; i++)
		get_fft_value (t3, i, &j);
gwcopy (t1,t3);
gwsafemul (t1,t2);
	for (i = 0; i < FFTLEN; i++)
		get_fft_value (t2, i, &j);
gwmul (t1, t2);
	for (i = 0; i < FFTLEN; i++)
		get_fft_value (t2, i, &j);
gwfftmul (t1, t2);
	for (i = 0; i < FFTLEN; i++)
		get_fft_value (t2, i, &j);
gwadd (t2, t2);
gwcopy (t1, t2);
gwfftfftmul (t1, t2, t2);
	for (i = 0; i < FFTLEN; i++)
		get_fft_value (t2, i, &j);
gwsquare (t3);
	for (i = 0; i < FFTLEN; i++)
		get_fft_value (t3, i, &j);
}
#endif

/*#define echk*/
#ifdef echk
for ( ; ; n++) {
ecm_setup1 (k, b, n, c);
ecm_setup2 ();
dbltogw (1.0, t1);
dbltogw (283457283654.0, t2);
gwmul (t1, t2);
dbltogw (23876287342234.0, t3);
gwfftmul (t1, t3);
ERRCHK = 1; MAXERR = 0.0;
for (int i = 1; i <= 1000; i++) {
gwsquare (t3);
gwsafemul (t2, t3);
gwsquare (t2);
}
n++;
}
#endif

/* Output a startup message */

	gwfft_description (fft_desc);
	sprintf (buf, "Using %s\n", fft_desc);
	OutputStr (buf);

/* Check for a continuation file */

	if (fileExists (filename)) {
		unsigned long save_B, save_B_processed, save_C_processed;

/* Allocate memory */

		x = gwalloc ();
		z = gwalloc ();
		gg = NULL;

/* Read in the save file */

		if (! ecm_restore (filename, &stage, &curve, &sigma, &save_B,
				   &save_B_processed, &save_C_processed,
				   x, z) ||

/* Handle the case where the save file is no good or we have a save file */
/* with a smaller bound #1 than the bound #1 we are presently working on. */
/* In either case restart the curve (and curve count) from scratch. */

		    B > save_B) {
			gwfree (x);
			gwfree (z);
			goto restart0;
		}

/* Compute Ad4 from sigma */

		curve_start_msg (curve, sigma, B, C);
		t1 = gwalloc ();
		t2 = gwalloc ();
		choose12 (k, b, n, c, t1, t2, sigma);
		gwfree (t1);
		gwfree (t2);

/* Continue in the middle of stage 1 */

		if (stage == ECM_STAGE1) {
			sieve_start = save_B_processed + 1;
			goto restart1;
		}
		
/* Allocate more memory */

		gg = gwalloc ();
		gwswap (z, gg);

/* We've finished stage 1, resume stage 2 */

		if (C > save_C_processed) {
			dbltogw (1.0, z);
			start_sieve (save_C_processed);
			prime = sieve ();
			goto restart3;
		}
		
/* We've finished stage 2, but haven't done the GCD yet */

		goto restart4;
	}

/* Loop processing the requested number of ECM curves */

restart0:
	last_output = fft_count = modinv_count = 0;

/* Allocate memory */

	x = gwalloc ();
	z = gwalloc ();
	gg = NULL;

/* Choose curve with order divisible by 16 and choose a point (x/z) on */
/* said curve. */

	do {
		unsigned long hi, lo;
		sigma = (rand () & 0x1F) * 65536.0 * 65536.0 * 65536.0;
		sigma += (rand () & 0xFFFF) * 65536.0 * 65536.0;
		if (CPU_FLAGS & CPU_RDTSC) rdtsc (&hi, &lo);
		sigma += lo ^ hi ^ ((unsigned long) rand () << 16);
	} while (sigma <= 5.0);
	if (specific_sigma > 5.0) sigma = specific_sigma;
	curve_start_msg (curve, sigma, B, C);
	choose12 (k, b, n, c, x, z, sigma);
	sieve_start = 2;

/* The stage 1 restart point */

restart1:
	stage = 1;
	start_timer (0);
	start_sieve (sieve_start);
	for ( ; ; ) {
		prime = sieve ();
		if (prime > B) break;

/* Apply as many powers of prime as long as prime^n <= B */
/* MEMUSED: 3 gwnums (x, z, AD4) + 10 for ell_mul */

		ell_mul (x, z, prime);
		if (prime <= SQRT_B) {
			unsigned long mult, max;
			mult = prime;
			max = B / prime;
			for ( ; ; ) {
				ell_mul (x, z, prime);
				mult *= prime;
				if (mult > max) break;
			}
		}

/* Print a message every so often */

		if (ITER_OUTPUT != 999999999 &&
		    fft_count >= last_output + 2 * ITER_OUTPUT) {
			sprintf (buf, "At prime %lu.  Time thusfar ", prime);
			OutputTimeStamp ();
			OutputStr (buf);
			end_timer (0);
			print_timer (0, TIMER_NL);
			start_timer (0);
			last_output = fft_count;
		}

/* Check for errors */

		if (gw_test_for_error ()) goto error;

/* Write a save file when the user interrupts the calculation and */
/* every DISK_WRITE_TIME minutes. */

		escaped = stopCheck ();
		time (&current_time);
		if (escaped || (QA_SAVE_FILES & 0x1) ||
		    current_time - start_time > write_time) {
			ecm_save (filename, ECM_STAGE1, curve, sigma,
				  B, prime, 0, x, z);
			if (escaped) {
				retval = FALSE;
				goto exit;
			}
			if (QA_SAVE_FILES & 0x1) {
				QA_SAVE_FILES &= ~0x1;
				retval = QA_SAVE_TEST;
				goto exit;
			}
			start_time = current_time;
		}
	}

/* Stage 1 complete */

	end_timer (0);
	sprintf (buf, "Stage 1 complete. %.0f transforms, %lu modular inverses. Time: ",
		 fft_count, modinv_count);
	OutputStr (buf);
	print_timer (0, TIMER_NL | TIMER_CLR);
	last_output = fft_count = modinv_count = 0;

/* Print out round off error */

	if (ERRCHK) {
		sprintf (buf, "Round off: %.10g\n", MAXERR);
		OutputStr (buf);
		MAXERR = 0.0;
	}

/* If we aren't doing a stage 2, then check to see if we found a factor. */
/* If we are doing a stage 2, then the stage 2 init will do this GCD for us. */

	if (C <= B) {
		start_timer (0);
		if ((QA_SAVE_FILES & 0x2) || !gcd (z)) {
			ecm_save (filename, ECM_STAGE1, curve, sigma,
				  B, B, 0, x, z);
			if (QA_SAVE_FILES & 0x2) {
				QA_SAVE_FILES &= ~0x2;
				retval = QA_SAVE_TEST;
				goto exit;
			}
			retval = FALSE;
			goto exit;
		}
		end_timer (0);
		OutputStr ("Stage 1 GCD complete. Time: ");
		print_timer (0, TIMER_NL | TIMER_CLR);
		if (FAC != NULL) goto bingo;

/* Alexander Kruppa wrote this code to normalize and output the x value */
/* along with N and sigma so that it can be used in Paul Zimmermann's */
/* superior GMP-ECM implementation of stage 2. */

		if (IniGetInt (INI_FILE, "GmpEcmHook", 0)) {
			char	*msg, *buf;
			int	msglen, i;
			giant	gx;
			char	*hex = "0123456789ABCDEF";

			normalize (x, z);

			gx = popg ((PARG >> 5) + 10);
			gwtogiant (x, gx);
			modg (N, gx);

			msglen = N->sign * 8 + 5;
			buf = (char *) malloc (msglen + msglen + 80);
			if (buf != NULL) {
				strcpy (buf, "N=");
				msg = buf + strlen (buf);
				for (i = 0; i < N->sign * 8; i++) {
				  msg[i+2] = hex[ ( N->n[N->sign - 1 - i/8] 
				                    >> ((7-i%8)*4)
				                  ) & 0xF ];
				}
				msg[i+2] = 0;
				for (i = 0; i < N->sign * 8 - 1; i++)
				  if (msg[i+2] != '0') break;
				msg[i] = '0'; msg[i+1] = 'x';
				strcpy (msg, msg+i);
				strcat (buf, "; QX=");
				msg = buf + strlen (buf);
				for (i = 0; i < gx->sign * 8; i++) {
				  msg[i+2] = hex[ ( gx->n[gx->sign - 1 - i/8] 
				                    >> ((7-i%8)*4)  
				                  ) & 0xF ];
				}
				msg[2+i] = 0;
				for (i = 0; i < gx->sign * 8 - 1; i++)
				  if (msg[i+2] != '0') break;
				msg[i] = '0'; msg[i+1] = 'x';
				strcpy (msg, msg+i);
				strcat (buf, "; SIGMA=");
				msg = buf + strlen (buf);
				sprintf (msg, "%.0f\n", sigma);
				writeResults (buf);
				free (buf);
			}
			pushg (1);
		}

/* Now do the next ECM curve */

		goto more_curves;
	}

/*
   Stage 2:  We support two types of stage 2's here.  One uses
   less memory and uses fewer extended GCDs, but is slower in accumulating
   each found prime.  Thanks to Richard Crandall and Paul Zimmermann
   for letting me liberally use their code and ideas here.
   x, z: coordinates of Q at the beginning of stage 2
*/

/* If we find a factor in stage 1 of a QA run then */
/* print the factor, divide it from N, and continue stage 2. */
/* On a QA continuation z is one, so do the GCD on x */

restart3:
	if (QA_IN_PROGRESS) {
		gwnum	p;
		p = (*(double *)z == 1.0) ? x : z;
		if (!gcd (p)) {
			retval = FALSE;
			goto exit;
		}
		if (FAC != NULL) {
			printFactor (k, b, n, c);
			divg (FAC, N);
			free (FAC);
			FAC = NULL;
		}
	}

/* Initialize variables for second stage */
/* Our goal is to fill up the nQx array with Q^1, Q^3, Q^5, ... */
/* normalized with only one modular inverse call. */

	start_timer (0);

/* Allocate memory for computing nQx values */
/* MEMUSED: 9 gwnums (x, z, AD4, 6 for nQx) */

	Q2x = gwalloc (); Q2z = gwalloc ();
	Qiminus2x = gwalloc (); Qiminus2z = gwalloc ();
	Qdiffx = gwalloc (); Qdiffz = gwalloc ();

/* Init values used in computing nQx.  We need Q^2, Q^1, and diff of Q^1. */
/* MEMUSED: 9 gwnums (x, z, AD4, 6 for computing nQx) + 2 temporaries */

	ell_dbl (x, z, Q2x, Q2z);
	ell_begin_fft (Q2x, Q2z, Q2x, Q2z);
	gwfft (x, x); gwfft (z, z);
	gwcopy (x, Qdiffx); gwcopy (z, Qdiffz);
	gwcopy (x, Qiminus2x); gwcopy (z, Qiminus2z);

/* Init the first nQx value with Q^1 */
/* MEMUSED: 9 gwnums (AD4, 6 for computing nQx, nQx[0], modinv_value) */

	nQx[0] = x;
	add_to_normalize_pool (nQx[0], z, 1);
	gwfree (z);

/* Compute the rest of the nQx values (Q^i for i >= 3) */
/* MEMUSED: 8 + nQx gwnums (AD4, 6 for computing nQx, nQx vals, modinv_val) */
/* MEMPEAK: 8 + nQx-1 + 2 for edd_add temporaries */

	for (i = 3; i < D; i = i + 2) {
		ell_add_special (Qiminus2x, Qiminus2z, Q2x, Q2z,
				 Qdiffx, Qdiffz, Qdiffx, Qdiffz);

		if (gw_test_for_error ()) goto error;

		if (stopCheck () || (QA_SAVE_FILES & 0x4)) {
			if (gg == NULL) {
				ecm_save (filename, ECM_STAGE1, curve, sigma,
					  B, B, 0, Qdiffx, Qdiffz);
			}
			if (QA_SAVE_FILES & 0x4) {
				QA_SAVE_FILES &= ~0x4;
				retval = QA_SAVE_TEST;
				goto exit;
			}
			retval = FALSE;
			goto exit;
		}

		gwfft (Qdiffx, Qdiffx); gwfft (Qdiffz, Qdiffz);
		if (relatively_prime (i, D)) {
			j = (i - 1) >> 1;
			nQx[j] = gwalloc ();
			gwcopy (Qdiffx, nQx[j]);
			add_to_normalize_pool (nQx[j], Qdiffz, 1);
		}
		gwswap (Qdiffx, Qiminus2x); gwswap (Qdiffz, Qiminus2z);
	}

/* Now compute Q^2D.  This will be used in computing Q^m values. */
/* Qiminus2 is Q^(D-1) and Qdiff if Q^(D-3).  Add Q^(D-1) and Q^2 to */
/* get Q^(D+1).  Then add Q^(D-1) and Q^(D+1) to get Q^2D.  Store Q^2D */
/* in Q2x and Q2z.  Normalize them so we can free Q2z later on. */
/* MEMUSED: 8 + nQx gwnums (AD4, 6 for computing nQx, nQx vals, modinv_val) */
/* MEMPEAK: 8 + nQx + 2 for edd_add temporaries */

	ell_add_special (Qiminus2x, Qiminus2z, Q2x, Q2z,
			 Qdiffx, Qdiffz, Qdiffx, Qdiffz);
	gwfftaddsub (Q2x, Q2z);			/* Recompute fft of Q2x,Q2z */
	ell_begin_fft (Qdiffx, Qdiffz, Qdiffx, Qdiffz);
	ell_add_special (Qiminus2x, Qiminus2z, Qdiffx, Qdiffz,
			 Q2x, Q2z, Q2x, Q2z);
	gwfft (Q2x, Q2x); gwfft (Q2z, Q2z);
	add_to_normalize_pool (Q2x, Q2z, 1);

/* Free most of the memory used in computing nQx values */
/* Keep two values we could free in case the upcoming normalize is */
/* aborted and we need to write a save file. */
/* MEMUSED: 5 + nQx gwnums (AD4, Q2x, Qiminus2x&z, nQx values, modinv_value) */

	gwfree (Q2z);
	gwfree (Qdiffx); gwfree (Qdiffz);

/* Normalize all the nQx values */
/* MEMUSED: 4 + nQx gwnums (AD4, Q2x, Qiminus2x&z, nQx values) */
/* MEMPEAK: 4 + nQx + 3 for UV, normalize, and pooled_modinv temporaries */

	if (! normalize_pool () || (QA_SAVE_FILES & 0x8)) {
		if (gg == NULL) {
			dbltogw (1.0, Q2x);
			gwfft (Q2x, Q2x);
			gwfftfftmul (Q2x, Qiminus2x, Qiminus2x);
			gwfftfftmul (Q2x, Qiminus2z, Qiminus2z);
			ecm_save (filename, ECM_STAGE1, curve, sigma,
				  B, B, 0, Qiminus2x, Qiminus2z);
		}
		if (QA_SAVE_FILES & 0x8) {
			QA_SAVE_FILES &= ~0x8;
			retval = QA_SAVE_TEST;
			goto exit;
		}
		retval = FALSE;
		goto exit;
	}

/* If we found a factor, we're done */

	if (FAC != NULL) goto bingo;

/* Free rest of the memory used in computing nQx values */
/* MEMUSED: 2 + nQx gwnums (AD4, Q2x, nQx values) */

	gwfree (Qiminus2x); gwfree (Qiminus2z);

/* Init code that computes Q^m */
/* MEMUSED: 6 + nQx gwnums (6 for computing mQx, nQx values) */
/* MEMPEAK: 6 + nQx + 6 for bin_ell_mul temporaries */

	m = (prime / D + 1) * D;
	mQ_init (nQx[0], m, Q2x);

/* Precompute the transforms of nQx */

	for (i = 0; i < D/2; i++) if (nQx[i] != NULL) gwfft (nQx[i], nQx[i]);

/* Now init the accumulator unless this value was read */
/* from a continuation file */
/* MEMUSED: 7 + nQx gwnums (6 for computing mQx, gg, nQx values) */

	if (gg == NULL) {
		gg = gwalloc ();
		dbltogw (1.0, gg);
	}

/* Initialization of stage 2 complete */

	end_timer (0);
	sprintf (buf, "Stage 2 init complete. %.0f transforms, %lu modular inverses. Time: ",
		 fft_count, modinv_count);
	OutputStr (buf);
	print_timer (0, TIMER_NL | TIMER_CLR);
	last_output = fft_count = modinv_count = 0;

/* Now do stage 2 */
/* Accumulate (mQx - nQx)(mQz + nQz) - mQx mQz + nQx nQz.		*/
/* Since nQz = 1, we have (the 4 FFT per prime continuation)		*/
/*		== (mQx - nQx)(mQz + 1) - mQx mQz + nQx			*/
/*		== mQx mQz - nQx mQz + mQx - nQx - mQx mQz + nQx	*/
/*		== mQx - nQx mQz					*/
/* If mQz also = 1 (the 2 FFT per prime continuation) then we accumulate*/
/*		== mQx - nQx						*/

	start_timer (0);
	stage = 2;
	for ( ; C > m-D; m += D+D) {
		gwnum	mQx, mQz;

/* Compute next Q^m value */
/* MEMUSED: 7 + nQx + E gwnums (6 for computing mQx, gg, nQx and E values) */
/* MEMPEAK: 7 + nQx + E + 2 for ell_add temporaries */

		if (!mQ_next (&mQx, &mQz) || (QA_SAVE_FILES & 0x10)) {
			dbltogw (1.0, Q2x);
			gwfft (Q2x, Q2x);
			gwfftfftmul (Q2x, nQx[0], nQx[0]);
			ecm_save (filename, ECM_STAGE2, curve, sigma,
				  B, B, prime, nQx[0], gg);
			if (QA_SAVE_FILES & 0x10) {
				QA_SAVE_FILES &= ~0x10;
				retval = QA_SAVE_TEST;
				goto exit;
			}
			retval = FALSE;
			goto exit;
		}
		if (FAC != NULL) goto bingo;
		memset (pairings, 0, (D + 15) >> 4);
		t1 = gwalloc ();

/* 2 FFT per prime continuation - deals with all normalized values */

		if (TWO_FFT_STAGE2) {
		    for ( ; ; prime = sieve ()) {
			if (prime < m) {	/* Do the m-D to m range */
				i = (m - prime) >> 1;
				bitset (pairings, i);
			} else if (prime < m+D) { /* Do the D to m+D range */
				i = (prime - m) >> 1;
				if (bittst (pairings, i)) continue;
			} else
				break;
			gwfftsub3 (mQx, nQx[i], t1);
			gwfftmul (t1, gg);
		    }
		}

/* 4 FFT per prime continuation - deals with only nQx values normalized */

		else {
		    for ( ; ; prime = sieve ()) {
			if (prime < m) {	/* Do the m-D to m range */
				i = (m - prime) >> 1;
				bitset (pairings, i);
			} else if (prime < m+D) { /* Do the D to m+D range */
				i = (prime - m) >> 1;
				if (bittst (pairings, i)) continue;
			} else
				break;
			gwfftfftmul (nQx[i], mQz, t1);
			gwfft (t1, t1);
			gwfftsub3 (mQx, t1, t1);
			gwfftmul (t1, gg);
		    }
		}
		gwfree (t1);

/* Print a message every so often */

		if (ITER_OUTPUT != 999999999 &&
		    fft_count >= last_output + 2 * ITER_OUTPUT) {
			sprintf (buf, "At prime %lu.  Time thusfar ", prime);
			OutputTimeStamp ();
			OutputStr (buf);
			end_timer (0);
			print_timer (0, TIMER_NL);
			start_timer (0);
			last_output = fft_count;
		}

/* Check for errors */

		if (gw_test_for_error ()) goto error;

/* Write a save file when the user interrupts the calculation and */
/* every DISK_WRITE_TIME minutes. */

		escaped = stopCheck ();
		time (&current_time);
		if (escaped || (QA_SAVE_FILES & 0x20) ||
		    current_time - start_time > write_time) {
			t1 = gwalloc ();
			dbltogw (1.0, t1);
			gwfft (t1, t1);
			gwfftfftmul (t1, nQx[0], t1);
			ecm_save (filename, ECM_STAGE2, curve, sigma,
				  B, B, prime, t1, gg);
			gwfree (t1);
			if (escaped) {
				retval = FALSE;
				goto exit;
			}
			if (QA_SAVE_FILES & 0x20) {
				QA_SAVE_FILES &= ~0x20;
				retval = QA_SAVE_TEST;
				goto exit;
			}
			start_time = current_time;
		}
	}
	mQ_term ();

/* Stage 2 is complete */

	end_timer (0);
	sprintf (buf, "Stage 2 complete. %.0f transforms, %lu modular inverses. Time: ",
		 fft_count, modinv_count);
	OutputStr (buf);
	print_timer (0, TIMER_NL | TIMER_CLR);
	last_output = fft_count = modinv_count = 0;

/* Print out round off error */

	if (ERRCHK) {
		sprintf (buf, "Round off: %.10g\n", MAXERR);
		OutputStr (buf);
		MAXERR = 0.0;
	}

/* See if we got lucky! */

restart4:
	start_timer (0);
	if ((QA_SAVE_FILES & 0x40) || !gcd (gg)) {
		ecm_save (filename, ECM_STAGE2, curve, sigma,
			  B, B, C, gg, gg);
		if (QA_SAVE_FILES & 0x40) {
			QA_SAVE_FILES &= ~0x40;
			retval = QA_SAVE_TEST;
			goto exit;
		}
		retval = FALSE;
		goto exit;
	}
	end_timer (0);
	OutputStr ("Stage 2 GCD complete. Time: ");
	print_timer (0, TIMER_NL | TIMER_CLR);
	if (FAC != NULL) goto bingo;

/* Do not loop if we are testing a specific curve */

more_curves:
	gwfreeall ();
	if (specific_sigma < 5.0 && ++curve <= curves_to_do)
		goto restart0;

/* Output line to results file indicating the number of curves run */

	sprintf (buf, "%s completed %ld ECM curves, B1=%lu, B2=%lu\n",
		 gwmodulo_as_string (), curves_to_do, B, C);
	writeResults (buf);
	/*spoolMessage (PRIMENET_RESULT_MESSAGE, buf);*/
	retval = TRUE;

/* Delete the exponent from the work-to-do-file */

	updateWorkToDo (k, b, n, c, WORK_ECM, orig_B);
	_unlink (filename);

/* Free memory and return */

exit:	ecm_cleanup ();
	return (retval);

/* Print a message if we found a factor! */

bingo:	sprintf (buf, "ECM found a factor in curve #%ld, stage #%d\n",
		 curve, stage);
	writeResults (buf);
	/*spoolMessage (PRIMENET_RESULT_MESSAGE, buf);*/
	sprintf (buf, "Sigma=%.0f, B1=%lu, B2=%lu.\n", sigma, B, C);
	writeResults (buf);
	/*spoolMessage (PRIMENET_RESULT_MESSAGE, buf);*/
	if (printFactor (k, b, n, c)) {
		_unlink (filename);
		if (!QA_IN_PROGRESS) {
			free (FAC);
			FAC = NULL;
		}
		retval = TRUE;
		goto exit;
	}
	clear_timer (0);

	if (! IniGetInt (INI_FILE, "ContinueECM", n < 5825)) {
		updateWorkToDo (k, b, n, c, WORK_FACTOR, 0);
		_unlink (filename);
		free (FAC);
		FAC = NULL;
		retval = TRUE;
		goto exit;
	}

	divg (FAC, N);
	free (FAC);
	FAC = NULL;
	goto more_curves;

/* Output a message saying we are restarting */

error:	OutputBoth ("SUMOUT error occurred.\n");

/* Sleep five minutes before restarting */

	ecm_cleanup ();
	if (! SleepFive ()) return (FALSE);

/* Restart from last save file */

	goto restart;
}

/* Read a file of ECM tests to run as part of a QA process */
/* The format of this file is: */
/*	k, n, c, sigma, B1, B2_start, B2_end, factor */
/* Use Advanced/Time 9991 to run the QA suite */

int ecm_QA ()
{
	FILE	*fd;
	int	savefiles;

/* Set the title */

	title ("QA");

/* Open QA file */

	fd = fopen ("qa_ecm", "r");
	if (fd == NULL) {
		OutputStr ("File named 'qa_ecm' could not be opened.\n");
		return (TRUE);
	}

/* Loop until the entire file is processed */

	for ( ; ; ) {
		double	k;
		unsigned long b, n, B1, B2_start, B2_end;
		signed long c;
		char	fac_str[80];
		double	sigma;
		giant	f;
		int	retval, success;

/* Read a line from the file */

		n = 0;
		fscanf (fd, "%f,%lu,%lu,%ld,%lf,%lu,%lu,%lu,%s\n",
			&k, &b, &n, &c, &sigma, &B1, &B2_start, &B2_end, &fac_str);
		if (n == 0) break;

/* If b is 1, set QA_TYPE */

		if (b == 1) {
			QA_TYPE = c;
			savefiles = B1;
			continue;
		}

/* Convert the factor we expect to find into a "giant" type */

		f = newgiant (strlen (fac_str));
		ctog (fac_str, f);

/*test various num_tmps
test 4 (or more?) stage 2 code paths
print out each test case (all relevant data)*/

/* Set some global variables indicating QA is in progress. */

		QA_IN_PROGRESS = TRUE;
		QA_SAVE_FILES = savefiles;

/* Do the ECM */

		if (B2_start < B1) B2_start = B1;
		do {
			retval = ecm (k, b, n, c, B1, B2_start, B2_end, 1, sigma);
			if (!retval) {
				fclose (fd);
				QA_IN_PROGRESS = FALSE;
				QA_TYPE = 0;
				QA_SAVE_FILES = 0;
				return (FALSE);
			}
		} while (retval == QA_SAVE_TEST);

/* See if we found the factor */

		if (FAC == NULL) success = FALSE;
		else {
			modg (f, FAC);
			success = isZero (FAC);
			free (FAC);
			FAC = NULL;
		}
		if (!success) {
			char	buf[200];
			sprintf (buf, "ERROR: Factor not found. Expected %s\n", fac_str);
			OutputBoth (buf);
		}
		free (f);
	}

/* Cleanup */

	QA_IN_PROGRESS = FALSE;
	QA_TYPE = 0;
	QA_SAVE_FILES = 0;
	fclose (fd);
	return (TRUE);
}

/**************************************************************
 *
 *	P-1 Functions
 *
 **************************************************************/

/* Perform setup functions.  This includes decding how big an FFT to */
/* use, allocating memory, calling the FFT setup code, etc. */

int pm1_setup1 (
	double	k,
	unsigned long b,
	unsigned long n,
	signed long c)
{
	unsigned long fftlen;
	int	res;

/* Setup the assembly code */

	fftlen = (k == 1.0 && b == 2 && c == -1) ?
		advanced_map_exponent_to_fftlen (n) : 0;
	res = gwsetup (k, b, n, c, fftlen);
	if (res) return (res);
	modinv_count = 0;
	N = NULL;
	nQx = NULL;
	eQx = NULL;

/* A kludge so that the error checking code is not as strict. */

	MAXDIFF *= IniGetInt (INI_FILE, "MaxDiffMultiplier", 1);
	return (0);
}

/* Perform cleanup functions. */

void pm1_cleanup (void)
{

/* Free memory */

	free (N);
	free (nQx);
	free (eQx);
	gwdone ();
}

/* Raises number to the given power */

void pm1_mul (
	gwnum	xx,
	unsigned long n)
{
	gwnum	orig_xx_fft;
	unsigned long c;

/* Find most significant bit and then ignore it */

	c = (unsigned long)(1<<31);
	while ((c&n) == 0) c >>= 1;
	c >>= 1;

/* Handle the second most significant bit */

	orig_xx_fft = gwalloc ();
	gwfft (xx, orig_xx_fft);
	gwfftfftmul (orig_xx_fft, orig_xx_fft, xx);
	if (c&n) gwfftmul (orig_xx_fft, xx);
	c >>= 1;

/* Do the rest of the bits */

	while (c) {
		gwsquare (xx);
		if (c&n) gwfftmul (orig_xx_fft, xx);
		c >>= 1;
	}
	gwfree (orig_xx_fft);
}

/* Code to init "finite differences" for computing successive */
/* values of x^(start+i*incr)^E */

void fd_init (
	unsigned long start,
	unsigned long incr,
	gwnum	x)		/* Caller must pass in the FFT of x */
{
	unsigned long i, j;
	giant	p;

/* Treat each eQx[i] as a binary value and compute (start+i*incr)^e */

	for (i = 0; i <= E; i++) {
		unsigned long val;
		p = newgiant (E * 2);
		val = start + i * incr;
		ultog (val, p);
		for (j = 2; j <= E; j++) ulmulg (val, p);
		eQx[i] = (gwnum) p;
	}		

/* Now do the finite differences */

	for (i = 1; i <= E; i++) {
		for (j = E; j >= i; j--) {
			subg ((giant) eQx[j-1], (giant) eQx[j]);
		}
	}

/* Now compute each x^difference */

	for (i = 0; i <= E; i++) {
		p = (giant) eQx[i];
		eQx[i] = gwalloc ();

/* Test for easy cases */

		ASSERTG (!isZero (p));
		if (isone (p)) {
			gwcopy (x, eQx[i]);
		}

/* Find most significant bit and then ignore it */

		else {
			int	len;

			len = bitlen (p);
			len--;

/* Perform the first squaring using the already FFTed value of x */
/* Then process the second and remaining bits of p */

			gwfftfftmul (x, x, eQx[i]);
			for ( ; ; ) {
				if (bitval (p, len-1)) gwfftmul (x, eQx[i]);
				len--;
				if (len == 0) break;
				gwsquare (eQx[i]);
			}

/* FFT the final result */

			gwfft (eQx[i], eQx[i]);
		}
		free (p);
	}
}

/* Code to compute next x^(start+i*incr)^E value */
/* Value is returned in eQx[0] - already FFTed */

void fd_next (void)
{
	unsigned long i;

	for (i = 0; i < E; i++) {
		gwfftfftmul (eQx[i], eQx[i+1], eQx[i]);
		gwfft (eQx[i], eQx[i]);
	}	
}

/* Terminate finite differences code */

void fd_term (void)
{
	unsigned long i;

/* Free each eQx[i] */

	for (i = 0; i <= E; i++) gwfree (eQx[i]);
}

/* Routines to create and read save files for a P-1 factoring job */

#define PM1_STAGE0	3	/* In stage 1, squaring small primes */
#define PM1_STAGE1	0	/* In stage 1, processing larger primes */
#define PM1_STAGE2	1	/* In stage 2 */
#define PM1_DONE	2	/* P-1 job complete */

struct pm1_state {
	unsigned long stage;	/* One of the 4 states listed above */
	unsigned long B_done;	/* We have completed calculating 3^e */
				/* to this bound #1 */
	unsigned long B;	/* We are trying to increase bound #1 */
				/* to this value */
	unsigned long C_done;	/* We have completed to this bound #2 */
	unsigned long C_start;	/* We are trying to increase bound #2 */
				/* from this starting point. */
	unsigned long C;	/* We are trying to increase bound #2 */
				/* to this end point */
	unsigned long numrels;	/* In a multi-pass stage 2, the number */
				/* of values relatively prime to D */
	unsigned long rels_done;/* In a multi-pass stage 2, the number */
				/* of relative primes already processed */
	unsigned long rels_this_pass;
				/* In a multi-pass stage 2, the number */
				/* of relative primes we are processing */
				/* this pass. */
	char	*bitarray;	/* Bit array for primes between C_start */
				/* and C */
	unsigned long bitarray_len;
				/* Number of bytes in the bit array */
	unsigned long pairs_set;/* Number of pairs originally set in */
				/* bitarray */
	unsigned long pairs_done;/* Number of pairs completed */
};

void pm1_save (
	char	*filename,
	struct pm1_state *state,
	unsigned long processed,
	gwnum	x,
	gwnum	gg)
{
	char	newfilename[16];
	int	fd;
	unsigned long magicnum, version;
	long	sum = 0, i;

/* If we are allowed to create multiple intermediate files, then */
/* write to a file called xNNNNNNN. */

	strcpy (newfilename, filename);
	if (TWO_BACKUP_FILES) newfilename[0] = 'x';

/* Create the intermediate file */

	fd = _open (newfilename, _O_BINARY|_O_WRONLY|_O_TRUNC|_O_CREAT, CREATE_FILE_ACCESS);
	if (fd < 0) return;

/* Write the file header */

	magicnum = 0x1a2b3c4d;
	if (_write (fd, &magicnum, sizeof (long)) != sizeof (long))
		goto writeerr;
	version = 4;
	if (_write (fd, &version, sizeof (long)) != sizeof (long))
		goto writeerr;

/* Write the file data */

	i = state->stage;
	if (! write_long (fd, i, &sum)) goto writeerr;
	if (! write_long (fd, state->B_done, &sum)) goto writeerr;
	if (! write_long (fd, state->B, &sum)) goto writeerr;
	if (! write_long (fd, state->C_done, &sum)) goto writeerr;
	if (! write_long (fd, state->C_start, &sum)) goto writeerr;
	if (! write_long (fd, state->C, &sum)) goto writeerr;
	if (! write_long (fd, processed, &sum)) goto writeerr;
	if (! write_long (fd, D, &sum)) goto writeerr;
	if (! write_long (fd, E, &sum)) goto writeerr;
	if (! write_long (fd, state->rels_done, &sum)) goto writeerr;
	if (! write_long (fd, state->bitarray_len, &sum)) goto writeerr;
	if (! write_array (fd, state->bitarray, state->bitarray_len, &sum))
		goto writeerr;
	if (! write_long (fd, state->pairs_set, &sum)) goto writeerr;
	if (! write_long (fd, state->pairs_done, &sum)) goto writeerr;

/* Write the data values */

	if (! write_gwnum (fd, x, &sum)) goto writeerr;
	if (gg != NULL && ! write_gwnum (fd, gg, &sum)) goto writeerr;

/* Write the checksum */

	if (_write (fd, &sum, sizeof (long)) != sizeof (long)) goto writeerr;
	_commit (fd);
	_close (fd);

/* Now rename the intermediate files */

	if (TWO_BACKUP_FILES) {
		_unlink (filename);
		rename (newfilename, filename);
	}
	return;

/* An error occured.  Close and delete the current file. */

writeerr:
	_close (fd);
	_unlink (newfilename);
}

/* Read a save file */

int pm1_restore (
	char	*filename,
	struct pm1_state *state,
	unsigned long *processed,
	gwnum	*x,
	gwnum	*gg)
{
	int	fd;
	unsigned long magicnum, version;
	long	sum = 0, i;

/* Open the intermediate file */

	fd = _open (filename, _O_BINARY | _O_RDONLY);
	if (fd < 0) goto error;

/* Read the file header */

	if (_read (fd, &magicnum, sizeof (long)) != sizeof (long))
		goto readerr;
	if (magicnum != 0x1a2b3c4d) goto readerr;

	if (_read (fd, &version, sizeof (long)) != sizeof (long))
		goto readerr;
	if (version < 3 || version > 4) goto readerr;

/* Read the file data */

	if (! read_long (fd, &state->stage, &sum)) goto readerr;

/* We can no longer process version 3 stage 2 save files.  We'll have to */
/* start stage 2 from scratch. */

	if (version == 3) {
		unsigned long numvals, relprime;
		if (! read_long (fd, &state->B_done, &sum)) goto readerr;
		if (! read_long (fd, &state->B, &sum)) goto readerr;
		if (! read_long (fd, &state->C_done, &sum)) goto readerr;
		if (! read_long (fd, &state->C_start, &sum)) goto readerr;
		if (! read_long (fd, &state->C, &sum)) goto readerr;
		if (! read_long (fd, processed, &sum)) goto readerr;
		if (! read_long (fd, &relprime, &sum)) goto readerr;
		if (! read_long (fd, &numvals, &sum)) goto readerr;
		state->bitarray_len = 0;
		state->bitarray = NULL;
		if (state->stage == PM1_STAGE2) {
			if (numvals == 0 || (relprime == 1 && numvals >= 12))
				state->C_done = *processed;
			state->stage = PM1_STAGE2;
		}
	}

	if (version == 4) {
		if (! read_long (fd, &state->B_done, &sum)) goto readerr;
		if (! read_long (fd, &state->B, &sum)) goto readerr;
		if (! read_long (fd, &state->C_done, &sum)) goto readerr;
		if (! read_long (fd, &state->C_start, &sum)) goto readerr;
		if (! read_long (fd, &state->C, &sum)) goto readerr;
		if (! read_long (fd, processed, &sum)) goto readerr;
		if (! read_long (fd, &D, &sum)) goto readerr;
		if (! read_long (fd, &E, &sum)) goto readerr;
		if (! read_long (fd, &state->rels_done, &sum)) goto readerr;
		if (! read_long (fd, &state->bitarray_len, &sum)) goto readerr;
		if (state->bitarray_len) {
			state->bitarray = (char *) malloc (state->bitarray_len);
			if (state->bitarray == NULL) goto readerr;
			if (! read_array (fd, state->bitarray,
					  state->bitarray_len, &sum))
				goto readerr;
		}
		if (! read_long (fd, &state->pairs_set, &sum)) goto readerr;
		if (! read_long (fd, &state->pairs_done, &sum)) goto readerr;
	}

/* Read the values */

	*x = gwalloc ();
	if (! read_gwnum (fd, *x, &sum)) goto readerr;

	*gg = NULL;
	if (state->stage == PM1_STAGE2) {
		*gg = gwalloc ();
		if (! read_gwnum (fd, *gg, &sum)) goto readerr;
	}

/* Read and compare the checksum */

	if (_read (fd, &i, sizeof (long)) != sizeof (long)) goto readerr;
	if (i != sum) goto readerr;
	_close (fd);
	return (TRUE);

/* An error occured.  Delete the current intermediate file. */
/* Set stage to -1 to indicate an error. */

readerr:
	OutputStr ("Error reading P-1 save file.\n");
	_close (fd);
error:
	_unlink (filename);
	return (FALSE);
}

/* Compute how many values we can allocate.  This function can calculate */
/* the value using either the maximum available memory or the currently */
/* available memory. */

unsigned long choose_pminus1_numvals (
	double	k,
	unsigned long b,
	unsigned long n,
	signed long c,
	int	use_max_mem)		/* True if calculation should use */
					/* maximum available memory */
{
	unsigned long memory;		/* Available memory in MB */
	double	temp, size;

/* Override numvals when QAing */

	if (QA_TYPE) return (QA_TYPE);

/* Default is to use up to 24MB of memory */

	memory = use_max_mem ? max_mem () : avail_mem ();
	if (memory == 0) memory = 24;

/* Compute the number of gwnum temporaries available */

	temp = (double) memory * 1000000.0 - 
		(double) gwmap_to_memused (k, b, n, c);
	size = (double) gwnum_size (FFTLEN);
	if (temp <= size) return (1);
	return ((unsigned long) (temp / size));
}

/* Calculate the number of values relatively prime to D.  D is a multiple */
/* of 30, 210, or 2310. */

unsigned long calc_numrels (
	unsigned long d)
{
	if (d >= 2310) return (d / 2310 * 480);
	if (d >= 210) return (d / 210 * 48);
	return (d / 30 * 8);
}

/* Compute the cost of a particular P-1 stage 2 plan. */

double cost_pminus1_plan (
	unsigned long B,		/* Stage 2 start point */
	unsigned long C,		/* Stage 2 end point */
	unsigned long d,
	unsigned long e,		/* Suyama power */
	unsigned long numvals,		/* Temp gwnums available for stage 2 */
	int	using_t3)		/* True if t3 avoids a gwfftadd3 */
{
	unsigned long numprimes, numpairings, numrels, passes, sets;
	double	cost;

/* Estimate the number of primes */

	numprimes = (unsigned long) (C / (log((double) C) - 1.0) - B / (log((double) B) - 1.0));

/* Calculate the number of values relatively prime to D */

	numrels = calc_numrels (d);

/* Estimate the number of prime pairs */

	if (e >= 2)
		numpairings = (unsigned long)
			(numprimes / 2.0 *
			 numprimes / ((C-B) * (double) numrels / d));

/* Compute how many passes this will take to process all the numbers */
/* relatively prime to D. */

	passes = (unsigned long)
		ceil ((double) numrels / (numvals - (e+1) - using_t3));

/* Compute the nQx setup costs.  To calculate nQx values, one fd_init is */
/* required on each pass with an average start point of D/2.  A single */
/* fd_init does E+1 powerings of roughly E*log2(startpoint) bits each */
/* at a cost of about 1.5 multiplies per bit.  In other words, */
/* (E+1) * E*log2(startpoint) * 1.5. */

	cost = passes * (e+1) * e*log((double)(d/2))/log((double)2.0) * 1.5 +

/* Then there are the D/2 calls to fd_next at E multiplies each. */

		d/2 * e +

/* Compute the eQx setup costs.  To calculate eQx values, one fd_init is */
/* required on each pass with a start point of B. */

		passes * (e+1) * e*log((double)B)/log((double)2.0) * 1.5;

/* If E=1 add the cost of (C-B)/D fd_next calls.  If E>=2, add the cost */
/* of (C-B)/(D+D) calls to fd_next (E multiplies). */

	if (e == 1)
		sets = (unsigned long) ceil ((double) (C-B) / d);
	else
		sets = (unsigned long) ceil ((double) (C-B) / (d+d));
	cost += passes * sets * e;

/* Finally, each prime pairing costs one multiply.  If E = 1, then there */
/* is no prime pairing.  If not using_t3 then each multiply costs more as */
/* there is an extra gwfftadd3 call. */

	if (e >= 2 && using_t3)
		cost += numprimes - numpairings;
	else if (e >= 2)
		cost += (numprimes - numpairings) * 1.1;
	else if (using_t3)
		cost += numprimes;
	else
		cost += numprimes * 1.1;

/* Return the resulting cost */

	return (cost);
}

/* Choose the best values for D and E.  One that reduces the number of */
/* multiplications, yet doesn't use too much memory. */

void choose_pminus1_plan (
	double	k,
	unsigned long b,
	unsigned long n,
	signed long c,
	struct pm1_state *state)	/* P-1 state structure */
{
	unsigned long B, C, numvals, d, e, i;
	double	cost, best_cost;

/* Handle case where there is no stage 2 */

	B = state->C_start;		/* Stage 2 starting point */
	C = state->C;			/* Stage 2 ending point */
	if (C <= B) {
		D = 0;
		E = 0;
		return;
	}

/* Calculate the number of temporaries we can use for nQx and eQx and t3. */
/* Base our decision on the maximum amount of memory available.  Also, */
/* the main loop uses one temp for gg, so subtract one from numvals. */

	numvals = choose_pminus1_numvals (k, b, n, c, 1);
	numvals--;

/* Try various values of D until we find the best one */

	best_cost = 1e99;
	for (d = 5 * 2310; d >= 30; ) {

/* Try various values of E and using_t3 until we find the best one */

		for (i = 0; i < 10; i++) {
			int	using_t3;

/* Calculate e and the using_t3 flag.  Don't cost out e values where we */
/* don't have enough memory. */

			if (i <= 1) e = 1;
			else if (i <= 3) e = 2;
			else if (i <= 5) e = 4;
			else if (i <= 7) e = 6;
			else e = 12;
			using_t3 = i & 1;
			if (numvals <= e + 1 + using_t3) break;

/* Calculate the cost of this stage 2 plan */

			cost = cost_pminus1_plan (B, C, d, e, numvals, using_t3);

/* Reward higher E values because they should find more factors */
/* These rewards are close to a complete guess.  Little studying has */
/* been done on how often higher e values will find a factor. */

			if (e == 4) cost *= .95;
			if (e == 6) cost *= .90;
			if (e == 12) cost *= .85;

/* Remember best cost and best d and e */

			if (cost < best_cost) {
				best_cost = cost;
				D = d;
				E = e;
			}
		}

/* Try next smaller value of d */

		if (d > 2310) d = d - 2310;
		else if (d > 210) d = d - 210;
		else d = d - 30;
	}

/* Print out our selection */

	if (QA_IN_PROGRESS) {
		char	buf[200];
		sprintf (buf, "numvals = %d, D = %d, E = %d\n", numvals, D, E);
		OutputBoth (buf);
	}
}

/* Choose the best implementation of the pminus1 plan given the current */
/* memory settings.  We may decide to wait for more memory to be available. */
/* We may choose to use the t3 temporary variable. */

void choose_pminus1_implementation (
	double	k,
	unsigned long b,
	unsigned long n,
	signed long c,
	struct pm1_state *state,	/* P-1 state variable */
	int	ll_testing,		/* True if we are LL testing */
	int	*do_later,		/* Should we wait until more */
					/* memory is available? */
	int	*using_t3_result)	/* Should we use a temporary for t3 */
{
	unsigned long B, C, numvals;
	int	using_t3;
	double	cost, best_cost;

/* Copy some state variables for easier access */

	B = state->C_start;		/* Stage 2 starting point */
	C = state->C;			/* Stage 2 ending point */

/* Compute the number of values relatively prime to D */

	state->numrels = calc_numrels (D);

/* Assume we won't be able to do this pass now */

	*do_later = TRUE;

/* Calculate the number of temporaries we can use for nQx and eQx and t3. */
/* Base our decision on the maximum amount of memory available.  Also, */
/* the main loop uses one temp for gg, so subtract one from numvals. */

	numvals = choose_pminus1_numvals (k, b, n, c, 0);
	numvals--;

/* Check if we are only supposed to run stage 2 when the maximum amount */
/* memory is available */

	if (numvals != choose_pminus1_numvals (k, b, n, c, 1) - 1 &&
	    ll_testing &&
	    IniGetInt (INI_FILE, "OnlyRunStage2WithMaxMemory", 0))
		return;

/* If not much memory is available right now, try shrinking E so that */
/* we can make some progress now.  The 10 in the formula below is arbitrary. */

	if (E > 2 && numvals < E + 10) E = 2;

/* Try with and without using_t3 to find the best cost. */

	best_cost = 1e99;
	for (using_t3 = 0; using_t3 < 2; using_t3++) {

/* If there isn't enough memory to run this scenario, then do not */
/* bother costing it out. */

		if (numvals <= E + 1 + using_t3) break;

/* Calculate the cost of this stage 2 plan */

		cost = cost_pminus1_plan (B, C, D, E, numvals, using_t3);

/* Remember best cost and best d and e */

		if (cost < best_cost) {
			best_cost = cost;
			*do_later = FALSE;
			*using_t3_result = using_t3;
			state->rels_this_pass = numvals - (E+1) - using_t3;
		}
	}

/* Adjust rels_this_pass down if it is too high */

	if (state->rels_this_pass > state->numrels - state->rels_done)
		state->rels_this_pass = state->numrels - state->rels_done;
}

/* Fill the bit array in such a way that it maximizes prime pairings. */

int fill_pminus1_bitarray (
	struct pm1_state *state)	/* P-1 state variable */
{
	unsigned long adjusted_C_start, prime, jprime, clear;
	unsigned long m, stage2incr, first_m, pair, i;
	unsigned long *j, *jset;
	unsigned long relp[] = {7,11,13,17,19,23,29,31,37,41,43,47,0};

/* Allocate the bitarray, pad it so that the stage 2 bit testing loop */
/* does not examine unallocated memory */

	state->bitarray_len = (state->C + D + D + 15) >> 4;
	state->bitarray = (char *) malloc (state->bitarray_len);
	if (state->bitarray == NULL) {
		OutputStr ("Out of memory\n");
errexit:	state->bitarray_len = 0;
		free (state->bitarray);
		state->bitarray = NULL;
		return (FALSE);
	}
	memset (state->bitarray, 0, state->bitarray_len);

/* Set one bit for each prime between C_start and C */

	start_sieve (state->C_start);
	for (prime = sieve (); prime <= state->C; prime = sieve ()) {
		bitset (state->bitarray, prime >> 1);
		if (stopCheck ()) goto errexit;
	}

/* Now "move" some of the primes around so that we both maximize pairings */
/* and move the effective C_start higher.  We do this by testing 13*prime */
/* or 17*prime, etc. instead of prime (as long as 13*prime < C). */

	if (D >= 2310) adjusted_C_start = state->C / 13, jset = relp + 2;
	else if (D >= 210) adjusted_C_start = state->C / 11, jset = relp + 1;
	else adjusted_C_start = state->C / 7, jset = relp;
	stage2incr = (E == 1) ? D : D + D;
	first_m = (adjusted_C_start / D) * D;
	for (prime = state->C_start | 1; prime < adjusted_C_start; prime+=2) {
		if (!bittst (state->bitarray, prime >> 1)) continue;
		clear = prime;
		for (j = jset; *j; j++) {
			jprime = *j * prime;
			if (jprime > state->C) break;
			bitclr (state->bitarray, clear >> 1);
			bitset (state->bitarray, jprime >> 1);
			clear = jprime;
			if (jprime < adjusted_C_start) continue;

/* Test if jprime pairs up */

			if (E == 1) break;
			m = (jprime - first_m) / stage2incr * stage2incr + first_m + D;
			if (jprime < m) pair = m + (m - jprime);
			else pair = m - (jprime - m);
			if (bittst (state->bitarray, pair >> 1)) break;
		}
		if (stopCheck ()) goto errexit;
	}

/* Count the number of pairs.  This is used to calculate the stage 2 */
/* percent complete. */

	m = (adjusted_C_start < state->C_start) ? state->C_start : adjusted_C_start;
	m = (m / D + 1) * D;
	for (state->pairs_set = 0; state->C > m-D; m += stage2incr) {
	    for (i = 1; i < D; i += 2) {
		if (! bittst (state->bitarray, (m - i) >> 1) &&
		    (E == 1 || ! bittst (state->bitarray, (m + i) >> 1)))
			continue;
		state->pairs_set++;
	    }
	}
	state->pairs_done = 0;

/* All done */

	return (TRUE);
}

/* Recursively compute exp used in initial 3^exp calculation of a P-1 */
/* factoring run.  Don't forget to include 2*n in exp when factoring */
/* Mersenne numbers since factors must be 1 mod 2n */

void calc_exp (
	double	k,		/* K in K*B^N+C */
	unsigned long b,	/* B in K*B^N+C */
	unsigned long n,	/* N in K*B^N+C */
	signed long c,		/* C in K*B^N+C */
	giant	g,
	unsigned long B1,	/* P-1 stage 1 bound */
	unsigned long *p,
	unsigned long lower,
	unsigned long upper)
{
	unsigned long len;

/* Compute the number of result words we are to calculate */

	len = upper - lower;

/* Use recursion to compute the exponent.  This will perform better */
/* because mulg will be handling arguments of equal size. */

	if (len >= 50) {
		giant	x;
		calc_exp (k, b, n, c, g, B1, p, lower, lower + (len >> 1));
		x = newgiant (len + len);
		calc_exp (k, b, n, c, x, B1, p, lower + (len >> 1), upper);
		mulg (x, g);
		free (x);
		return;
	}

/* For Mersenne numbers, 2^n-1, make sure we include 2n in the calculated */
/* exponent (since factors are of the form 2kn+1).  For Fermat numbers, */
/* 2^n+1 (n is a power of 2), make sure the exponent is included in the */
/* calculated exponent as factors are of the form kn+1.  Otherwise, do */
/* nothing special -- start  with one. */

	if (lower == 0 && k == 1.0 && b == 2 && c == -1) itog (2*n, g);
	else if (lower == 0 && k == 1.0 && b == 2 && c == 1) itog (n, g);
	else setone (g);

/* Find all the primes in the range and use as many powers as possible */

	for ( ; *p <= B1 && (unsigned long) g->sign < len; *p = sieve ()) {
		unsigned long val, max;
		val = *p;
		max = B1 / *p;
		while (val <= max) val *= *p;
		ulmulg (val, g);
	}
}

/* Calculate how much of stage 2 is complete. */

double calc_stage2_pct (
	struct pm1_state *state)
{
	return ((double) state->pairs_done / (double) state->pairs_set);
}

/* Main P-1 entry point */

int pminus1 (
	double	k,		/* K in K*B^N+C */
	unsigned long b,	/* B in K*B^N+C */
	unsigned long n,	/* N in K*B^N+C */
	signed long c,		/* C in K*B^N+C */
	unsigned long B,	/* Stage 1 bound */
	unsigned long C_start,	/* Stage 2 starting point (usually equals B) */
	unsigned long C,	/* Stage 2 ending point */
	int	ll_testing)	/* Set to how_far_factored if we are */
				/* pre-factoring for a future LL test */
{
	struct pm1_state state;
	giant	exp;
	unsigned long stage_0_limit;
	unsigned long SQRT_B;
	unsigned long numrels, first_rel, last_rel;
	unsigned long i, j, m, stage2incr, prime, len, bit_number;
	unsigned long error_recovery_mode = 0;
	gwnum	x, gg, t3;
	char	filename[16], buf[100];
	int	res, retval, stage, escaped, saving, near_fft_limit, echk;
	long	write_time = DISK_WRITE_TIME * 60;
	time_t	start_time, current_time;
	double	pct, last_output, last_contact, last_output_r;
	int	using_t3;	/* Indicates we are using the gwnum t3 */
				/* to avoid a gwfftadd3 in stage 2 */
	int	do_later;	/* Set if low available memory suggests */
				/* we'd be better off waiting for nighttime */

/* Clear all timers */

restart:
	clear_timers ();
	escaped = FALSE;

/* Init filename */

	tempFileName (filename, n);
	strcat (filename, EXTENSION);
	filename[0] = 'm';
	if (c == 1) filename[0]--;

/* Get the current time */

	time (&start_time);

/* Choose a default value for the second bound if none was specified */

	if (C == 0) C = (B >= 42900000) ? 4290000000UL : B * 100;
	if (C < B) C = B;
	else if (B < 30) {
		OutputStr ("Using minimum bound #1 of 30\n");
		B = 30;
	}

/* Perform setup functions */

	memset (&state, 0, sizeof (state));
	res = pm1_setup1 (k, b, n, c);
	if (res) {
		sprintf (buf, "Cannot initialize FFT code, errcode=%d\n", res);
		OutputBoth (buf);
		return (FALSE);
	}
	last_contact = last_output = last_output_r = fft_count = 0;
	EXP_BEING_WORKED_ON = n;
	EXP_BEING_FACTORED = 1;

/* If we are near the maximum exponent this fft length can test, then we */
/* will roundoff check all multiplies */

	near_fft_limit = exponent_near_fft_limit ();
	gwsetnormroutine (0, ERRCHK || near_fft_limit, 0);

/* Compute the number we are factoring */

	if (!setN (k, b, n, c)) {
		pm1_cleanup ();
		return (FALSE);
	}

/* Output startup message */

	title ("P-1");
	if (!error_recovery_mode) {
		char	fft_desc[100];
		sprintf (buf, "P-1 on %s with B1=%lu, B2=%lu\n",
			 gwmodulo_as_string (), B, C);
		OutputStr (buf);
		gwfft_description (fft_desc);
		sprintf (buf, "Using %s\n", fft_desc);
		OutputStr (buf);
	}

/* Check for a continuation file */

	if (fileExists (filename)) {
		unsigned long processed;

/* Read in the save file.  If there was an error reading the file then */
/* restart the P-1 factoring job from scratch. */

		if (! pm1_restore (filename, &state, &processed, &x, &gg))
			goto from_scratch;

/* Handle stage 0 save files.  If the B values do not match, then use */
/* the bound given in the save file -- this may well result in an */
/* increased execution time if the saved B is larger than the B passed */
/* in to this routine. */

		if (state.stage == PM1_STAGE0) {
			bit_number = processed;
			goto restart0;
		}

/* To avoid an infinite loop of repeatable roundoff errors, we square */
/* the value read in from the P-1 save file.  This won't affect our final */
/* results, but will change the FFT data. */

		if (error_recovery_mode) {
			gwsquare (x);
			pm1_save (filename, &state, processed, x, gg);
			error_recovery_mode = 0;
		}

/* Handle stage 1 save files */

		if (state.stage == PM1_STAGE1) {
			if (B <= processed) {
				state.B_done = processed;
				state.C_done = processed;
				state.B = processed;
				goto restart2;
			}
			if (B < state.B) state.B = B;
			start_sieve (processed + 1);
			prime = sieve ();
			goto restart1;
		}

/* Handle stage 2 save files */

		if (state.stage == PM1_STAGE2) {

/* If B is larger than the one in the save file, then go back and */
/* do some more stage 1 processing.  Since this is very upsetting to */
/* an LL tester that has already begun stage 2 only do this for the */
/* non-LL tester. */

			if (B > state.B_done) {
				if (!ll_testing) {
					gwfree (gg);
					goto more_B;
				}
				B = state.B_done;
				C_start = state.B_done;
				sprintf (buf, "Ignoring suggested B1 value, using B1=%lu from the save file\n", B);
				OutputStr (buf);
			}

/* If we've already done enough stage 2, go do the stage 2 GCD */

			if (C <= state.C_done) {
				stage = 2;
				goto restart4;
			}

/* If we never really started stage 2, then do so now */

			if (state.bitarray_len == 0) goto more_C;

/* If LL testing and bound #2 has changed then use the original bound #2 */

			if (ll_testing && C > state.C) {
				C = state.C;
				sprintf (buf, "Ignoring suggested B2 value, using B2=%lu from the save file\n", C);
				OutputStr (buf);
			}

/* Resume stage 2 */

			goto restart3b;
		}

/* Handle case where we have a completed save file (the PM1_DONE state) */

		if (B > state.B_done) goto more_B;
		if (C > state.C_done) goto restart3a;

/* Note: if C_start != B then the user is using the undocumented feature */
/* of doing stage 2 in pieces.  Assume he knows what he is doing */

		if (C_start != B) {
			state.C_done = state.B_done;
			goto restart3a;
		}

/* The save file indicates we've tested to these bounds already */

		sprintf (buf, "%s already tested to B1=%lu and B2=%lu.\n",
			 gwmodulo_as_string (), state.B_done, state.C_done);
		OutputBoth (buf);
		goto done;
	}

/* Start this P-1 run from scratch starting with x = 3 */

from_scratch:
	bit_number = 0;
	x = gwalloc ();
	dbltogw (3.0, x);
	state.B_done = 0;
	state.B = B;

/* First restart point.  Compute the big exponent (a multiple of small */
/* primes).  Then compute 3^exponent.  The exponent always contains 2*p. */
/* We only compute 1.5 * B bits (up to 1.5 million).  The rest of the */
/* exponent will be done one prime at a time in the second part of stage 1. */
/* This stage uses 2 transforms per exponent bit. */

restart0:
	state.stage = PM1_STAGE0;
	start_timer (0);
	start_timer (1);
	start_sieve (2);
	prime = sieve ();
	stage_0_limit = (state.B > 1000000) ? 1000000 : state.B;
	i = ((unsigned long) (stage_0_limit * 1.5) >> 5) + 4;
	exp = newgiant (i << 1);
	calc_exp (k, b, n, c, exp, state.B, &prime, 0, i);

/* Find most significant bit and then ignore it */

	len = bitlen (exp);
	len--;

/* Now take the exponent and raise x to that power */

	gwsetmulbyconst (3);
	while (bit_number < len) {

/* Contact the server every now and then */

		if (fft_count > last_contact + 200) {
			pct = (double) bit_number / (double) len;
			if (prime < B) pct *= (double) prime / (double) B;
			if (B != C) EXP_PERCENT_COMPLETE = 0.5 - pct * 0.5;
			if (!communicateWithServer ()) escaped = 1;
			if (!pauseWhileRunning ()) escaped = 1;
			last_contact = fft_count;
		}

/* To avoid an infinite loop of repeatable roundoff errors, carefully */
/* get us past the offending iteration. */

		if (error_recovery_mode && bit_number == error_recovery_mode) {
			careful_iteration (x, NULL, n);
			if (bitval (exp, len - bit_number - 1)) {
				gwnum	three;
				three = gwalloc ();
				dbltogw (3.0, three);
				gwsetnormroutine (0, 0, 0);
				gwmul (three, x);
				gwfree (three);
			}
			error_recovery_mode = 0;
			saving = TRUE;
		}

/* Set various flags.  They control whether error-checking or the next FFT */
/* can be started. */

		else {
			escaped |= stopCheck ();
			saving = FALSE;
			echk = escaped || ERRCHK || near_fft_limit;
			if ((bit_number & 127) == 64) {
				time (&current_time);
				saving = (current_time - start_time > write_time);
				echk = 1;
			}

/* Either square x or square x and multiply it by three. */

			gwstartnextfft (!escaped && !saving &&
					!(QA_SAVE_FILES & 0x1) &&
					!error_recovery_mode &&
					bit_number+1 != len);
			if (bitval (exp, len - bit_number - 1)) {
				gwsetnormroutine (0, echk, 1);
				gwsquare (x);
			} else {
				gwsetnormroutine (0, echk, 0);
				gwsquare (x);
			}
		}

/* Test for an error */

		if (gw_test_for_error () || MAXERR >= 0.40625) goto error;
		bit_number++;

/* Every N squarings, output a progress report */

		if (ITER_OUTPUT != 999999999 &&
		    fft_count >= last_output + 2 * ITER_OUTPUT) {
			char	mask[80];
			pct = (double) bit_number / (double) len;
			if (prime < B) pct *= (double) prime / (double) B;
			pct = trunc_percent (pct * 100.0);
			sprintf (mask, "%%.%df%%%% P-1 stage 1", PRECISION);
			sprintf (buf, mask, pct);
			title (buf);
			sprintf (mask,
				 "%%s stage 1 is %%.%df%%%% complete. Time: ",
				 PRECISION);
			sprintf (buf, mask, gwmodulo_as_string (), pct);
			OutputTimeStamp ();
			OutputStr (buf);
			end_timer (0);
			print_timer (0, TIMER_NL | TIMER_OPT_CLR);
			start_timer (0);
			last_output = fft_count;
		}

/* Every N squarings, output a progress report to the results file */

		if ((ITER_OUTPUT_RES != 999999999 &&
		     fft_count >= last_output_r + 2 * ITER_OUTPUT_RES) ||
		    (NO_GUI && escaped)) {
			char	mask[80];
			pct = (double) bit_number / (double) len;
			if (prime < B) pct *= (double) prime / (double) B;
			pct = trunc_percent (pct * 100.0);
			sprintf (mask,
				 "%%s stage 1 is %%.%df%%%% complete.\n",
				 PRECISION);
			sprintf (buf, mask, gwmodulo_as_string (), pct);
			writeResults (buf);
			last_output_r = fft_count;
		}

/* Check for escape and/or if its time to write a save file */

		if (escaped || saving || (QA_SAVE_FILES & 0x1)) {
			pm1_save (filename, &state, bit_number, x, NULL);
			if (escaped) {
				free (exp);
				retval = FALSE;
				goto exit;
			}
			if (QA_SAVE_FILES & 0x1) {
				free (exp);
				QA_SAVE_FILES &= ~0x1;
				retval = QA_SAVE_TEST;
				goto exit;
			}
			time (&current_time);
			start_time = current_time;
		}

/* Use this opportunity to perform other miscellaneous tasks that may */
/* be required by this particular OS port */

		doMiscTasks ();
	}

/* If roundoff error recovery returned to restart0, but the roundoff error */
/* occurs after the above loop, then square the value and create a save file. */
/* This won't affect our final results, but will change the FFT data. */

	if (error_recovery_mode) {
		gwsquare (x);
		pm1_save (filename, &state, bit_number, x, NULL);
		error_recovery_mode = 0;
	}

/* Do stage 0 cleanup */

	gwsetnormroutine (0, ERRCHK || near_fft_limit, 0);
	free (exp);
	end_timer (0);
	end_timer (1);

/* This situation will probably never happen, but will handle it anyway */

	if (B > stage_0_limit && B < state.B) state.B = B;

/* Second restart point.  Do the larger primes of stage 1. */
/* This stage uses 2.5 transforms per exponent bit. */

restart1:
	start_timer (0);
	start_timer (1);
	state.stage = PM1_STAGE1;
	SQRT_B = (unsigned long) sqrt ((double) state.B);
	for ( ; prime <= state.B; prime = sieve ()) {

/* Apply as many powers of prime as long as prime^n <= B */

		if (prime > state.B_done) pm1_mul (x, prime);
		if (prime <= SQRT_B) {
			unsigned long mult, max;
			mult = prime;
			max = state.B / prime;
			for ( ; ; ) {
				mult *= prime;
				if (mult > state.B_done) pm1_mul (x, prime);
				if (mult > max) break;
			}
		}

/* Test for an error */

		if (gw_test_for_error () || MAXERR >= 0.40625) goto error;

/* Test for user interrupt */

		escaped = stopCheck ();

/* Contact the server every now and then */

		if (fft_count > last_contact + 200) {
			pct = (double) prime / (double) B;
			if (B != C) EXP_PERCENT_COMPLETE = pct * 0.5;
			if (!communicateWithServer ()) escaped = 1;
			if (!pauseWhileRunning ()) escaped = 1;
			last_contact = fft_count;
		}

/* Every N primes, output a progress report */

		if (ITER_OUTPUT != 999999999 &&
		    fft_count >= last_output + 2 * ITER_OUTPUT) {
			char	mask[80];
			pct = (double) prime / (double) B;
			pct = trunc_percent (pct * 100.0);
			sprintf (mask, "%%.%df%%%% P-1 stage 1", PRECISION);
			sprintf (buf, mask, pct);
			title (buf);
			sprintf (mask,
				 "%%s stage 1 is %%.%df%%%% complete. Time: ",
				 PRECISION);
			sprintf (buf, mask, gwmodulo_as_string (), pct);
			OutputTimeStamp ();
			OutputStr (buf);
			end_timer (0);
			print_timer (0, TIMER_NL | TIMER_OPT_CLR);
			start_timer (0);
			last_output = fft_count;
		}

/* Every N primes, output a progress report to the results file */

		if ((ITER_OUTPUT_RES != 999999999 &&
		     fft_count >= last_output_r + 2 * ITER_OUTPUT_RES) ||
		    (NO_GUI && escaped)) {
			char	mask[80];
			pct = (double) prime / (double) B;
			pct = trunc_percent (pct * 100.0);
			sprintf (mask,
				 "%%s stage 1 is %%.%df%%%% complete.\n",
				 PRECISION);
			sprintf (buf, mask, gwmodulo_as_string (), pct);
			writeResults (buf);
			last_output_r = fft_count;
		}

/* Check for escape and/or if its time to write a save file */

		time (&current_time);
		if (escaped || (QA_SAVE_FILES & 0x2) ||
		    current_time - start_time > write_time) {
			pm1_save (filename, &state, prime, x, NULL);
			if (escaped) {
				retval = FALSE;
				goto exit;
			}
			if (QA_SAVE_FILES & 0x2) {
				QA_SAVE_FILES &= ~0x2;
				retval = QA_SAVE_TEST;
				goto exit;
			}
			start_time = current_time;
		}

/* Use this opportunity to perform other miscellaneous tasks that may */
/* be required by this particular OS port */

		doMiscTasks ();
	}
	state.B_done = state.B;
	state.C_done = state.B;
	end_timer (0);
	end_timer (1);

/* Check for the rare case where we need to do even more stage 1 */
/* This happens when a save file was created with a smaller bound #1 */
/* than the bound #1 passed into this routine */

	if (B > state.B) {
more_B:		state.B = B;
		start_sieve (2);
		prime = sieve ();
		goto restart1;
	}

/* Stage 1 complete, print a message */

	sprintf (buf, "%s stage 1 complete. %.0f transforms. Time: ",
		 gwmodulo_as_string (), fft_count);
	OutputStr (buf);
	print_timer (1, TIMER_NL | TIMER_CLR);
	clear_timers ();
	last_contact = last_output = last_output_r = fft_count = 0;

/* Print out round off error */

	if (ERRCHK) {
		sprintf (buf, "Round off: %.10g\n", MAXERR);
		OutputStr (buf);
		MAXERR = 0.0;
	}

/* Check to see if we found a factor - do GCD (x-1, N) */

restart2:
	if (C <= B ||
	    (!QA_IN_PROGRESS && IniGetInt (INI_FILE, "Stage1GCD", 1))) {
		if (ll_testing)
			OutputStr ("Starting stage 1 GCD - please be patient.\n");
		start_timer (0);
		gwaddsmall (x, -1);
		retval = gcd (x);
		gwaddsmall (x, 1);
		if (! retval) {
			pm1_save (filename, &state, B, x, NULL);
			goto exit;
		}
		end_timer (0);
		OutputStr ("Stage 1 GCD complete. Time: ");
		print_timer (0, TIMER_NL | TIMER_CLR);
		stage = 1;
		if (FAC != NULL) goto bingo;
	}

/* Skip second stage if so requested */

	if (C <= B) goto msg_and_exit;

/*
   Stage 2:  Use ideas from Crandall, Zimmermann, and Montgomery on each
   prime below C.  This code is more efficient the more memory you can
   give it.
   x: value at the end of stage 1 
*/

/* Initialize variables for second stage.  We set gg to x-1 in case the */
/* user opted to skip the GCD after stage 1. */

restart3a:
	gg = gwalloc ();
	gwcopy (x, gg);
	gwaddsmall (gg, -1);
	state.stage = PM1_STAGE2;

/* Choose a good value for D and E - one that reduces the number of */
/* multiplications, yet doesn't use too much memory.  This plan is based */
/* on the maximum available memory.  Once the plan is selected it cannot */
/* be changed in multi-pass stage 2 runs. */

more_C:	state.C_start = (C_start > state.C_done) ? C_start : state.C_done;
	state.C = C;
	choose_pminus1_plan (k, b, n, c, &state);
	if (D == 0) {	/* We'll never have enough memory for stage 2 */
		OutputStr ("Insufficient memory to ever run stage 2.\n");
		C = state.B_done;
		goto restart4;
	}
	if (! fill_pminus1_bitarray (&state)) {
		pm1_save (filename, &state, 0, x, gg);
		retval = FALSE;
		goto exit;
	}
	state.rels_done = 0;

/* Restart here when in the middle of pass stage 2 or */
/* Move to the next pass of a multi-pass stage 2 run */

restart3b:

/* Choose the best plan implementation given the currently available memory. */
/* This implementation could be anything from "wait until we have more */
/* memory" to deciding whether using_t3 should be set. */

	choose_pminus1_implementation (k, b, n, c, &state, ll_testing, &do_later, &using_t3);
	if (do_later || (QA_SAVE_FILES & 0x4)) {
		pm1_save (filename, &state, 0, x, gg);
		if (QA_SAVE_FILES & 0x4) {
			QA_SAVE_FILES &= ~0x4;
			retval = QA_SAVE_TEST;
			goto exit;
		}
		goto nomem;
	}

/* Here is where we restart the next pass of a multi-pass stage 2 */

	start_timer (0);
	start_timer (1);

/* On first pass, allocate P-1 stage 2 memory */

	if (nQx == NULL) {
		nQx = (gwnum *) malloc ((D>>1) * sizeof (gwnum));
		if (nQx == NULL) goto lowmem;
		eQx = (gwnum *) malloc ((E+1) * sizeof (gwnum));
		if (eQx == NULL) goto lowmem;
	}

/* Clear the nQx array for this pass */

	memset (nQx, 0, (D>>1) * sizeof (gwnum));

/* Compute x^(1^e), x^(3^e), ..., x^((D-1)^e) */

	for (i = 1, j = 0; ; i += 2) {
		if (! relatively_prime (i, D)) continue;
		if (++j > state.rels_done) break;
	}
	first_rel = i;
	gwfft (x, x);				/* fd_init requires fft of x */
	fd_init (i, 2, x);
	for (numrels = 0; ; ) {			/* Compute x^(i^e) */
		if (relatively_prime (i, D)) {
			j = (i - 1) >> 1;
			nQx[j] = gwalloc ();
			if (nQx[j] == NULL) {
				gwfftfftmul (x, x, x);	/* Unfft x for save */
				goto lowmem;
			}
			gwcopy (eQx[0], nQx[j]);
			numrels++;
			last_rel = i;
		}
		i = i + 2;
		if (i >= D) break;
		if (numrels == state.rels_this_pass) break;
		fd_next ();
		if (gw_test_for_error () || MAXERR >= 0.40625) goto error;
		if (stopCheck () || (QA_SAVE_FILES & 0x8)) {
			fd_term ();
			gwfftfftmul (x, x, x);	/* Unfft x - generates x^2 */
			pm1_save (filename, &state, 0, x, gg);
			if (QA_SAVE_FILES & 0x8) {
				QA_SAVE_FILES &= ~0x8;
				retval = QA_SAVE_TEST;
				goto exit;
			}
			retval = FALSE;
			goto exit;
		}
	}
	fd_term ();

/* Compute m = CEIL(start/D)*D, the first group we work on in stage 2 */

	if (D >= 2310) m = state.C / 13;
	else if (D >= 210) m = state.C / 11;
	else m = state.C / 7;
	if (m < state.C_start) m = state.C_start;
	m = (m / D + 1) * D;
	stage2incr = (E == 1) ? D : D + D;

/* Scan the bit array until we find the first group with a bit set. */
/* When continuing from a save file there could be many groups that */
/* have already been completed. */

	for ( ; state.C > m-D; m += stage2incr) {
	    for (i = first_rel; i <= last_rel; i += 2) {
		if (nQx[i>>1] == NULL) continue;
		if (! bittst (state.bitarray, (m - i) >> 1) &&
		    (E == 1 || ! bittst (state.bitarray, (m + i) >> 1)))
			continue;
		goto found_a_bit;
	    }
	}
found_a_bit:;

/* Initialize for computing successive x^(m^e) */

	fd_init (m, stage2incr, x);

/* Unfft x for use in save files.  Actually this generates x^2 which */
/* is just fine - no stage 2 factors will be missed (in fact it could */
/* find more factors) */

	gwfftfftmul (x, x, x);

/* Now touch all the nQx and eQx values so that when gg is used, x is */
/* swapped out rather than a value we will need in the near future. */
/* In other words, make the gwnum x the least-recently-used. */

	for (i = 0; i <= E; i++) gwtouch (eQx[i]);
	for (i = first_rel; i < last_rel; i += 2) {
		j = i >> 1;
		if (nQx[j] != NULL) gwtouch (nQx[j]);
	}

/* Now do a pass of stage 2 */

	stage = 2;

/* When E >= 2, we can do prime pairing and each loop iteration */
/* handles the range m-D to m+D.  When E = 1, each iteration handles */
/* the range m-D to m. */

	if (using_t3) {
		t3 = gwalloc ();
		if (t3 == NULL) goto lowmem;
	}
	for ( ; state.C > m-D; m += stage2incr) {
	    int	inner_loop_done = FALSE;
	    int	last_pass = (m + stage2incr - D >= state.C);
	    time (&current_time);
	    saving = (current_time - start_time > write_time ||
		      (QA_SAVE_FILES & 0x10));

/* Test all the relprimes between m-D and m */

	    for (i = first_rel; ; i += 2) {

/* Move onto the next m value when we are done with all the relprimes */

		if (i > last_rel) {	/* Compute next x^(m^e) */
			if (!last_pass) fd_next ();
			inner_loop_done = TRUE;
			escaped = stopCheck ();
			goto errchk;
		}

/* Skip this relprime if we aren't processing it this pass */ 

		j = i >> 1;
		if (nQx[j] == NULL) continue;

/* Skip this relprime if neither m - i nor its pair m + i are set */
/* in the bitarray. */

		if (! bittst (state.bitarray, (m - i) >> 1) &&
		    (E == 1 || ! bittst (state.bitarray, (m + i) >> 1)))
			continue;

/* Mul this eQx - nQx value into gg */

		escaped = stopCheck ();
		gwstartnextfft (!escaped && !saving);
		if (using_t3) {
			gwfftsub3 (eQx[0], nQx[j], t3);
			gwfftmul (t3, gg);
		} else {
			gwfftsub3 (eQx[0], nQx[j], eQx[0]);
			gwfftmul (eQx[0], gg);
			gwfftadd3 (eQx[0], nQx[j], eQx[0]);
		}

/* Clear this bit or bits in case a save file is written.  Subtract one */
/* from the bit count so that calc_stage2_pct is accurate. */

		bitclr (state.bitarray, (m - i) >> 1);
		if (E >= 2) bitclr (state.bitarray, (m + i) >> 1);
		state.pairs_done++;

/* Test for errors */

errchk:		if (gw_test_for_error () || MAXERR >= 0.40625) goto error;

/* Contact the server every now and then */

		if (fft_count > last_contact + 200) {
			pct = calc_stage2_pct (&state);
			EXP_PERCENT_COMPLETE = 0.5 + pct * 0.5;
			if (!communicateWithServer ()) escaped = 1;
			if (!pauseWhileRunning ()) escaped = 1;
			last_contact = fft_count;
		}

/* Write out a message every now and then */

		if (ITER_OUTPUT != 999999999 &&
		    fft_count >= last_output + 2 * ITER_OUTPUT) {
			char	mask[80];
			pct = calc_stage2_pct (&state);
			pct = trunc_percent (pct * 100.0);
			sprintf (mask, "%%.%df%%%% P-1 stage 2", PRECISION);
			sprintf (buf, mask, pct);
			title (buf);
			sprintf (mask,
				"%%s stage 2 is %%.%df%%%% complete. Time: ",
				PRECISION);
			sprintf (buf, mask, gwmodulo_as_string (), pct);
			OutputTimeStamp ();
			OutputStr (buf);
			end_timer (0);
			print_timer (0, TIMER_NL | TIMER_OPT_CLR);
			start_timer (0);
			last_output = fft_count;
		}

/* Write out a message to the results file every now and then */

		if ((ITER_OUTPUT_RES != 999999999 &&
		     fft_count >= last_output_r + 2 * ITER_OUTPUT_RES) ||
		    (NO_GUI && escaped)) {
			char	mask[80];
			pct = calc_stage2_pct (&state);
			pct = trunc_percent (pct * 100.0);
			sprintf (mask,
				"%%s stage 2 is %%.%df%%%% complete.\n",
				PRECISION);
			sprintf (buf, mask, gwmodulo_as_string (), pct);
			writeResults (buf);
			last_output_r = fft_count;
		}

/* Periodicly write a save file.  If we escaped, free eQx memory so */
/* that pm1_save can reuse it to convert x and gg to binary.  If we */
/* have been using t3 as a temporary, free that for the same reason. */
/* "Touch" gg so that in low memory situations, the reading in of x */
/* swaps out one of the eQx or nQx values rather than gg. */

		if (escaped || saving) {
			if (escaped) fd_term ();
			if (using_t3) gwfree (t3);
			gwtouch (gg);
			pm1_save (filename, &state, 0, x, gg);
			if (escaped) {
				retval = FALSE;
				goto exit;
			}
			if (QA_SAVE_FILES & 0x10) {
				QA_SAVE_FILES &= ~0x10;
				retval = QA_SAVE_TEST;
				goto exit;
			}
			start_time = current_time;
			saving = FALSE;
			if (using_t3) t3 = gwalloc ();
		}

/* Use this opportunity to perform other miscellaneous tasks that may */
/* be required by this particular OS port */

		doMiscTasks ();

/* Leave inner loop to work on the next m value */

		if (inner_loop_done) break;
	    }
	}
	if (using_t3) gwfree (t3);
	fd_term ();

/* Free up the nQx values for the next pass */

	for (i = first_rel; i <= last_rel; i += 2) {
		j = i >> 1;
		if (nQx[j] != NULL) gwfree (nQx[j]);
	}

/* Since we set gwstartnextfft above, we must do another harmless squaring */
/* here to make sure gg has not been partially FFTed.  We cannot convert gg */
/* to an integer for GCD if it has been partially FFTed. */

	gwstartnextfft (FALSE);
	gwsquare (gg);

/* Check to see if another pass is required */

	end_timer (0);
	end_timer (1);
	state.rels_done += state.rels_this_pass;
	if (state.rels_done < state.numrels) goto restart3b;
	free (state.bitarray);
	state.bitarray = NULL;
	state.bitarray_len = 0;

/* Check for the rare case where we need to do even more stage 2. */
/* This happens when a save file was created with a smaller bound #2 */
/* than the bound #2 passed into this routine. */

	state.C_done = state.C;
	if (C > state.C_done) goto more_C;

/* Stage 2 is complete */

	sprintf (buf, "%s stage 2 complete. %.0f transforms. Time: ",
		 gwmodulo_as_string (), fft_count);
	OutputStr (buf);
	print_timer (1, TIMER_NL | TIMER_CLR);
	clear_timers ();

/* Print out round off error */

	if (ERRCHK) {
		sprintf (buf, "Round off: %.10g\n", MAXERR);
		OutputStr (buf);
		MAXERR = 0.0;
	}

/* See if we got lucky! */

restart4:
	if (ll_testing)
		OutputStr ("Starting stage 2 GCD - please be patient.\n");
	start_timer (0);
	if ((QA_SAVE_FILES & 0x20) || !gcd (gg)) {
		pm1_save (filename, &state, C, x, gg);
		if (QA_SAVE_FILES & 0x20) {
			QA_SAVE_FILES &= ~0x20;
			retval = QA_SAVE_TEST;
			goto exit;
		}
		retval = FALSE;
		goto exit;
	}
	state.stage = PM1_DONE;
	end_timer (0);
	OutputStr ("Stage 2 GCD complete. Time: ");
	print_timer (0, TIMER_NL | TIMER_CLR);
	if (FAC != NULL) goto bingo;

/* Output line to results file indicating P-1 run */

msg_and_exit:
	sprintf (buf, "%s completed P-1, B1=%lu", gwmodulo_as_string (), B);
	if (C > B) {
		if (E <= 2)
			sprintf (buf+strlen(buf), ", B2=%lu", C);
		else
			sprintf (buf+strlen(buf), ", B2=%lu, E=%lu", C, E);
	}
	sprintf (buf+strlen(buf), ", Wc%d: %08lX\n", PORT, SEC5 (n, B, C));
	writeResults (buf);
	if (ll_testing) spoolMessage (PRIMENET_RESULT_MESSAGE, buf);

/* Send a kludgy message to server that P-1 factoring is complete. */
/* We add 0.5 to how_far_factored to indicate this state. */

	{
		struct primenetAssignmentResult pkt;
		memset (&pkt, 0, sizeof (pkt));
		pkt.exponent = n;
		pkt.resultType = PRIMENET_RESULT_NOFACTOR;
		pkt.resultInfo.how_far_factored = (double) ll_testing + 0.5;
		spoolMessage (PRIMENET_ASSIGNMENT_RESULT, &pkt);
	}

/* Create save file so that we can expand bound 1 or bound 2 at a later date */
/* If this is pre-factoring for an LL test, then delete the large save file */

	if (ll_testing)
		_unlink (filename);
	else
		pm1_save (filename, &state, 0, x, NULL);

/* Return sucsessful completion */

done:	retval = TRUE;

/* Update the exponent in the work-to-do-file */

	updateWorkToDo (k, b, n, c, WORK_PMINUS1, 0);

/* Free memory and return */

exit:	pm1_cleanup ();
	return (retval);

/* Low on memory, reduce memory settings and try again */

lowmem:	fd_term ();
	pm1_save (filename, &state, 0, x, gg);
	pm1_cleanup ();
	i = (unsigned int) (0.80 * max_mem ());
	if (DAY_MEMORY >= i) DAY_MEMORY = i;
	if (NIGHT_MEMORY >= i) NIGHT_MEMORY = i;
	sprintf (buf, "Memory allocation error.  Day and night memory settings changed to %dMB and %dMB.\n",
		 DAY_MEMORY, NIGHT_MEMORY);
	OutputBoth (buf);
	goto restart;

/* Return stop code indicating we don't have enough memory right now */

nomem:	OutputStr ("Not enough memory available to run stage 2 now.\n");
	OutputStr ("Will try again at a later time.\n");
	STOP_REASON = STOP_NOT_ENOUGH_MEM;
	retval = FALSE;
	goto exit;

/* Print a message if we found a factor! */

bingo:	if (stage == 1)
		sprintf (buf, "P-1 found a factor in stage #1, B1=%lu.\n", B);
	else
		sprintf (buf,
			 "P-1 found a factor in stage #2, B1=%lu, B2=%lu.\n",
			 B, C);
	writeResults (buf);
	printFactor (k, b, n, c);
	if (isone (FAC)) goto msg_and_exit;
	if (QA_IN_PROGRESS) {
		_unlink (filename);
		retval = TRUE;
		goto exit;
	}

/* If LL testing, output result to the server and free all save files */

	if (ll_testing) {
		struct primenetAssignmentResult pkt;
		memset (&pkt, 0, sizeof (pkt));
		pkt.exponent = n;
		pkt.resultType = PRIMENET_RESULT_FACTOR;
		gtoc (FAC, pkt.resultInfo.factor,
		      sizeof (pkt.resultInfo.factor));
		spoolMessage (PRIMENET_ASSIGNMENT_RESULT, &pkt);
		_unlink (filename);
		filename[0] = 'p';
		_unlink (filename);
		filename[0] = 'q';
		_unlink (filename);
	}

/* Otherwise create save file so that we can expand bound 1 or bound 2 */
/* at a later date. */

	else
		pm1_save (filename, &state, 0, x, NULL);

/* Remove the exponent from the worktodo.ini file */

	updateWorkToDo (k, b, n, c, WORK_FACTOR, 0);

/* Cleanup and return */

	free (FAC);
	FAC = NULL;
	retval = TRUE;
	goto exit;

/* Output an error message saying we are restarting. */
/* Sleep five minutes before restarting from last save file. */
/* Once errors start occurring, decrease the disk write time */

error:	pm1_cleanup ();
	if (near_fft_limit && MAXERR >= 0.40625) {
		sprintf (buf, "Possible roundoff error (%.8g), backtracking to last save file.\n", MAXERR);
		OutputStr (buf);
	} else {
		OutputBoth ("SUMOUT error occurred.\n");
		if (! SleepFive ()) return (FALSE);
	}
	error_recovery_mode = bit_number ? bit_number : 1;
	write_time >>= 1;
	if (write_time < 300) write_time = 300;
	goto restart;
}

/* Read a file of P-1 tests to run as part of a QA process */
/* The format of this file is: */
/*	k, n, c, B1, B2_start, B2_end, factor */
/* Use Advanced/Time 9992 to run the QA suite */

int pminus1_QA ()
{
	FILE	*fd;
	int	savefiles;

/* Set the title */

	title ("QA");

/* Open QA file */

	fd = fopen ("qa_pm1", "r");
	if (fd == NULL) {
		OutputStr ("File named 'qa_pm1' could not be opened.\n");
		return (TRUE);
	}

/* Loop until the entire file is processed */

	for ( ; ; ) {
		double	k;
		unsigned long b, n, B1, B2_start, B2_end;
		signed long c;
		char	fac_str[80];
		giant	f;
		int	retval, success;

/* Read a line from the file */

		n = 0;
		fscanf (fd, "%f,%lu,%lu,%ld,%lu,%lu,%lu,%s\n",
			&k, &b, &n, &c, &B1, &B2_start, &B2_end, &fac_str);
		if (n == 0) break;

/* If p is 1, set QA_TYPE */

		if (n == 1) {
			QA_TYPE = c;
			savefiles = B1;
			continue;
		}

/* Convert the factor we expect to find into a "giant" type */

		f = newgiant (strlen (fac_str));
		ctog (fac_str, f);

/*test various num_tmps
test 4 (or more?) stage 2 code paths
print out each test case (all relevant data)*/

/* Set some global variables indicating QA is in progress. */

		QA_IN_PROGRESS = TRUE;
		QA_SAVE_FILES = savefiles;

/* Do the P-1 */

		if (B2_start < B1) B2_start = B1;
		do {
			retval = pminus1 (k, b, n, c, B1, B2_start, B2_end, FALSE);
			if (!retval) {
				QA_IN_PROGRESS = FALSE;
				QA_TYPE = 0;
				QA_SAVE_FILES = 0;
				fclose (fd);
				return (FALSE);
			}
		} while (retval == QA_SAVE_TEST);

/* See if we found the factor */

		if (FAC == NULL) success = FALSE;
		else {
			modg (f, FAC);
			success = isZero (FAC);
			free (FAC);
			FAC = NULL;
		}
		if (!success) {
			char	buf[200];
			sprintf (buf, "ERROR: Factor not found. Expected %s\n", fac_str);
			OutputBoth (buf);
		}
		free (f);
	}

/* Cleanup */

	QA_IN_PROGRESS = FALSE;
	QA_TYPE = 0;
	QA_SAVE_FILES = 0;
	fclose (fd);
	return (TRUE);
}
