#define ALL_BENCH		1//0
#define UNIQUE			6*8//0
#define UNIQUE_SMALL_FIRST	1//0
#define CHECK_OFTEN		0
#define MIN_K_BITS		1
#define MAX_K_BITS		5//49
#define MAX_C_BITS_FOR_SMALL_K	4
#define MAX_C_BITS_FOR_LARGE_K	3//20
#define MAX_B			2
#define MIN_N			300000//300
#define MAX_N			500000000//12000000

/* Some gwnum.c internal definitions */

extern "C" double bit_length;		/* Bit length of k*b^n */
extern "C" double fft_max_bits_per_word;/* Maximum bits per data word that */
					/* this FFT size can support */
extern "C" double virtual_bits_per_word ();
extern "C" int gwinfo (double, unsigned long, unsigned long, signed long,
		       unsigned long);
extern "C" unsigned long* INFT[4];

/* Unique FFT size code.  Return TRUE if FFT has been done before. */

static	char	unique_minus[400] = {0};
static	char	unique_plus[400] = {0};
static	int	unique_count = UNIQUE;

int unique_tested (unsigned long fftlen, signed long c)
{
	int	i, count;
	char	*arr;

	if (UNIQUE == 0) return (FALSE);

	if (UNIQUE_SMALL_FIRST &&
	    fftlen > 4000000 && unique_count > UNIQUE / 2) return (TRUE);

	for (count = 0; (fftlen & 1) == 0; count++) fftlen >>= 1;

	arr = (c > 0) ? (char *) unique_plus : (char *) unique_minus;
	if (arr[count * 8 + fftlen]) return (TRUE);
	arr[count * 8 + fftlen] = 1;  unique_count--;
	return (FALSE);
}

/* Generate a random k of the given size */

double gen_k (int bits)
{
	double	k;
	int	i;

/* To reduce number of generic reductions due to gcd (k,c) > 1, only */
/* return small k values not divisible by 3, 5, and 7 */

	do {
		k = 1.0;
		for (i = bits; --i; )
			k = k * 2.0 + ((i == 1) ? 1 : (rand () & 1));
	} while (k > 7.0 && k < 4000000000.0 &&
		 ((unsigned long) k % 3 == 0 ||
		  (unsigned long) k % 5 == 0 ||
		  (unsigned long) k % 7 == 0));
	return (k);
}

/* Generate a random c of the given size */

int gen_c (int bits)
{
	int	c, i;

/* To reduce number of generic reductions due to gcd (k,c) > 1, only */
/* return c values not divisible by 3, 5, and 7 */

	do {
		c = 1;
		for (i = bits; --i; )
			c = c * 2 + ((i == 1) ? 1 : (rand () & 1));
	} while (c > 7 && (c % 3 == 0 || c % 5 == 0 || c % 7 == 0));
	if (rand () & 1) c = -c;
	return (c);
}

/* Generate a random n of the given size */

unsigned int gen_n ()
{
	unsigned int n, nbits;
static	unsigned int min_n_bits = 0;
static	unsigned int max_n_bits = 0;
static	unsigned int range_n_bits = 0;

	if (max_n_bits == 0) {
		min_n_bits = (unsigned int)
			floor (log ((double) MIN_N) / log (2.0));
		max_n_bits = (unsigned int)
			ceil (log ((double) MAX_N) / log (2.0));
		range_n_bits = max_n_bits - min_n_bits + 1;
	}
	do {
		nbits = (unsigned int) rand () % range_n_bits + min_n_bits;
		n = 1;
		while (--nbits) n = n * 2 + (rand () & 1);
	} while (n < MIN_N || n > MAX_N);
	return (n);
}


/* Generate random data to start the test */

void gen_data (gwnum x, giant g)
{
	unsigned long i, len;

/* Generate the random number */

	len = (((unsigned long) bit_length) >> 5) + 1;
	for (i = 0; i < len; i++) {
		g->n[i] = ((unsigned long) rand() << 20) +
			  ((unsigned long) rand() << 10) +
			  (unsigned long) rand();
	}
	g->sign = len;
	specialmodg (g);
	gianttogw (g, x);
}

/* Set and print random seed */

void set_seed ()
{
	int	seed;
	char	buf[100];

	seed = time (NULL);
	srand (seed);
	sprintf (buf, "Random seed is %ld\n", seed);
	OutputBoth (buf);
}

/* Compare gwnum to a giant */

void compare (gwnum x, giant g)
{
	giant	tmp;

	tmp = popg ((PARG >> 4) + 13);
	gwtogiant (x, tmp);
	if (gcompg (g, tmp) != 0)
		OutputBoth ("Test failed.\n");
	pushg (1);
}

/* Thoroughly test the current setup.  This is just like test_it except */
/* that rather than using giants code to test the results, we compare */
/* the final result to every possible FFT implementation for this FFT */
/* length.  This is useful in that more FFT code is tested, the tests are */
/* faster -- especially in the really large FFT cases where giants is not */
/* practical. */

void test_it_all (double k, unsigned long b, unsigned long n, signed long c)
{
	unsigned long *info, fftlen;
	gwnum	x, x2, x3, x4;
	giant	g, g2, g3;
	int	i, ii, res;
	double	diff, maxdiff;
	char	buf[200];

	k = 1.0;
	c = (c < 0) ? -1 : 1;
	g = g2 = g3 = NULL;

/* Find out what fftlen we would normally use */

	res = gwinfo (k, b, n, c, 0);
	if (res == 1) return;		// Don't do zero pad FFTs
	if (res) {
		sprintf (buf,
			 "Cannot QA k=%.0f, b=%ld, n=%ld, c=%ld: res=%d\n",
			 k, b, n, c, res);
		OutputBoth (buf);
		return;
	}
	info = INFT[0];
	fftlen = info[1];

	if (fftlen < 160000) return;	// Only large FFTs have multiple impls.

	if (unique_tested (fftlen, c)) return;	// Only test new fft lengths
	
	sprintf (buf, "QA of k=%.0f, b=%ld, n=%ld, c=%ld with FFTlen=%ldK\n",
		 k, b, n, c, fftlen >> 10);
	OutputBoth (buf);

/* Loop over all possible FFT implementations */

	for (ii = 1; ; ii++) {
		bench_pick_nth_fft = ii;
		res = gwsetup (k, b, n, c, fftlen);
		bench_pick_nth_fft = 0;
		if (res) break;
		if (GENERAL_MOD) {	// Can't happen
			gwdone ();	// BUG - we ought to undo unique bit
			return;		// Don't do generic mod FFTs
		}

/* Alloc and init numbers */

		x = gwalloc ();
		if (x == NULL) goto nomem;
		x2 = gwalloc ();
		if (x2 == NULL) goto nomem;
		x3 = gwalloc ();
		if (x3 == NULL) goto nomem;
		x4 = gwalloc ();
		if (x4 == NULL) goto nomem;

		if (ii == 1) {
			g = newgiant (((unsigned long) bit_length >> 4) + 20);
			if (g == NULL) goto nomem;
			g2 = newgiant (((unsigned long) bit_length >> 4) + 20);
			if (g2 == NULL) goto nomem;
			g3 = newgiant (((unsigned long) bit_length >> 4) + 20);
			if (g3 == NULL) goto nomem;
			gen_data (x, g);
		} else
			gianttogw (g, x);

/* Test 50 squarings */	

		gwcopy (x, x2);
		maxdiff = 0.0;
		gwsetnormroutine (0, 1, 0);	/* Enable error checking */
		for (i = 0; i < 50; i++) {
			gwstartnextfft ((i & 3) == 2);
						/* Test POSTFFT sometimes */
						/* Test gwsetaddin without */
						/* and with POSTFFT set */
			if ((i == 45 || i == 46) && abs (CARG) == 1)
				gwsetaddin (-31);
			gwsquare (x);
			diff = fabs (gwsuminp (x) - gwsumout (x));
			if (diff > maxdiff) maxdiff = diff;
			if ((i == 45 || i == 46) && abs (CARG) == 1)
				gwsetaddin (0);
		}
		if (MAXDIFF < 1e50)
			sprintf (buf, "Squares complete. MaxErr=%.8g, SumoutDiff=%.8g/%.8g(%d to 1)\n", MAXERR, maxdiff, MAXDIFF, (int) (MAXDIFF / maxdiff));
		else
			sprintf (buf, "Squares complete. MaxErr=%.10g\n", MAXERR);
		OutputBoth (buf);

/* Test mul by const */

		gwsetmulbyconst (3);
		gwsetnormroutine (0, 1, 1);
		gwsquare (x);
		gwsetnormroutine (0, 1, 0);
		diff = fabs (gwsuminp (x) - gwsumout (x));
		if (diff > maxdiff) maxdiff = diff;

/* Test square carefully */

		if (abs (CARG) == 1) gwsetaddin (-42);
		gwsquare_carefully (x);
		diff = fabs (gwsuminp (x) - gwsumout (x));
		if (diff > maxdiff) maxdiff = diff;
		if (abs (CARG) == 1) gwsetaddin (0);

/* Test gwaddquick, gwsubquick */

		gwadd3quick (x, x2, x3);
		gwsub3quick (x, x2, x4);

/* Test gwadd and gwsub */

		gwadd (x, x); gwadd (x, x); gwadd (x, x);
		gwsub (x3, x);
		gwadd (x4, x);
		gwadd3 (x3, x4, x2);
		gwsub3 (x3, x, x4);
		gwadd (x2, x);
		gwadd (x4, x);

/* Test gwaddsub */

		gwaddsub (x, x2);	// compute x+x2 and x-x2
		gwaddsub4 (x, x2, x3, x4);	// compute x+x2 and x-x2
		gwadd (x2, x);
		gwadd (x3, x);
		gwadd (x4, x);

/* Test gwaddsmall */

		gwaddsmall (x, 111);

/* Do some multiplies to make sure that the adds and subtracts above */
/* normalized properly. */

		gwfft (x, x);
		gwfftfftmul (x, x, x);
		diff = fabs (gwsuminp (x) - gwsumout (x));
		if (diff > maxdiff) maxdiff = diff;

		gwfft (x, x2); gwcopy (x2, x); gwfftadd3 (x, x2, x4);
		gwfftmul (x4, x3);
		diff = fabs (gwsuminp (x3) - gwsumout (x3));
		if (diff > maxdiff) maxdiff = diff;
		gwfft (x3, x4);
		gwfftfftmul (x4, x2, x);
		diff = fabs (gwsuminp (x) - gwsumout (x));
		if (diff > maxdiff) maxdiff = diff;

/* Print final stats */

		if (GWERROR) OutputBoth ("GWERROR set during calculations.\n");
		if (maxdiff > MAXDIFF) OutputBoth ("Sumout failed during test.\n");
		if (MAXDIFF < 1e50)
			sprintf (buf, "Test complete. MaxErr=%.8g, SumoutDiff=%.8g/%.8g(%d to 1)\n", MAXERR, maxdiff, MAXDIFF, (int) (MAXDIFF / maxdiff));
		else
			sprintf (buf, "Test complete. MaxErr=%.10g\n", MAXERR);
		OutputBoth (buf);

/* Do the final compare */

		if (ii == 1) {
			gwtogiant (x, g2);
		} else {
			gwtogiant (x, g3);
			if (gcompg (g2, g3)) {
				strcpy (buf, "Mismatched result.\n");
				OutputBoth (buf);
			} else {
				strcpy (buf, "Results match!\n");
				OutputBoth (buf);
			}
		}
		OutputBoth ("\n");

/* Do next FFT implementation */

		gwdone ();
	}
bye:	free (g);
	free (g2);
	free (g3);
	return;

nomem:	OutputBoth ("Out of memory\n");
	gwdone ();
	goto bye;
}


/* Thoroughly test the current setup */

void test_it (int square_only)
{
	gwnum	x, x2, x3, x4;
	giant	g, g2, g3, g4;
	int	i;
	double	diff, maxdiff = 0.0;
	char	buf[200];

/* Init random number generator */

//	set_seed ();

/* Alloc and init numbers */

	x = gwalloc ();
	x2 = gwalloc ();
	x3 = gwalloc ();
	x4 = gwalloc ();
	g = popg (((unsigned long) bit_length >> 4) + 20);
	g2 = popg (((unsigned long) bit_length >> 4) + 20);
	g3 = popg (((unsigned long) bit_length >> 4) + 20);
	g4 = popg (((unsigned long) bit_length >> 4) + 20);
	gen_data (x, g);
	gwcopy (x, x2); gtog (g, g2);

/* Test 50 squarings */	

	gwsetnormroutine (0, 1, 0);		/* Enable error checking */
	for (i = 0; i < 50; i++) {
		gwstartnextfft ((i & 3) == 2);	/* Test POSTFFT sometimes */
						/* Test gwsetaddin without */
						/* and with POSTFFT set */
		if ((i == 45 || i == 46) && abs (CARG) == 1)
			gwsetaddin (-31);
		gwsquare (x);
		diff = fabs (gwsuminp (x) - gwsumout (x));
		if (diff > maxdiff) maxdiff = diff;
		squareg (g);
		if ((i == 45 || i == 46) && abs (CARG) == 1) {
			iaddg (-31, g);
			gwsetaddin (0);
		}
		specialmodg (g);
		if (CHECK_OFTEN && (i & 3) != 2) compare (x, g);
	}
	if (square_only) goto done;
	if (MAXDIFF < 1e50)
		sprintf (buf, "Squares complete. MaxErr=%.8g, SumoutDiff=%.8g/%.8g(%d to 1)\n", MAXERR, maxdiff, MAXDIFF, (int) (MAXDIFF / maxdiff));
	else
		sprintf (buf, "Squares complete. MaxErr=%.10g\n", MAXERR);
	OutputBoth (buf);

/* Test mul by const */

	gwsetmulbyconst (3);
	gwsetnormroutine (0, 1, 1);
	gwsquare (x);
	gwsetnormroutine (0, 1, 0);
	diff = fabs (gwsuminp (x) - gwsumout (x));
	if (diff > maxdiff) maxdiff = diff;
	squareg (g); imulg (3, g);
	specialmodg (g);
	if (CHECK_OFTEN) compare (x, g);

/* Test square carefully */

	if (abs (CARG) == 1) gwsetaddin (-42);
	gwsquare_carefully (x);
	diff = fabs (gwsuminp (x) - gwsumout (x));
	if (diff > maxdiff) maxdiff = diff;
	squareg (g);
	if (abs (CARG) == 1) { iaddg (-42, g); gwsetaddin (0); }
	specialmodg (g);
	if (CHECK_OFTEN) compare (x, g);

/* Test gwaddquick, gwsubquick */

	gwadd3quick (x, x2, x3); gtog (g, g3); addg (g2, g3);
	gwsub3quick (x, x2, x4); gtog (g, g4); subg (g2, g4);
	if (CHECK_OFTEN) {
		specialmodg (g3); compare (x3, g3);
		specialmodg (g4); compare (x4, g4);
	}

/* Test gwadd and gwsub */

	gwadd (x, x); gwadd (x, x); gwadd (x, x); imulg (8, g);
	gwsub (x3, x); subg (g3, g);
	gwadd (x4, x); addg (g4, g);
	gwadd3 (x3, x4, x2); gtog (g3, g2); addg (g4, g2);
	gwsub3 (x3, x, x4); gtog (g3, g4); subg (g, g4);
	if (CHECK_OFTEN) {
		specialmodg (g); compare (x, g);
		specialmodg (g2); compare (x2, g2);
		specialmodg (g4); compare (x4, g4);
	}
	gwadd (x2, x); addg (g2, g);
	gwadd (x4, x); addg (g4, g);

/* Test gwaddsub */

	gwaddsub (x, x2);	// compute x+x2 and x-x2
	subg (g, g2);		// x2-x
	addg (g, g);		// x+x
	addg (g2, g);		// x+x2
	negg (g2);		// x-x2
	gwaddsub4 (x, x2, x3, x4);	// compute x+x2 and x-x2
	gtog (g, g3); addg (g2, g3);
	gtog (g, g4); subg (g2, g4);
	if (CHECK_OFTEN) {
		specialmodg (g); compare (x, g);
		specialmodg (g2); compare (x2, g2);
		specialmodg (g3); compare (x3, g3);
		specialmodg (g4); compare (x4, g4);
	}
	gwadd (x2, x); addg (g2, g);
	gwadd (x3, x); addg (g3, g);
	gwadd (x4, x); addg (g4, g);

/* Test gwaddsmall */

	gwaddsmall (x, 111);
	iaddg (111, g); specialmodg (g);
	if (CHECK_OFTEN) compare (x, g);

/* Do some multiplies to make sure that the adds and subtracts above */
/* normalized properly. */

	gwfft (x, x);
	gwfftfftmul (x, x, x);
	diff = fabs (gwsuminp (x) - gwsumout (x));
	if (diff > maxdiff) maxdiff = diff;
	squareg (g); specialmodg (g);
	if (CHECK_OFTEN) compare (x, g);

	gwfft (x, x2); gwcopy (x2, x); gwfftadd3 (x, x2, x4);
	gwfftmul (x4, x3); addg (g3, g3); mulg (g, g3); specialmodg (g3);
	diff = fabs (gwsuminp (x3) - gwsumout (x3));
	if (diff > maxdiff) maxdiff = diff;
	if (CHECK_OFTEN) compare (x3, g3);
	gwfft (x3, x4);
	gwfftfftmul (x4, x2, x);
	diff = fabs (gwsuminp (x) - gwsumout (x));
	if (diff > maxdiff) maxdiff = diff;
	mulg (g3, g); specialmodg (g);
	if (CHECK_OFTEN) compare (x, g);

/* Do the final compare */

done:	if (!CHECK_OFTEN) compare (x, g);
	if (GWERROR) OutputBoth ("GWERROR set during calculations.\n");

/* Print final stats */

	if (maxdiff > MAXDIFF) OutputBoth ("Sumout failed during test.\n");
	if (MAXDIFF < 1e50)
		sprintf (buf, "Test complete. MaxErr=%.8g, SumoutDiff=%.8g/%.8g(%d to 1)\n", MAXERR, maxdiff, MAXDIFF, (int) (MAXDIFF / maxdiff));
	else
		sprintf (buf, "Test complete. MaxErr=%.10g\n", MAXERR);
	OutputBoth (buf);
	OutputBoth ("\n");
}


void test_methodically (void)
{
	int	kbits, cbits;
	double	k;
	unsigned int b;
	int	c;
	unsigned int n;

// Test general mods
// Test rational mods
// Test near fftlen max and test near fft min

	b = 2;
//	for k = 1 bit to 48 bits
//		for c = 1 bit to +/-24 bits
//			for n = ? to ? to test each fftlen

	set_seed ();
	k = gen_k (kbits);
	c = gen_c (cbits);
	gwsetup (k, b, n, c, 0);
// print number as string
// print fftlen
	test_it (1);
	gwdone ();
}

void test_randomly (void)
{
	int	kbits, cbits;
	double	k;
	unsigned int b;
	int	c;
	unsigned int n;
	char	buf[100], fft_desc[100];

	set_seed ();
	b = 2;
	for ( ; ; ) {
		if (MAX_K_BITS == MIN_K_BITS) kbits = 0;
		else kbits = (unsigned int) rand () % (MAX_K_BITS - MIN_K_BITS + 1);
		kbits += MIN_K_BITS;
		if (kbits <= 20)
			cbits = (unsigned int) rand () % MAX_C_BITS_FOR_SMALL_K + 1;
		else
			cbits = (unsigned int) rand () % MAX_C_BITS_FOR_LARGE_K + 1;
		k = gen_k (kbits);
		c = gen_c (cbits);
		n = gen_n ();

		if (ALL_BENCH) {
			test_it_all (k, b, n, c);
		}

		else {
			switch (rand () % 5) {
			case 0:	CPU_L2_CACHE_SIZE = 128; break;
			case 1:	CPU_L2_CACHE_SIZE = 256; break;
			case 2:	CPU_L2_CACHE_SIZE = 512; break;
			case 3:	CPU_L2_CACHE_SIZE = 1024; break;
			case 4:	CPU_L2_CACHE_SIZE = -1; break;
			}
			gwsetup (k, b, n, c, 0);
			sprintf (buf, "Starting QA run on %s, kbits=%ld, cbits=%ld\n",
				 gwmodulo_as_string (), kbits, cbits);
			OutputBoth (buf);
			gwfft_description (fft_desc);
			if (FFTLEN < 140000)
				sprintf (buf,
					 "Using %s, bits_per_word=%.5g/%.5g\n",
					 fft_desc, virtual_bits_per_word (),
					 fft_max_bits_per_word);
			else
				sprintf (buf, "Using %s, bits_per_word=%.5g/%.5g, L2_cache_size=%d\n",
					 fft_desc, virtual_bits_per_word (),
					 fft_max_bits_per_word, CPU_L2_CACHE_SIZE);
			OutputBoth (buf);
			test_it (0);
			gwdone ();
		}

		if (UNIQUE) {
			if (unique_count == 0) return;
		}
	}
}
