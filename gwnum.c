/*----------------------------------------------------------------------
| This file contains the C routines and global variables that are used
| in the multi-precision arithmetic routines.  That is, all routines
| that deal with the gwnum data type.
+---------------------------------------------------------------------*/

/* global variables */

EXTERNC unsigned long PARG=0;	/* The exponent we are testing */
EXTERNC unsigned long FFTLEN=0;	/* The FFT size we are using */
EXTERNC unsigned long PLUS1=0;	/* True if factoring 2^P+1 */
EXTERNC unsigned long GWERROR=0;/* True if an error is detected */
EXTERNC double MAXERR = 0.0;	/* Convolution error in a multiplication */
EXTERNC double MAXDIFF = 0.0;	/* Maximum allowable difference between */
				/* sum of inputs and outputs */
EXTERNC void (*GWPROCPTRS[11])()={NULL}; /* Ptrs to assembly routines */
EXTERNC double *NORM_ARRAY1 = NULL; /* Array #1 of normalization mults */
EXTERNC double *NORM_ARRAY2 = NULL; /* Array #2 of normalization mults */
EXTERNC unsigned long INFP=0;	/* For assembly language arg passing */
EXTERNC unsigned long INFF=0;	/* For assembly language arg passing */
EXTERNC unsigned long INFT=0;	/* For assembly language arg passing */
EXTERNC void *SRCARG = NULL;	/* For assembly language arg passing */
EXTERNC void *SRC2ARG = NULL;	/* For assembly language arg passing */
EXTERNC void *DESTARG = NULL;	/* For assembly language arg passing */
EXTERNC void *DEST2ARG = NULL;	/* For assembly language arg passing */
unsigned long fft_count = 0;	/* Count of forward and inverse FFTs */
void	*gwnum_memory;		/* Allocated memory */

gwnum	*gwnum_alloc = NULL;	/* Array of allocated gwnums */
unsigned int gwnum_alloc_count = 0; /* Count of allocated gwnums */
unsigned int gwnum_alloc_array_size = 0; /* Size of gwnum_alloc array */
gwnum	*gwnum_free = NULL;	/* Array of available gwnums */
unsigned int gwnum_free_count = 0; /* Count of available gwnums */

EXTERNC unsigned long CARRYH=0;	/* For multi-precision asm routines */
EXTERNC unsigned long CARRYL=0;
EXTERNC unsigned long RES=0;

/* Allocate memory and initialize assembly code for arithmetic */
/* modulo 2^N-1 or 2^N+1 */

void gwsetup (
	unsigned long p,	/* Exponent to test */
	unsigned long fftlen,	/* Specific FFT size to use (or zero) */
	int	fft_type)	/* 0 for mod 2^N-1, +1 for mod 2^N+1 */
{
	unsigned long mem_needed;

	PARG = p;
	PLUS1 = fft_type;
	FFTLEN = fftlen;
	mem_needed = gwsetup1 ();
	gwnum_memory = malloc (32 + mem_needed);
	SRCARG = (void *) (((unsigned long) gwnum_memory + 31) & 0xFFFFFFE0);
	gwsetup2 ();
#ifdef MEM_MEASURE
{
int x;
__asm mov x, eax
char buf[80];
sprintf (buf, "%d, mem: %d\n", FFTLEN, x);
OutputBoth(buf);
}
#endif
	MAXERR = 0.0;
	GWERROR = 0;

/* Compute maximum allowable difference for error checking */
/* This error check is disabled for mod 2^N+1 arithmetic */

	if (PLUS1)
		MAXDIFF = 1.0E80;

/* We have observed that the difference seems to vary based on */
/* two times the number of bits per double + log (FFTLEN).  This makes */
/* sense as this is the number of bits in each result word of the FFT. */
/* Call this value in the previous sentence "bits".  The maximum */
/* difference ever observed is 2 ^ (bits - 49.97).  For safety sake, */
/* we will tolerate twice that amount. */

	else {
		double bits;
		bits = (double) p / (double) FFTLEN;
		if (bits < 12.0) bits = 12.0;
		bits = 2.0 * bits + log ((double) FFTLEN) / log (2.0);
		MAXDIFF = pow (2.0, bits - 48.97);
	}

/* Clear counters */

	fft_count = 0;

/* Default size of gwnum_alloc array is 50 */

	gwnum_alloc = NULL;
	gwnum_alloc_count = 0;
	gwnum_alloc_array_size = 50;
	gwnum_free = NULL;
	gwnum_free_count = 0;
}

/* Cleanup any memory allocated for multi-precision math */

void gwdone ()
{
	unsigned int i;
	free (gwnum_memory);
	free (gwnum_free);
	if (gwnum_alloc != NULL) {
		for (i = 0; i < gwnum_alloc_count; i++)
			free (*(char**)((char *) gwnum_alloc[i] - 32));
		free (gwnum_alloc);
	}
}

/* Routine to allocate aligned memory for our big numbers */
/* Memory is allocated on 32-byte boundaries, with an additional */
/* 32 bytes prior to the data for storing useful stuff */

gwnum gwalloc ()
{
	unsigned long size;
	char	*p, *q;

/* Return cached gwnum if possible */

	if (gwnum_free_count)
		return (gwnum_free[--gwnum_free_count]);

/* Allocate arrays if necessary */

	if (gwnum_alloc == NULL) {
		gwnum_free = (gwnum *)
			malloc (gwnum_alloc_array_size * sizeof (gwnum));
		if (gwnum_free == NULL) return (NULL);
		gwnum_alloc = (gwnum *)
			malloc (gwnum_alloc_array_size * sizeof (gwnum));
		if (gwnum_alloc == NULL) return (NULL);
	} else if (gwnum_alloc_count == gwnum_alloc_array_size) {
		gwnum_alloc_array_size += gwnum_alloc_array_size >> 1;
		gwnum_free = (gwnum *)
			realloc (gwnum_free,
				 gwnum_alloc_array_size * sizeof (gwnum));
		if (gwnum_free == NULL) return (NULL);
		gwnum_alloc = (gwnum *)
			realloc (gwnum_alloc,
				 gwnum_alloc_array_size * sizeof (gwnum));
		if (gwnum_alloc == NULL) return (NULL);
	}

/* Use addr function on the last FFT value to compute the size. */

	size = gwnum_size (FFTLEN);
	p = (char *) malloc (32 + 32 + size);
	if (p == NULL) return (NULL);
	q = (char *) (((unsigned long) p + 63) & 0xFFFFFFE0);
	* (char **) (q - 32) = p;
	* (unsigned long *) (q - 8) = size;
	* (double *) (q - 16) = 0.0;
	* (double *) (q - 24) = 0.0;

/* Save pointer for easier cleanup */

	gwnum_alloc[gwnum_alloc_count++] = (gwnum) q;

/* Return the gwnum */

	return ((gwnum) q);
}

/* Free one of our special numbers */

void gwfree (
	gwnum	q)
{
	gwnum_free[gwnum_free_count++] = q;
}

/* Specialized routines that let the giants code share the free */
/* memory pool used by gwnums. */

void gwfree_temporarily (
	gwnum	q)
{
	gwfree (q);
}
void gwrealloc_temporarily (
	gwnum	q)
{
	unsigned long i, j;

	for (i = j = 0; i < gwnum_free_count; i++)
		if (gwnum_free[i] != q) gwnum_free[j++] = gwnum_free[i];
	gwnum_free_count = j;
}

/* Free all of our special numbers */

void gwfreeall ()
{
	unsigned int i;
	if (gwnum_alloc == NULL) return;
	for (i = 0; i < gwnum_alloc_count; i++)
		gwnum_free[i] = gwnum_alloc[i];
	gwnum_free_count = gwnum_alloc_count;
}

/* To optimize use of the L1 cache we scramble the FFT data. */
/* Note:  The Intel L1 data cache is 8KB two-way set associative with */
/* 32 byte cache lines.  Later CPUs have more cache, but we are prepared */
/* for the worst case.  This tiny cache will require us to perform */
/* three "passes" to perform a large FFT. Each pass must minimize */
/* L1 cache line conflicts - that is have no data at the same address */
/* modulo 4096 */

/* 1) We'd like to do as much work as possible in the final pass (called */
/*    pass 2 in a lot of the code).  Since some cache space is required */
/*    for sine/cosine data, we only use half of the L1 cache for FFT data. */
/*    4KB = 512 values = 256 complex values.  Thus, the final pass will */
/*    perform 8 FFT levels.  Also note that it will be advantageous to */
/*    have the real and imaginary values in the same cache line.  Thus, */
/*    the first cache line contains the 0th, 128th, 256th, and 384th FFT */
/*    data values.  Where the 0th and 128th values comprise a single */
/*    complex number as does the 256th and 384th. */

/* 2) To eliminate cache line conflicts in the middle pass (called pass 1 */
/*    in a lot of this code), 32 bytes is wasted after 4KB of FFT data */
/*    If we did not do this every pass 1 value would try to occupy the */
/*    same L1 cache line! */

/* Putting it all together, for FFTLEN=2^16 you get this memory layout:	*/
/*	0	128	256	384		(32 bytes)		*/
/*	1	129	257	385		(32 bytes)		*/
/*		   etc.							*/
/*	127	255	383	511		(32 bytes)		*/
/*		(32 wasted bytes)					*/
/*	512	640	768	896		(32 bytes)		*/
/*	513	641	769	897		(32 bytes)		*/
/*		   etc.							*/
/*	639	767	895	1023		(32 bytes)		*/
/*		(32 wasted bytes)					*/
/*	1024	1152	1280	1408		(32 bytes)		*/
/*		   etc.							*/

/* Well.... I implemented the above only to discover I had dreadful */
/* performance in pass 1.  How can that be?  The problem is that each  */
/* cache line in pass 1 comes from a different 4KB page.  Therefore, */
/* pass 1 accessed 128 different pages.  This is a problem because the */
/* Pentium chip has only 64 TLBs (translation lookaside buffers) to map */
/* logical page addresses into physical addresses.  So we need to shuffle */
/* the data further so that pass 1 data is on fewer pages while */
/* pass 2 data is spread over more pages. */

/* 1st 4KB page		2nd page	...	18th page	*/
/* 0 128 256 384	waste			waste		*/
/* 512 640 768 896	8 136 264 392		waste		*/
/* ...								*/
/* 7680 ...					waste		*/
/* 1 129 257 385				8192 8320 ...	*/
/* 513 ...            	9 137 265 393		8704 ...	*/
/*               ...						*/
/* 7 ...							*/
/* 519 ...							*/
/*	         ...						*/
/* 7687 ...							*/

/* That is, waste 32 bytes after each 512 FFT data values (4KB). */
/* Except after 8192 FFT data values go to the next 4KB page and waste */
/* the first 16*32 bytes.  If you look carefully at the above, you'll see */
/* that in pass 2 the FFT data (values 0 through 511) comes from the first */
/* 16 4KB pages (actually the waste bytes make this 17 4KB pages).  Similarly,
/* the pass 1 data (values 0 up to 65536 stepping by 256) comes from 16 */
/* different 4KB pages  ---  and there are no L1 cache line conflicts!!! */
/* Furthermore, when accessing pages, the pages are an odd number apart */
/* (1 page apart in pass 2, 17 pages apart in pass 1).  This is good in */
/* distributing the pages uniformly among the 4-way set-associative */
/* TLB cache. */

/* How does the above scheme work for the three pass case?  As you might */
/* imagine, more adjustments are necessary.  When doing a 1M FFT we will */
/* work in three passes.  Pass 2 looks at 0 up to 512 step 1, pass 1 */
/* looks at 0 up to 65536 step 256, and pass 0 looks at 0 up to 1048576 */
/* step 32768.  This corresponds to 5 levels in pass 0, 7 in pass 1, 8 in */
/* pass 2.  Notice above that both pass 1 and pass 2 look at the values */
/* 0 and 128 thus they should be on the same 4KB page.  Likewise, values */
/* 0 and 32768 are both used in pass 0 and pass 1 and should be on the */
/* same 4KB page.  After analyzing the various FFT sizes and TLB hit */
/* patterns, I settled on this memory layout: */

/* 1st 4KB page		*/
/* 0 128 256 384	*/
/* 512 640 768 896	*/
/* 16K 16K+128 ...	*/
/* 16K+512...		*/
/* ...			*/
/* 7*16K+512 ...	*/
/* 1 129 257 385	*/
/* 513 ...            	*/
/* ...			*/
/* 7*16K+512+7 ...	*/

/* This is much like the previous layout except that instead of 16 cache */
/* lines that are 512 apart, there are only 2 cache lines that are 512 apart */
/* and 8 that are 16K apart. */

/* To eliminate the cache line conflicts in pass 2, 32 bytes are wasted */
/* every 4KB.  To eliminate cache line conflicts in pass 1 and to keep */
/* the TLB hits uniform, after 16 pages we waste the rest of the 17th */
/* page and the first 8 cache lines of the 18th page.  To eliminate cache */
/* line conflicts in pass 0 and keep the TLB hits uniform, after 16 sets of */
/* 17 pages we move to the next 4KB page and waste the first 16 cache lines */
/* then after wasting 8 sets of 16, we waste another 32 bytes. */

/* It is now getting hard to visualize the FFT, so this program will print */
/* out the cache lines and TLB distributions for FFTS above 64K. */

#ifdef INCLUDED_PROGRAM
#include <stdio.h> 
unsigned long FFTLEN = 0;
/* Copy the addr function here */
void xmain (int incr, int endpt) { 
long	i, x, tlbs[16], lines[128]; 
for (i = 0; i <= 15; i++) tlbs[i] = 0; 
for (i = 0; i <= 127; i++) lines[i] = 0; 
printf ("\n\nTest fftlen: %d, incr: %d, endpt: %d\n", FFTLEN, incr, endpt);
for (i = 0; i < endpt; i += incr) { 
	x = (long) addr((long*)(32*19), i); 
	printf ("i: %d, addr: %d, page: %d, tlb line: %d, cache line: %d\n", 
		i, x, x >> 12, (x >> 12) & 15, (x >> 5) & 127); 
	tlbs[(x >> 12) & 15]++; lines[(x >> 5) & 127]++; 
} 
printf ("\n\nTLBS:"); for (i=0; i<=15; i++) printf (" %d", tlbs[i]); 
printf ("\n\nCache Lines:"); 
for (i = 0; i <= 127; i++) printf (" %d", lines[i]); 
printf ("\n"); 
} 
int main (int argc, char **argv) { 
FFTLEN = 65536 * 2;  xmain (1, 512); xmain (256, 32768); xmain (16384, FFTLEN);
FFTLEN = 65536 * 4;  xmain (1, 512); xmain (256, 65536); xmain (32768, FFTLEN);
FFTLEN = 65536 * 8;  xmain (1, 512); xmain (256, 32768); xmain (16384, FFTLEN);
FFTLEN = 65536 * 16; xmain (1, 512); xmain (256, 65536); xmain (32768, FFTLEN);
FFTLEN = 65536 * 32; xmain (1, 512); xmain (256, 32768); xmain (16384, FFTLEN);
FFTLEN = 65536 * 64; xmain (1, 512); xmain (256, 65536); xmain (32768, FFTLEN);
}
#endif

/* Below is a table of FFT sizes, FFT levels done in each of the three */
/* passes, L1 cache lines used, logical pages touched, and actual pages */
/* touched.  The logical and actual pages touched can be different because */
/* the waste bytes cause "spillage" of data from one 4KB page onto the next. */

/* FFT    FFT levels	L1 cache    Logical pages  Actual pages	*/
/* size   in each pass	lines used  accessed	   accessed	*/
/* ----   ------------	----------  -------------  ------------	*/
/* 4096K  7/7/8		128/128/128 32/16/16	   33/17/18	*/
/* 2048K  7/6/8		128/64/128  16/16/16	   17/17/18	*/
/* 1024K  5/7/8		32/128/128  8/16/16	   9/17/18	*/
/* 512K	  5/6/8		32/64/128   4/16/16	   5/17/18	*/
/* 256K	  3/7/8		8/128/128   2/16/16	   3/17/18	*/
/* 128K	  3/6/8		8/64/128    1/16/16	   2/17/18	*/
/* 64K	  8/8		128/128	    8/16	   9/18		*/
/* 32K	  7/8		64/128	    4/16	   5/18		*/
/* 16K	  6/8		32/128	    2/16	   3/18		*/
/* 8K	  7/6		128/32	    16/1 (flat memory model)	*/
/* 4K	  6/6		64/32	    8/1 (flat memory model)	*/
/* 2K	  5/6		32/32	    4/1 (flat memory model)	*/
/* 1K	  10		256	    2 (flat memory model)	*/
/* 512	  9		128	    1 (flat memory model)	*/
/* 256	  8		64	    1 (flat memory model)	*/

/* NOTE: I once had the brilliant idea of interleaving the sin/cos data */
/* with the FFT data.  That is, the data occupies the even cache */
/* lines and the sin/cos data is in the odd cache lines.  At first */
/* this seems counter productive, as only 4K of FFT data will now fit */
/* in the 8K L1 cache.  However, if you look at how an FFT operates */
/* you'll see loading FFT data, multiply by sin/cos data, store FFT data, */
/* load next block of FFT data, multiply by sin/cos data, store FFT data, */
/* etc.  By storing the sin/cos data in the odd cache lines, loading the */
/* next block of FFT data will toss out the previous block of FFT data */
/* rather than the reusable sin/cos data. */
/* For some reason, however, interleaving resulted in slower performance. */

unsigned long addr_offset (unsigned long fftlen, unsigned long i)
{
	unsigned long addr, i1, i2, i3, i4, i5, i6;

/* Small FFTs use a flat memory model. */

	if (fftlen <= 1024)
		addr = i * 8;

/* Medium-size FFTs use a near flat memory model.  Waste 32 bytes every 512 */
/* data values. */

	else if (fftlen <= 8192)
		addr = i * 8 + (i >> 9) * 32;

/* Large-size FFTs use a convoluted memory model.  Waste bytes as described */
/* above.  Break the element number into parts for scrambling. */
/* From the assembly code, we have these distances between elements */
/* 	dist1 =	512				*/
/*	dist8 = (4096+32)			*/
/*	dist128 = 8				*/
/*	dist8192 = (17*4096+16*32)		*/

	else if (fftlen <= 65536) {
		i1 = i & 7;  i >>= 3;			/* Bottom 3 bits */
		i2 = i & 15; i >>= 4;			/* Next 4 bits */
		i3 = i & 63; i >>= 6;			/* Next 6 bits */
		addr =  (i1 << 9) + i2 * (4096+32) +
			(i3 << 3) + i * (17*4096+16*32);
	}

/* Extra-large-size FFTs use a very convoluted memory model.  Waste bytes */
/* as described above.  Break the element number into parts for scrambling. */
/* From the assembly code, we have these distances between elements */
/* 	dist1 =	512				*/
/*	dist8 = (4096+32)			*/
/*	dist128 = 8				*/
/*	dist1024 = (17*4096+8*32)		*/
/*	dist16K = 64				*/
/*	dist128K = (16*17*4096+4096+16*32)	*/
/*	dist1M = (8*16*17*4096+8*4096+4096+32)	*/

	else {
		i1 = i & 7;  i >>= 3;			/* Bottom 3 bits */
		i2 = i & 15; i >>= 4;			/* Next 4 bits */
		i3 = i & 7;  i >>= 3;			/* Next 3 bits */
		i4 = i & 15; i >>= 4;			/* Next 4 bits */
		i5 = i & 7;  i >>= 3;			/* Next 3 bits */
		i6 = i & 7;  i >>= 3;			/* Next 3 bits */
		addr =  (i1 << 9) + i2 * (4096+32) +
			(i3 << 3) + i4 * (17*4096+8*32) +
			(i5 << 6) + i6 * (16*17*4096+4096+16*32) +
			i * (8*16*17*4096+8*4096+4096+32);
	}

/* Return the offset */

	return (addr);
}

/* Return the address of ith element in the FFT array */

double *addr (gwnum g, unsigned long i)
{
	return ((double *) ((unsigned long) g + addr_offset (FFTLEN, i)));
}

/* Return the size of a gwnum */

unsigned long gwnum_size (unsigned long fftlen)
{
	return (addr_offset (fftlen, fftlen - 1) + sizeof (double));
}

/* Each FFT word is multiplied by a two-to-phi value.  These */
/* routines set and get the FFT value without the two-to-phi */
/* multiplier. */

void get_fft_value (
	gwnum	g,
	unsigned long i,
	long	*retval)
{
	int	n;
	double	ttmp;

/* Compute the multiplier. */

	if (FFTLEN > 8192) {
		ttmp = FFTLEN * NORM_ARRAY1[(i>>7)*2+1];
		if (i & 127) ttmp *= NORM_ARRAY2[(i&127)*2-1];
	} else if (FFTLEN > 128) {
		ttmp = FFTLEN * NORM_ARRAY1[(i>>5)*2+1];
		if (i & 31) ttmp *= NORM_ARRAY2[(i&31)*2-1];
	} else
		ttmp = FFTLEN * NORM_ARRAY2[i*2+1];

/* Use frexp to scale the value so that it is between 0.5 and 1. */
/* frexp generates 0.5 when we really want 1.0 */

	ttmp = frexp (ttmp, &n);
	if (ttmp == 0.5) ttmp = 1.0;

/* Multiply by two-to-minus-phi to generate an integer. */

	ttmp *= * addr (g, i);
	if (ttmp < -0.5)
		*retval = (long) (ttmp - 0.5);
	else
		*retval = (long) (ttmp + 0.5);
}

void set_fft_value (
	gwnum	g,
	unsigned long i,
	long	val)
{
	int	n;
	double	ttp;

/* Compute the multiplier.  Massage it so that it is between 1 and 2 */

	if (FFTLEN > 8192) {
		ttp = NORM_ARRAY1[(i>>7)*2];
		if (i & 127) ttp *= NORM_ARRAY2[(i&127)*2-2];
		ttp = frexp (ttp, &n) * 2.0;
	} else if (FFTLEN > 128) {
		ttp = NORM_ARRAY1[(i>>5)*2];
		if (i & 31) ttp *= NORM_ARRAY2[(i&31)*2-2];
		ttp = frexp (ttp, &n) * 2.0;
	} else
		ttp = NORM_ARRAY2[i*2];

/* Multiply by two-to-minus-phi to generate the proper float. */

	* addr (g, i) = val * ttp;
}

/* Convert a double to a gwnum */

void dbltogw (double d, gwnum g)
{
	unsigned long i;
	double	base1, base2;
	int	bits;

	bits = (int) (PARG / FFTLEN);
	base1 = (double) (1L << bits);
	base2 = (double) (1L << (bits+1));
	for (i = 0; i < FFTLEN; i++) {
		if (d < base1) {
			set_fft_value (g, i, (long) d);
			d = 0.0;
		} else {
			double	base, rem;
			base = is_big_word (i) ? base2 : base1;
			rem = fmod (d, base);
			set_fft_value (g, i, (long) rem);
			d = (d - rem) / base;
		}
	}
	((long *) g)[-1] = 0;	/* Clear needs-normalize counter */
}

/* Routine to handle a multiplication and modulo operation where */
/* the intermediate multiplication result can be more than 32 bits. */

unsigned long mulmod (
	unsigned long a,
	unsigned long b,
	unsigned long c)
{
	SRCARG = (void*) a;
	SRC2ARG = (void*) b;
	DESTARG = (void*) c;
	emulmod ();
	return ((unsigned long) DESTARG);
}

/* Some words in the FFT data contain floor(p/N), some words contain */
/* floor(p/N)+1 bits.  This function returns TRUE in the latter case. */

int is_big_word (
	unsigned long i)
{
	unsigned long b, s, t;

/* Compute the number of big FFT words. */

	b = PARG % FFTLEN;
	if (b == 0) return (FALSE);

/* Compute the number of small FFT words. */

	s = FFTLEN - b;

/* The big words are uniformly distributed.  Use special */
/* arithmetic to avoid overflows in the b * i operation. */

	t = mulmod (b, i, FFTLEN);
	return (t == 0 || t > s);
}

/* Routine map a bit number into an FFT word and bit within that word */

void bitaddr (
	unsigned long bit,
	unsigned long *word,
	unsigned long *bit_in_word)
{
	unsigned long b, c;

/* What word is the bit in? */

	*word = (unsigned long) ((double) bit * FFTLEN / PARG);

/* Compute the number of bits in the word */

	b = PARG / FFTLEN;
	if (is_big_word (*word)) b++;

/* Compute the bit within the word. */

	c = mulmod (bit, FFTLEN, PARG);
	*bit_in_word = c * b / PARG;
}

/* Given an exponent, determine the fft length */

unsigned long map_exponent_to_fftlen (
	unsigned long p,
	int	fft_type)
{
	unsigned long *info;

/* Get pointer to fft info and return the FFT length */

	gwinfo (p, 0, fft_type);
	info = (unsigned long *) INFT;
	return (info[1]);
}

/* Given an fft length, determine the maximum allowable exponent */

unsigned long map_fftlen_to_max_exponent (
	unsigned long fftlen,
	int	fft_type)
{
	unsigned long *info;

/* Get pointer to fft info and return the FFT length */

	gwinfo (0, fftlen, fft_type);
	info = (unsigned long *) INFT;
	return (info[0]);
}

/* Given an fft length, determine how much memory is used for */
/* normalization and sin/cos tables */

unsigned long map_fftlen_to_memused (
	unsigned long fftlen,
	int	fft_type)
{
	unsigned long *info;

/* Get pointer to fft info and return the FFT length */

	gwinfo (0, fftlen, fft_type);
	info = (unsigned long *) INFT;
	return (info[3]);
}

/* Make a guess as to how long a squaring will take. */

double map_fftlen_to_timing (
	unsigned long fftlen,
	int	fft_type,
	int	cpu_type,
	unsigned long cpu_speed)
{
	double	timing;
	unsigned long *info;

/* Get pointer to fft info */

	gwinfo (0, fftlen, fft_type);
	info = (unsigned long *) INFT;

/* Use my PII-400 timings as a guide. */

	timing = ((float *) info)[2];

/* Since the program is about 10% memory bound, the program will not */
/* speed up linearly with increase in chip speed.  Note, no attempt is */
/* made between 66 MHz memory and 100 MHz memory - we're just returning */
/* an educated guess here. */

	timing = 0.10 * timing + 0.90 * timing * 400.0 / cpu_speed;
	if (cpu_type <= 4) timing = timing * REL_486_SPEED;
	if (cpu_type == 5) timing = timing * REL_PENT_SPEED;
	if (cpu_type == 7) timing = timing * REL_K6_SPEED;
	if (cpu_type == 11) timing = timing * REL_K7_SPEED;
	return (timing);
}
