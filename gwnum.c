/*----------------------------------------------------------------------
| This file contains the C routines and global variables that are used
| in the multi-precision arithmetic routines.  That is, all routines
| that deal with the gwnum data type.
+---------------------------------------------------------------------*/

/* global variables */

EXTERNC unsigned long PARG=0;	/* The exponent we are testing */
EXTERNC unsigned long FFTLEN=0;	/* The FFT size we are using */
EXTERNC unsigned long NUMBIG=0;	/* Number of big words in the FFT */
EXTERNC unsigned long NUMLIT=0;	/* Number of little words in the FFT */
EXTERNC double FFTLEN_INV=0.0;	/* The inverse of FFTLEN */
EXTERNC unsigned long BITS_PER_WORD=0;/* Bits in a little word */
EXTERNC unsigned long PLUS1=0;	/* True if operating modulo 2^P+1 */
EXTERNC unsigned long GWERROR=0;/* True if an error is detected */
EXTERNC double MAXERR = 0.0;	/* Convolution error in a multiplication */
EXTERNC double MAXDIFF = 0.0;	/* Maximum allowable difference between */
				/* sum of inputs and outputs */
EXTERNC double PROTHVALS[13]={0.0};/* Values used in proth mod assembly routines*/
EXTERNC unsigned long FFTZERO[8] = {0};/* Number of fft values to NOT zero */
				/* during a post-multiply normalization. */
EXTERNC unsigned long COPYZERO[8] = {0};/* Ptrs to help in gwcopyzero */
EXTERNC void (*GWPROCPTRS[24])()={NULL}; /* Ptrs to assembly routines */
EXTERNC void (*NORMRTN)() = NULL; /* The post-multiply normalization routine */
EXTERNC unsigned long POSTFFT = 0;/* True if assembly code can start the */
				/* FFT process on the result of a multiply */
EXTERNC unsigned long ADDIN_ROW = 0;/* For adding a constant after multiply */
EXTERNC unsigned long ADDIN_OFFSET = 0;
EXTERNC double ADDIN_VALUE = 0.0;

EXTERNC unsigned long INFP=0;	/* For assembly language arg passing */
EXTERNC unsigned long INFF=0;	/* For assembly language arg passing */
EXTERNC unsigned long INFT=0;	/* For assembly language arg passing */
EXTERNC void *SRCARG = NULL;	/* For assembly language arg passing */
EXTERNC void *SRC2ARG = NULL;	/* For assembly language arg passing */
EXTERNC void *DESTARG = NULL;	/* For assembly language arg passing */
EXTERNC void *DEST2ARG = NULL;	/* For assembly language arg passing */
unsigned long fft_count = 0;	/* Count of forward and inverse FFTs */
void	*gwnum_memory;		/* Allocated memory */
unsigned long GW_ALIGNMENT = 0;	/* How to align allocated gwnums */

gwnum	*gwnum_alloc = NULL;	/* Array of allocated gwnums */
unsigned int gwnum_alloc_count = 0; /* Count of allocated gwnums */
unsigned int gwnum_alloc_array_size = 0; /* Size of gwnum_alloc array */
gwnum	*gwnum_free = NULL;	/* Array of available gwnums */
unsigned int gwnum_free_count = 0; /* Count of available gwnums */

EXTERNC unsigned long CARRYH=0;	/* For multi-precision asm routines */
EXTERNC unsigned long CARRYL=0;
EXTERNC unsigned long RES=0;

/* Assembly helper routines */

EXTERNC void gwsetup2 (void);
EXTERNC void gwinfo1 (void);
EXTERNC void emulmod (void);
EXTERNC void etwo_to_pow (void);
EXTERNC void etwo_to_pow_over_fftlen (void);
EXTERNC void esincos (void);
EXTERNC void esincos3 (void);
EXTERNC void fpu_init (void);

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

/* Find the power of two greater than or equal to N. */

unsigned long pow_two_above_or_equal (
	unsigned long n)
{
static	unsigned long save_n = 0;
static	unsigned long save_result = 0;

	if (n != save_n) {
		save_n = n;
		save_result = 1;
		for (n = n - 1; n; n = n >> 1) save_result = save_result << 1;
	}
	return (save_result);
}


/* Routine to compute 2 ^ (N/FFTLEN) */

double two_to_pow (
	long	n)
{
	double	result;
	SRCARG = (void*) n;
	DESTARG = (void*) &result;
	etwo_to_pow ();
	return (result);
}

/* This routine builds a sin/cos table - used by gwsetup */

double *build_sin_cos_table (
	double	*table,		/* Pointer to the table to fill in */
	unsigned long N,
	int	hermetian_skip,	/* True if some sin/cos values are skipped */
	int	type)		/* 0 = old style - a plain old array */
				/* 1 = SSE2 - data is duplicated */
				/* 2 = SSE2 - data is interleaved */
{
	unsigned long i;

/* Special case the really small sin/cos tables.  If N is between 9 and 16
/* or between 33 and 64, then the assembly code is only doing one FFT level. */
/* In this case, the code just uses the middle sin/cos values of a 2N sized */
/* table.  We could optimize this inefficient memory usage at a later date. */

	if (N <= 8) return (table);
	if (N >= 9 && N <= 16) N = N * 2;
	if (N >= 33 && N <= 64 && type == 1) N = N * 2;

/* Loop to build table. */

	for (i = hermetian_skip ? ((N & 4) ? 4 : 8) : 0; i < N; i += 4) {
		unsigned long shifted_i, shifted_N, flipped_i;
		double	sincos[6];

/* Flip the bits in i.  Our prime-factor-FFT makes this a little complex. */
/* The algorithm below works, but I've long since forgotten why. */

		shifted_i = i; shifted_N = N; flipped_i = 0;
		while ((shifted_N & 1) == 0) {
			flipped_i <<= 1;
			if (shifted_i & 1) flipped_i++;
			shifted_i >>= 1;
			shifted_N >>= 1;
		}
		flipped_i = (flipped_i * shifted_N) + shifted_i;

/* When the FFT is working on real data Hermetian symettry allows us to */
/* eliminate half of the FFT data and consequently half of the sin/cos data */
/* Case 1:  If shifted source is > shifted N/2, then we */
/* do not need these sin/cos values. */
/* Case 2:  If shifted source is zero, loop to find the top */
/* two bits.  Skip the number if the top two bits equal 3. */

		if (hermetian_skip) {
			if (shifted_i > shifted_N / 2) continue;
			if (shifted_i == 0) {
				unsigned long j;
				for (j = i; j > 3; j >>= 1);
				if (j == 3) continue;
			}
		}

/* Compute the 3 sin/cos values */

		SRCARG = (void *) flipped_i;
		SRC2ARG = (void *) N;
		DESTARG = (void *) &sincos;
		esincos3 ();

/* Copy the sin/cos values in the appropriate way */

		if (type == 0) {
			memcpy (table, sincos, sizeof (sincos));
			table += 6;
		} else if (type == 1) {
			table[0] = table[1] = sincos[0];
			table[2] = table[3] = sincos[1];
			table[4] = table[5] = sincos[2];
			table[6] = table[7] = sincos[3];
			table[8] = table[9] = sincos[4];
			table[10] = table[11] = sincos[5];
			table += 12;
		} else if (type == 2) {
			table[0] = sincos[0];
			table[2] = sincos[1];
			table[4] = sincos[2];
			table[6] = sincos[3];
			table[8] = sincos[4];
			table[10] = sincos[5];
			type++;
		} else {
			table[1] = sincos[0];
			table[3] = sincos[1];
			table[5] = sincos[2];
			table[7] = sincos[3];
			table[9] = sincos[4];
			table[11] = sincos[5];
			type--;
			table += 12;
		}
	}
	return (table);
}

/* This routine builds a pass 2 premultiplier table - used by gwsetup */

double *build_premult_table (
	double	*table,		/* Pointer to the table to fill in */
	unsigned long pass2_size)
{
	unsigned long i, N, incr, type;

/* Build a premultiplier table for the second pass assuming either a */
/* eleven level or seven level second pass. */

	N = FFTLEN;
	if (N <= 32768) incr = 256;
	else incr = 2048;

/* Mod 2^N+1 arithmetic uses a half-length FFT */

	if (PLUS1) N >>= 1;

/* Mod 2^N+1 arithmetic starts at first data set, */
/* mod 2^N-1 skips some data sets */

	if (PLUS1) i = 0;
	else i = incr * 4;

/* Loop to build table. */

	type = 0;
	for ( ; i < N; i += incr) {
		unsigned long shifted_i, shifted_N, flipped_i, j, k, l;
		double	sincos[2];

/* Flip the bits in i.  Our prime-factor-FFT makes this a little complex. */
/* The algorithm below works, but I've long since forgotten why. */

		shifted_i = i; shifted_N = N; flipped_i = 0;
		while ((shifted_N & 1) == 0) {
			flipped_i <<= 1;
			if (shifted_i & 1) flipped_i++;
			shifted_i >>= 1;
			shifted_N >>= 1;
		}
		flipped_i = (flipped_i * shifted_N) + shifted_i;

/* When the FFT is working on real data Hermetian symettry allows us to */
/* eliminate half of the FFT data and consequently half of the sin/cos data */
/* Case 1:  If shifted source is > shifted N/2, then we */
/* do not need these sin/cos values. */
/* Case 2:  If shifted source is zero, loop to find the top */
/* two bits.  Skip the number if the top two bits equal 3. */

		if (!PLUS1) {
			if (shifted_i > shifted_N / 2) continue;
			if (shifted_i == 0) {
				unsigned long j;
				for (j = i; j > 3; j >>= 1);	
				if (j == 3) continue;
			}
		}

/* Generate the group multipliers */

		j = 0;
		for (k = 0; k < incr / 4; k += 4) {

/* There are 4 multipliers in a XMM_PMD set */

			for (l = 0; l < 4; l++) {

/* Compute the sin/cos value */

				SRCARG = (void *) ((j + l * incr/4 * flipped_i) % N);
				SRC2ARG = (void *) N;
				DESTARG = (void *) &sincos;
				esincos ();

/* Save the premultiplier value */

				table[l*4+type] = sincos[0];
				table[l*4+2+type] = sincos[1];
			}
			table += 16;

/* Next multiplier */

			j = j + 4 * flipped_i;
		}
	
/* Generate the 4 column multipliers * 4 sin/cos values */

		j = 0;
		for (k = 0; k < 4; k++) {
			for (l = 0; l < 4; l++) {

/* Compute the sin/cos value */

				SRCARG = (void *) ((j + l * N/16) % N);
				SRC2ARG = (void *) N;
				DESTARG = (void *) &sincos;
				esincos ();

/* Save the premultiplier value */

				table[l*4+type] = sincos[0];
				table[l*4+2+type] = sincos[1];
			}
			table += 16;
			j = j + flipped_i;
		}

		if (type == 0) table -= (incr / 4 + 16) * 4;
		type = 1 - type;
 	}

	return (table);
}

/* This routine builds a normalization table - used by SSE2 normalizaion */
/* routines */

double *build_norm_table (
	double	*table,		/* Pointer to the table to fill in */
	int	col)		/* TRUE if building column, not group, table */
{
	unsigned long i, k, half_filled, num_cols;

/* Handle one-pass FFTs first, there are no group multipliers */

	if (FFTLEN <= 8192) {
		if (!col) return (table);

/* Compute the offset of the last word.  In PFA FFTs there is data after */
/* this point, but the cache lines are only half-filled */

		half_filled = addr_offset (FFTLEN, FFTLEN-1) / sizeof (double) + 1;

/* If this is a rational FFT (all FFT words contain the same number of */
/* bits), then create two dummy cache lines containing the value 2/FFTLEN */
/* as the two-to-minus-phi multiplier and 1.0 as the two-to-phi multiplier. */
/* The normalization code requires this. */

		if (NUMLIT == 0) {
			table[0] = table[1] = table[4] = table[5] =
			table[8] = table[9] = table[12] = table[13] = 2.0 / FFTLEN;
			table[2] = table[3] = table[6] = table[7] =
			table[10] = table[11] = table[14] = table[15] = 1.0;
			return (table + 16);
		}

/* Loop to build table */

		for (i = 0; i < FFTLEN; i++) {
			unsigned long j, table_entry;
			double	ttp, ttmp;
			long	n;

/* Call asm routines to compute the two multipliers */

			n = mulmod (i, NUMLIT, FFTLEN);
			ttp = two_to_pow (n);
			SRCARG = (void*) (-n);
			DESTARG = (void*) &ttmp;
			etwo_to_pow_over_fftlen ();

/* Find where this data appears in the FFT array and in the table we are building. */

			j = addr_offset (FFTLEN, i) / sizeof (double);
			table_entry = j >> 1;
			if (j >= half_filled) table_entry -= (j - half_filled) / 4;

/* Now set the entry for the MSW or LSW in an SSE2 pair */

			table[table_entry*4+(j&1)] = ttmp;
			table[table_entry*4+2+(j&1)] = ttp;
		}
		return (table + FFTLEN + FFTLEN);
	}

/* Two pass FFTs are handled here */

	num_cols = (FFTLEN <= 32768) ? 128 : 1024;
	if (col) {

/* If this is a rational FFT (all FFT words contain the same number of */
/* bits), then create two dummy cache lines containing the value 2/FFTLEN */
/* as the two-to-minus-phi multiplier and 1.0 as the two-to-phi multiplier. */
/* The normalization code requires this. */

		if (NUMLIT == 0) {
			table[0] = table[1] = table[4] = table[5] =
			table[8] = table[9] = table[12] = table[13] = 2.0 / FFTLEN;
			table[2] = table[3] = table[6] = table[7] =
			table[10] = table[11] = table[14] = table[15] = 1.0;
			return (table + 16);
		}

/* Loop to build table */

		for (i = 0; i < num_cols; i++) {
			double	ttp, ttmp;
			long	n;

/* Call asm routines to compute the two multipliers */

			n = mulmod (i, NUMLIT, FFTLEN);
			ttp = two_to_pow (n);
			SRCARG = (void*) (-n);
			DESTARG = (void*) &ttmp;
			etwo_to_pow_over_fftlen ();

/* Now set the entry for BOTH the MSW and LSW in an SSE2 pair */

			table[i*4] = ttmp;
			table[i*4+1] = ttmp;
			table[i*4+2] = ttp;
			table[i*4+3] = ttp;
		}
		return (table + num_cols * 4);
	}

/* If this is a rational FFT (all FFT words contain the same number of */
/* bits), then create two dummy cache lines containing the value 2/FFTLEN */
/* as the two-to-minus-phi multiplier and 1.0 as the two-to-phi multiplier. */
/* The normalization code requires this. */

	else {
		unsigned long pfa, h, hlimit, haddin, m, mmult;
		
		if (NUMLIT == 0) {
			table[0] = table[1] = table[4] = table[5] =
			table[8] = table[9] = table[12] = table[13] = 2.0 / FFTLEN;
			table[2] = table[3] = table[6] = table[7] =
			table[10] = table[11] = table[14] = table[15] = 1.0;
			return (table + 16);
		}

/* Determine if this is a PFA 5, 6, 7, or 8 */

		for (pfa = FFTLEN; pfa > 8; pfa >>= 1);

/* Loop to build table */

		hlimit = FFTLEN / 4 / (2*num_cols);
		for (h = 0; h < hlimit; h++) {
			if (pfa == 5) {
				if (h < hlimit / 5) {
					haddin = h * 2 * num_cols;
					mmult = FFTLEN / 20;
				} else {
					haddin = FFTLEN/5 + (h - hlimit/5) * 2 * num_cols;
					mmult = FFTLEN / 5;
				}
			} else if (pfa == 6) {
				if (h < hlimit / 3) {
					haddin = h * 2 * num_cols;
					mmult = FFTLEN / 12;
				} else {
					haddin = FFTLEN/3 + (h - hlimit/3) * 2 * num_cols;
					mmult = FFTLEN / 6;
				}
			} else if (pfa == 7) {
				if (h < hlimit / 7) {
					haddin = h * 2 * num_cols;
					mmult = FFTLEN / 28;
				} else if (h < 3 * hlimit / 7) {
					haddin = FFTLEN/7 + (h - hlimit/7) * 2 * num_cols;
					mmult = FFTLEN / 14;
				} else {
					haddin = 3*FFTLEN/7 + (h - 3*hlimit/7) * 2 * num_cols;
					mmult = FFTLEN / 7;
				}
			} else {
				haddin = h * 2 * num_cols;
				mmult = FFTLEN / 4;
			}
		for (m = 0; m < 4; m++) {
		for (k = 0; k < 2; k++) {
			double	ttp, ttmp;
			long	n;

/* Call asm routines to compute the two multipliers */

			n = mulmod (haddin + m * mmult + k * num_cols, NUMLIT, FFTLEN);
			ttp = two_to_pow (n);
			SRCARG = (void*) (-n);
			DESTARG = (void*) &ttmp;
			etwo_to_pow ();

/* Now set the entry for BOTH the MSW and LSW in an SSE2 pair */

			table[k] = ttmp;
			table[2+k] = ttp;
		}
			table += 4;
		}
		}
		return (table);
	}
}

/* This routine builds a big/little flags table - used by SSE2 normalizaion */
/* routines */

double *build_biglit_table (
	double	*table)		/* Pointer to the table to fill in */
{
	unsigned char *p;
	unsigned long h, i, j, k, m, half_filled, gap, pass1_cache_lines;
	unsigned long pfa, hlimit, haddin, mmult;

/* Handle one pass FFTs differently */

	if (FFTLEN <= 8192) {

/* Compute the offset of the last word.  In PFA FFTs there is data after */
/* this point, but the cache lines are only half-filled */

		half_filled = addr_offset (FFTLEN, FFTLEN-1) / sizeof (double) + 1;

/* Loop to build table */

		p = (unsigned char *) table;
		for (i = 0; i < FFTLEN; i++) {
			unsigned long table_entry;

/* Find where this data appears in the FFT array and in the table we are building. */

			j = addr_offset (FFTLEN, i) / sizeof (double);
			table_entry = j >> 1;
			if (j >= half_filled) table_entry -= (j - half_filled) / 4;

/* Now set the biglit table entry for a LSW in an SSE2 pair */

			if ((j & 1) == 0) {
				p[table_entry] = is_big_word (i) * 16;
			}

/* Otherwise, set the biglit table entry for a MSW in an SSE2 pair */

			else {
				if (is_big_word (i)) p[table_entry] += 32;
			}
		}
		return ((double *) (p + FFTLEN / 2));
	}

/* Determine if this is a PFA 5, 6, 7, or 8 */

	for (pfa = FFTLEN; pfa > 8; pfa >>= 1);

/* Determine the gap between XMM high and low words */

	if (FFTLEN <= 32768) gap = 128;
	else gap = 1024;

	if (FFTLEN <= 1048576) pass1_cache_lines = 8;
	else if (FFTLEN <= 2097152) pass1_cache_lines = 4;
	else pass1_cache_lines = 2;

/* Loop to build table in exactly the same order that it will be */
/* used by the assembly code.  This is especially ugly in the PFA cases */

	p = (unsigned char *) table;
	hlimit = FFTLEN / 4 / (2*gap);
	for (i = 0; i < gap; i += pass1_cache_lines) {
	for (h = 0; h < hlimit; h++) {
		if (pfa == 5) {
			if (h < hlimit / 5) {
				haddin = h * 2 * gap;
				mmult = FFTLEN / 20;
			} else {
				haddin = FFTLEN/5 + (h - hlimit/5) * 2 * gap;
				mmult = FFTLEN / 5;
			}
		} else if (pfa == 6) {
			if (h < hlimit / 3) {
				haddin = h * 2 * gap;
				mmult = FFTLEN / 12;
			} else {
				haddin = FFTLEN/3 + (h - hlimit/3) * 2 * gap;
				mmult = FFTLEN / 6;
			}
		} else if (pfa == 7) {
			if (h < hlimit / 7) {
				haddin = h * 2 * gap;
				mmult = FFTLEN / 28;
			} else if (h < 3 * hlimit / 7) {
				haddin = FFTLEN/7 + (h - hlimit/7) * 2 * gap;
				mmult = FFTLEN / 14;
			} else {
				haddin = 3*FFTLEN/7 + (h - 3*hlimit/7) * 2 * gap;
				mmult = FFTLEN / 7;
			}
		} else {
			haddin = h * 2 * gap;
			mmult = FFTLEN / 4;
		}
	for (j = 0; j < pass1_cache_lines; j++) {
	for (m = 0; m < 4; m++) {
	for (k = 0; k < 2 * gap; k += gap) {
		unsigned long word;

/* Now set the big/little flag for a LSW in an SSE2 pair */
/* Otherwise, set the big/little flag for a MSW in an SSE2 pair */

		word = haddin + i + j + m * mmult + k;
		if (k == 0) *p = is_big_word (word) * 16;
		else if (is_big_word (word)) *p += 32;

/* Set the ttp and ttmp fudge flags for two pass FFTs */
/* The fudge flag is set if the col mult * the grp mult will be greater than 2 */

		if (mulmod (word, NUMLIT, FFTLEN) < mulmod (word & (gap-1), NUMLIT, FFTLEN)) {
			if (k == 0) *p += 64;
			else *p += 128;
		}
	}
	p++;
	}
	}
	}
	}
	return ((double *) p);
}

/* Allocate memory and initialize assembly code for arithmetic */
/* modulo 2^N-1 or 2^N+1 */

void gwsetup (
	unsigned long p,	/* Exponent to test */
	unsigned long fftlen,	/* Specific FFT size to use (or zero) */
	int	fft_type)	/* 0 for mod 2^N-1, +1 for mod 2^N+1 */
{
	unsigned long mem_needed;
	unsigned long *info;

/* The SSE2 code has trouble with real small and real large exponents. */
/* Turn off the SSE2 flag if testing one of these (a simplistic workaround). */

	if (p < 256 || p > 78360000) CPU_FLAGS &= ~CPU_SSE2;

/* The SSE2 code for 2^N+1 testing has not been written yet. */
/* Use the old code for now */

	if (fft_type == 1) CPU_FLAGS &= ~CPU_SSE2;

/* Get pointer to fft info and allocate needed memory */

	fpu_init ();
	gwinfo (p, fftlen, fft_type);
	info = (unsigned long *) INFT;
	mem_needed = info[3];
	gwnum_memory = malloc (mem_needed + 4096);
		
/* Setup some useful global variables */

	PARG = p;
	PLUS1 = fft_type;
	FFTLEN = info[1];
	FFTLEN_INV = 1.0 / (double) FFTLEN;
	BITS_PER_WORD = p / FFTLEN;
	NUMBIG = p % FFTLEN;
	if (NUMBIG == 0) {	/* Asm code is happier with NUMBIG non-zero */
		NUMBIG = FFTLEN;
		BITS_PER_WORD--;
	}
	NUMLIT = FFTLEN - NUMBIG;

/* SSE2 code does much of its initialization in C code. */
/* I sure wish we had done that with the old code!! */
/* Align the allocated memory on a 4KB boundary */

	if (CPU_FLAGS & CPU_SSE2) {
		double *tables;
		unsigned long pass1_size, pass2_size;
		tables = (double *)
			(((unsigned long) gwnum_memory + 4095) & ~4095);
//double *t1=tables;

/* Determine how big the pass 2 size is.  This affects how we build */
/* many of the sin/cos tables. */

		if (FFTLEN <= 8192) pass2_size = 1;
		else if (FFTLEN <= 32768) pass2_size = 256;
		else pass2_size = 2048;

/* Build sin/cos tables used in pass 1.  If FFTLEN is a power of two, */
/* many of the sin/cos tables can be shared. */

		pass1_size = FFTLEN / pass2_size;

		((double **)GWPROCPTRS)[2] = tables;
		tables = build_sin_cos_table (tables, pass1_size, 1, pass2_size == 1 ? 3 : 1);

		if (pass2_size > 1 && pass1_size == pow_two_above_or_equal (pass1_size))
			GWPROCPTRS[3] = GWPROCPTRS[2];
		else {
			((double **)GWPROCPTRS)[3] = tables;
			tables = build_sin_cos_table (tables, pass1_size/4, 1, 1);
		}

		if (pass1_size == pow_two_above_or_equal (pass1_size)) {
			GWPROCPTRS[4] = GWPROCPTRS[3];
			GWPROCPTRS[5] = GWPROCPTRS[3];
			GWPROCPTRS[6] = GWPROCPTRS[3];
		} else {
			((double **)GWPROCPTRS)[4] = tables;
			tables = build_sin_cos_table (tables, pass1_size/16, 1, 1);
			((double **)GWPROCPTRS)[5] = tables;
			tables = build_sin_cos_table (tables, pass1_size/64, 1, 1);
			((double **)GWPROCPTRS)[6] = tables;
			tables = build_sin_cos_table (tables, pass1_size/256, 1, 1);
		}

/* Build sin/cos and premultiplier tables used in pass 2 of two pass FFTs */

		if (pass2_size > 1) {
			((double **)GWPROCPTRS)[0] = tables;
			tables = build_premult_table (tables, pass2_size);
			((double **)GWPROCPTRS)[1] = tables;
			tables = build_sin_cos_table (tables, pass2_size, 0, 1);

			((double **)GWPROCPTRS)[7] = tables;
			tables = build_sin_cos_table (tables, pass2_size * 4, 1, 3);

			((double **)GWPROCPTRS)[8] = tables;
			tables = build_sin_cos_table (tables, pass2_size, 1, 1);

//			if (pass1_size == pow_two_above_or_equal (pass1_size)) {
				GWPROCPTRS[9] = GWPROCPTRS[8];
				GWPROCPTRS[10] = GWPROCPTRS[8];
				GWPROCPTRS[11] = GWPROCPTRS[8];
//			} else {
//				((double **)GWPROCPTRS)[9] = tables;
//				tables = build_sin_cos_table (tables, pass2_size/4, 1, 1);
//				((double **)GWPROCPTRS)[10] = tables;
//				tables = build_sin_cos_table (tables, pass2_size/16, 1, 1);
//				((double **)GWPROCPTRS)[11] = tables;
//				tables = build_sin_cos_table (tables, pass2_size/64, 1, 1);
//			}
		}

/* Build the normalization tables.  The first table is the group normalization */
/* multipliers.  The second table is the column normalization multipliers. */

		((double **)GWPROCPTRS)[12] = tables;
		tables = build_norm_table (tables, 0);
		((double **)GWPROCPTRS)[13] = tables;
		tables = build_norm_table (tables, 1);

/* Build the table of big vs. little flags. */

		((double **)GWPROCPTRS)[14] = tables;
		tables = build_biglit_table (tables);

/* Allocate a table for carries.  Init with XMM_BIGVAL. */

		if (pass2_size > 1) {
			int	i, carry_table_size;
			double	xmm_bigval;
			((double **)GWPROCPTRS)[15] = tables;
			carry_table_size = FFTLEN / (pass2_size / 2);
			xmm_bigval = 3.0 * 131072.0 * 131072.0 * 131072.0;
			for (i = 0; i < carry_table_size; i++) *tables++ = xmm_bigval;
			tables += carry_table_size;
		}
//{
//char buf[80];
//sprintf (buf, "%d, mem: %d\n", FFTLEN, (int) tables - (int) t1);
//OutputBoth(buf);
//}
	}

/* Now call assembly routine to finish off the initialization */
/* Align the allocated memory on a 4KB boundary */

	SRCARG = (void *) (((unsigned long) gwnum_memory + 4095) & ~4095);
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

/* Point to default normalization routines, no proth mod by default */

	gwsetnormroutine (0, 0, 0);
	gwstartnextfft (FALSE);
	gwsetaddin (0, 0);
	* (double **) &PROTHVALS[12] = NULL;

/* Clear globals */

	MAXERR = 0.0;
	GWERROR = 0;
	FFTZERO[0] = 0;
	COPYZERO[0] = 0;

/* Compute maximum allowable difference for error checking */
/* This error check is disabled for mod 2^N+1 arithmetic */

	if (PLUS1)
		MAXDIFF = 1.0E80;

/* We have observed that the difference seems to vary based on the size */
/* the FFT result word.  This is two times the number of bits per double. */
/* Subtract 1 from bits per double because one bit is the sign bit. */
/* Add in a percentage of the log(FFTLEN) to account for carries. */
/* We use a different threshold for SSE2 which uses 64-bit instead of */
/* 80-bit doubles during the FFT */

	else {
		double bits_per_double, total_bits, loglen;
		bits_per_double = (double) p / (double) FFTLEN - 1.0;
		if (bits_per_double < 15.0) bits_per_double = 15.0;
		loglen = log ((double) FFTLEN) / log (2.0);
		loglen *= 0.69;
		total_bits = bits_per_double * 2.0 + loglen * 2.0;
		MAXDIFF = pow (2.0, total_bits -
				((CPU_FLAGS & CPU_SSE2) ? 47.08 : 50.65));
	}

/* Clear counters */

	fft_count = 0;

/* Default size of gwnum_alloc array is 50 */

	gwnum_alloc = NULL;
	gwnum_alloc_count = 0;
	gwnum_alloc_array_size = 50;
	gwnum_free = NULL;
	gwnum_free_count = 0;

/* Compute alignment for allocated data */

	if (CPU_TYPE <= 10)
		GW_ALIGNMENT = 32;		/* P3 and earlier */
	else if (CPU_TYPE == 11)
		GW_ALIGNMENT = 64;		/* Athlon */
	else if (CPU_TYPE >= 12) {		/* P4 and later */
		if (FFTLEN <= 8192)		/* One pass */
			GW_ALIGNMENT = 128;	/* Cache line alignment */
		else				/* Two passes */
			GW_ALIGNMENT = 4096;	/* Page alignment */
	}
}

/* Cleanup any memory allocated for multi-precision math */

void gwdone (void)
{
	unsigned int i;
	free (* (double **) &PROTHVALS[12]);
	free (gwnum_memory);
	free (gwnum_free);
	if (gwnum_alloc != NULL) {
		for (i = 0; i < gwnum_alloc_count; i++)
			free (*(char**)((char *) gwnum_alloc[i] - 32));
		free (gwnum_alloc);
	}
	FFTLEN = 0;

/* The SSE2 code has trouble with real small exponents, real large */
/* exponents and 2^N+1 modulo.  Undo the workaround used in gwsetup. */

	if (PARG < 256 || PARG > 78360000 || PLUS1) setCpuFlags ();
}

/* Routine to allocate aligned memory for our big numbers */
/* Memory is allocated on 128-byte boundaries, with an additional */
/* 32 bytes prior to the data for storing useful stuff */

gwnum gwalloc (void)
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
/* Allocate 32 extra bytes for header information and allocate */
/* extra bytes to assure the data is aligned on a cache line */

	size = gwnum_size (FFTLEN);
	p = (char *) malloc (size + 32 + GW_ALIGNMENT);
	if (p == NULL) return (NULL);
	q = (char *) (((unsigned long) p + 32 + GW_ALIGNMENT - 1) & ~(GW_ALIGNMENT - 1));

/* Initialize the header */

	* (char **) (q - 32) = p;		/* Ptr to free */
	* (unsigned long *) (q - 8) = size;	/* Size in bytes */
	* (unsigned long *) (q - 4) = 0;	/* Unnormalized adds count */
	* (unsigned long *) (q - 28) = 0;	/* Has-been-pre-ffted flag */
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

void gwfreeall (void)
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

/* P4 uses a different memory layout - more suitable to SSE2 */

	if (CPU_FLAGS & CPU_SSE2) {

/* Small FFTs use one pass, not very convoluted.  This the example for	*/
/* a length 2048 FFT:							*/
/*	0	512	1	513	1024	1536	1025	1537	*/
/*	2	...							*/
/*	...								*/
/*	510								*/
/* and PFA-style FFTs are a little tricker.  For example, a 1536 FFT:	*/
/*	0	512	1	513	1024	1280	1025	1281	*/
/*	...								*/
/*	254	766	255	767	1276	1534	1277	1535	*/
/*	256	768	257	769	zero	zero	zero	zero	*/
/*	...								*/
/*	510								*/

		if (fftlen <= 8192) {
			unsigned long pow2len, cachelines, pfafudge;
			pow2len = pow_two_above_or_equal (fftlen);
			pfafudge = (fftlen - pow2len/2) / 2;
			if (i >= pow2len/2 + pfafudge)
				i += (pow2len/4 - pfafudge);
			cachelines = pow2len >> 3;
			i1 = i & 1; i >>= 1;
			i2 = i & (cachelines - 1); i /= cachelines;
			i3 = i & 1; i >>= 1;
			addr = (((((i2 << 1) + i) << 1) + i1) << 1) + i3;
			addr = addr * sizeof (double);
		}

/* Larger FFTs use two passes.  This the example for a length 64K FFT:	*/
/*	0	1K	16K	17K	32K	33K	48K	49K	*/
/*	1	...							*/
/*	...								*/
/*	1023	...							*/
/*	2K	...							*/
/*	...								*/

		else if (fftlen <= 32768) {
			unsigned long sets, pfa, temp;
			sets = FFTLEN >> 10;
			i1 = i & 127; i >>= 7;
			i2 = i & 1; i >>= 1;
			i3 = 0;
			for (pfa = FFTLEN; pfa > 8; pfa >>= 1);
			if (pfa & 1) {		/* pfa is 5 or 7 */
				temp = sets / pfa;
				if (i < temp * 4) {
					i3 = i % temp; i /= temp; pfa = 0;
				} else {
					i3 = temp; i -= temp * 4; sets -= temp;
				}
			}
			if (pfa & 2) {		/* pfa was 6 or 7 */
				temp = sets / 3;
				if (i < temp * 4) {
					i3 += i % temp; i /= temp; pfa = 0;
				} else {
					i3 += temp; i -= temp * 4; sets -= temp;
				}
			}
			if (pfa && sets) {
				i3 += i & (sets - 1); i /= sets;
			}
			addr = (((((i3 * 130) + i1) << 2) + i) << 1) + i2;
			addr = addr * sizeof (double);
		} else {
			unsigned long sets, pfa, temp;
			sets = FFTLEN >> 13;
			i1 = i & 1023; i >>= 10;
			i2 = i & 1; i >>= 1;
			i3 = 0;
			for (pfa = FFTLEN; pfa > 8; pfa >>= 1);
			if (pfa & 1) {		/* pfa is 5 or 7 */
				temp = sets / pfa;
				if (i < temp * 4) {
					i3 = i % temp; i /= temp; pfa = 0;
				} else {
					i3 = temp; i -= temp * 4; sets -= temp;
				}
			}
			if (pfa & 2) {		/* pfa was 6 or 7 */
				temp = sets / 3;
				if (i < temp * 4) {
					i3 += i % temp; i /= temp; pfa = 0;
				} else {
					i3 += temp; i -= temp * 4; sets -= temp;
				}
			}
			if (pfa && sets) {
				i3 += i & (sets - 1); i /= sets;
			}
			addr = (((((i3 * 1090) + i1) << 2) + i) << 1) + i2;
			addr = addr * sizeof (double);
		}
	}

/* Small FFTs use a flat memory model. */

	else if (fftlen <= 1024)
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

/* Return the size of a gwnum.  Note the one pass SSE2 normalization code */
/* requires allocating the same space as the power-of-two fftlen */

unsigned long gwnum_size (unsigned long fftlen)
{
	if ((CPU_FLAGS & CPU_SSE2) && fftlen <= 8192)
		fftlen = pow_two_above_or_equal (fftlen);
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
	double	ttmp;

/* Handle the rational FFT case quickly */

	if (NUMLIT == 0) {
		*retval = (long) * addr (g, i);
		return;
	}

/* Compute the multiplier */

	ttmp = two_to_pow (- (long) mulmod (i, NUMLIT, FFTLEN));

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
	double	ttp;

/* Handle the rational FFT case quickly */

	if (NUMLIT == 0 || val == 0.0) {
		* addr (g, i) = val;
		return;
	}

/* Compute the multiplier */

	ttp = two_to_pow (mulmod (i, NUMLIT, FFTLEN));

/* Multiply by two-to-phi to generate the proper double. */

	* addr (g, i) = val * ttp;
}

/* Convert a double to a gwnum */

void dbltogw (double d, gwnum g)
{
	unsigned long i;
	double	base1, base2;

	base1 = (double) (1L << BITS_PER_WORD);
	base2 = (double) (1L << (BITS_PER_WORD+1));
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
	((long *) g)[-7] = 0;	/* Clear has been FFTed flag */
}

/* Some words in the FFT data contain floor(p/N), some words contain */
/* floor(p/N)+1 bits.  This function returns TRUE in the latter case. */

int is_big_word (
	unsigned long i)
{
	unsigned long t;

/* If all words are the same size, then all words are big words. */

	if (NUMLIT == 0) return (TRUE);

/* The big words are uniformly distributed.  Use special */
/* arithmetic to avoid overflows in the b * i operation. */

	t = mulmod (NUMBIG, i, FFTLEN);
	return (t == 0 || t > NUMLIT);
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

	b = BITS_PER_WORD;
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

/* Use my PII-400 or P4-1400 timings as a guide. */

	timing = ((float *) info)[2];

/* Since the program is about 10% memory bound, the program will not */
/* speed up linearly with increase in chip speed.  Note, no attempt is */
/* made to differentiate between 66 MHz memory and 100 MHz memory - we're */
/* just returning an educated guess here. */

	if (CPU_FLAGS & CPU_SSE2) {
		timing = 0.10 * timing + 0.90 * timing * 1400.0 / cpu_speed;
	} else {
		timing = 0.10 * timing + 0.90 * timing * 400.0 / cpu_speed;
		if (cpu_type <= 4) timing = timing * REL_486_SPEED;
		if (cpu_type == 5) timing = timing * REL_PENT_SPEED;
		if (cpu_type == 7) timing = timing * REL_K6_SPEED;
		if (cpu_type == 11) timing = timing * REL_K7_SPEED;
	}
	return (timing);
}


/* Internal routine to help gwsetzero and gwcopyzero */

void calc8ptrs (
	unsigned long n,
	unsigned long *ptrs)
{
	unsigned long i, j, k;

/* This is a grossly inefficient way to do this.  However, it should */
/* be called rarely. */

	for (i = 0; i < 8; i++) ptrs[i] = 0;
	for (i = 0; i < n; i++) {
		j = addr_offset (FFTLEN, i);
		k = (j & 63) >> 3;
		if (j >= ptrs[k]) ptrs[k] = j - (k << 3) + 64;
	}
}


/* Routine that helps assembly code with zeroing the upper FFT words of */
/* a multiply */

void gwsetzero (
	unsigned long zerocnt)
{
static	unsigned long saved_n = 0;

/* Old FFT code required only a count of words to NOT zero */

	if (! (CPU_FLAGS & CPU_SSE2)) {
		FFTZERO[0] = FFTLEN - zerocnt;
		return;
	}

/* The new P4 code requires 8 pointers to tell us where each column */
/* transitions from not-zeroing to zeroing. */

	if (FFTZERO[0] == 0 || zerocnt != saved_n) {
		saved_n = zerocnt;
		calc8ptrs (FFTLEN - zerocnt, (unsigned long *) FFTZERO);
	}
}


/* Routine that sets up and calls assembly code to copy a gwnum from */
/* source to dest while zeroing some lower FFT words */

void gwcopyzero (
	gwnum	s,
	gwnum	d,
	unsigned long n)
{
static	unsigned long saved_n = 0;

	SRCARG = s;
	DESTARG = d;
	SRC2ARG = (void*)(n);
	if ((CPU_FLAGS & CPU_SSE2) && (COPYZERO[0] == 0 || n != saved_n)) {
		saved_n = n;
		calc8ptrs (n, (unsigned long *) COPYZERO);
	}
	gw_copyzero ();
}


/* Routine that tells the assembly code to add a small value to the */
/* results of each multiply */

void gwsetaddin (
	unsigned long word,
	long	val)
{

/* Remember the offset to the FFT data value */

	ADDIN_OFFSET = addr_offset (FFTLEN, word);

/* If this is a two-pass SSE2 FFT, then we need to tell the assembly code */
/* the affected "row", that is which set of pass 1 data the add-in will */
/* take place */

	if (CPU_FLAGS & CPU_SSE2) {
		unsigned long row, pass1_cache_lines;
		if (FFTLEN <= 1048576) pass1_cache_lines = 8;
		else if (FFTLEN <= 2097152) pass1_cache_lines = 4;
		else pass1_cache_lines = 2;
		if (FFTLEN <= 8192) {
			row = ADDIN_OFFSET & 31;
			if (row == 8) ADDIN_OFFSET += 8;
			if (row == 16) ADDIN_OFFSET -= 8;
		} else if (FFTLEN <= 32768) {
			row = (word & 127) / pass1_cache_lines;
			ADDIN_ROW = 128 / pass1_cache_lines - row;
			ADDIN_OFFSET -= row * pass1_cache_lines * 64;
		} else {
			row = (word & 1023) / pass1_cache_lines;
			ADDIN_ROW = 1024 / pass1_cache_lines - row;
			ADDIN_OFFSET -= row * pass1_cache_lines * 64;
		}
	}

/* Set the addin value - multiply it by two-to-phi and FFTLEN/2 */

	ADDIN_VALUE = (double) val * FFTLEN * 0.5 *
			two_to_pow (mulmod (word, NUMLIT, FFTLEN));
}

/* Mod k*2^n+/-1 routines.  Only works if bits-per-word is an integer */

void gwprothsetup (
	unsigned long k,
	unsigned long n,
	int	inc)			/* Plus or minus one */
{
	unsigned long rem, temp, bits;

	bits = BITS_PER_WORD + 1;	/* All words are big words */

	PROTHVALS[0] = (double) k;
	PROTHVALS[1] = 1.0 / (double) k;

	rem = (PARG - n) % bits;
	PROTHVALS[2] = (double) (1L << rem) * -inc;

	*((unsigned long *) &PROTHVALS[3]) = addr_offset (FFTLEN, FFTLEN - 1);
	*((unsigned long *) &PROTHVALS[4]) = (PARG - n) / bits;
	*((unsigned long *) &PROTHVALS[5]) = ((PARG - n) / bits + 9) / 8;

	temp = (n + bits - 1) / bits;
	*((unsigned long *) &PROTHVALS[6]) = addr_offset (FFTLEN, temp + 4);
	*((unsigned long *) &PROTHVALS[7]) = addr_offset (FFTLEN, temp + 3);
	*((unsigned long *) &PROTHVALS[8]) = addr_offset (FFTLEN, temp + 2);
	*((unsigned long *) &PROTHVALS[9]) = addr_offset (FFTLEN, temp + 1);
	*((unsigned long *) &PROTHVALS[10]) = addr_offset (FFTLEN, temp + 0);

	temp = (long) malloc (((PARG - n) / bits + 24) * sizeof (double) + 8);
	*((double **) &PROTHVALS[12]) = (double *) temp;
	*((double **) &PROTHVALS[11]) = (double *) ((temp + 7) & 0xFFFFFFF8);
}
