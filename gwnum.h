/*----------------------------------------------------------------------
| This file contains the headers and definitions that are used
| in the multi-precision IBDWT arithmetic routines.  That is, all routines
| that deal with the gwnum data type.
+---------------------------------------------------------------------*/

/* A psuedo declaration for our big numbers.  The actual pointers to */
/* these big numbers are to the data array.  The 32 bytes prior to the */
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

typedef struct {
/*	char	pad[32];	   Used to track unnormalized add/sub */
				/* and original address */
	double	data[512];	/* The big number broken into chunks */
				/* This array is variably sized. */
} *gwnum;

#define MAX_PRIME	79300000L	/* Maximum number of bits */
#define MAX_PRIME_SSE2	77910000L	/* SSE2 bit limit */
#define MAX_FFTLEN	4194304		/* 4096K FFT */

/* global variables */

EXTERNC unsigned long PARG;	/* The exponent we are testing */
EXTERNC unsigned long FFTLEN;	/* The FFT size we are using */
EXTERNC unsigned long NUMBIG;	/* Number of big words in the FFT */
EXTERNC unsigned long NUMLIT;	/* Number of little words in the FFT */
EXTERNC unsigned long PLUS1;	/* True if factoring 2^P+1 */
EXTERNC unsigned long GWERROR;	/* True if an error is detected */
EXTERNC unsigned long FFTZERO[8];/* Number of fft values to NOT zero during */
				/* post-multiply normalization. */
EXTERNC unsigned long COPYZERO[8];/* Ptrs to help in gwcopyzero */
EXTERNC double MAXERR;		/* Convolution error in a multiplication */
EXTERNC double MAXDIFF;		/* Maximum allowable difference between */
				/* sum of inputs and outputs */
EXTERNC double PROTHVALS[13];	/* Values used in proth mod assembly routines*/
EXTERNC void (*GWPROCPTRS[24])(void);/* Ptrs to assembly routines */
EXTERNC unsigned long INFP;	/* For assembly language arg passing */
EXTERNC unsigned long INFF;	/* For assembly language arg passing */
EXTERNC unsigned long INFT;	/* For assembly language arg passing */
EXTERNC void *SRCARG;		/* For assembly language arg passing */
EXTERNC void *SRC2ARG;		/* For assembly language arg passing */
EXTERNC void *DESTARG;		/* For assembly language arg passing */
EXTERNC void *DEST2ARG;		/* For assembly language arg passing */
extern void *GW_BIGBUF;		/* Optional buffer to allocate gwnums in */
extern unsigned long GW_BIGBUF_SIZE;	/* Size of the optional buffer */
extern gwnum *gwnum_alloc;		/* Array of allocated gwnums */
extern unsigned int gwnum_alloc_count;	/* Count of allocated gwnums */
extern unsigned int gwnum_alloc_array_size;/* Size of gwnum_alloc array */
extern gwnum *gwnum_free;		/* Array of available gwnums */
extern unsigned int gwnum_free_count;	/* Count of available gwnums */
extern unsigned long GW_ALIGNMENT;	/* How to align allocated gwnums */

/* Types of FFTs supported */

#define GW_MERSENNE_MOD		0	/* Have the FFT work mod 2^N-1 */
#define GW_FERMAT_MOD		1	/* Have the FFT work mod 2^N+1 */

/* gwnum routines */

void gwsetup (unsigned long, unsigned long, int);
void gwdone (void);
gwnum gwalloc (void);
void gwfree (gwnum);
void gwfreeall (void);
void dbltogw (double, gwnum);
double *addr (gwnum, unsigned long);
unsigned long gwnum_size (unsigned long);
void get_fft_value (gwnum, unsigned long, long *);
void set_fft_value (gwnum, unsigned long, long);
int is_valid_fft_value (gwnum, unsigned long);
int is_big_word (unsigned long);
void bitaddr (unsigned long, unsigned long *, unsigned long *);
#define gw_set_max_allocs(n)	if (gwnum_alloc==NULL) gwnum_alloc_array_size=n
#define gw_test_for_error()		GWERROR
#define gw_test_illegal_sumout()	(GWERROR & 1)
#define gw_test_mismatched_sums()	(GWERROR & 2)
#define gwsuminp(g)			((g)->data[-2])
#define gwsumout(g)			((g)->data[-3])

unsigned long map_exponent_to_fftlen (unsigned long, int);
unsigned long map_fftlen_to_max_exponent (unsigned long, int);
double map_fftlen_to_timing (unsigned long, int, int, double);
unsigned long map_fftlen_to_memused (unsigned long, int);

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
/* gwsetup	Initializes the gw code for a specific Mersenne exponent */
/* gwswap	Quickly swaps two gw numbers */
/* gwcopy(s,d)	Copies gwnum s to d */
/* gwcopyzero(s,d,n) Copies gwnum s to d and zeroes last N words */
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

#define gw_fft()	(*GWPROCPTRS[0])()
#define gw_square()	(*GWPROCPTRS[1])()
#define gw_mul()	(*GWPROCPTRS[2])()
#define gw_mulf()	(*GWPROCPTRS[3])()
#define gw_copy()	(*GWPROCPTRS[4])()
#define gw_add()	(*GWPROCPTRS[5])()
#define gw_addq()	(*GWPROCPTRS[6])()
#define gw_sub()	(*GWPROCPTRS[7])()
#define gw_subq()	(*GWPROCPTRS[8])()
#define gw_addsub()	(*GWPROCPTRS[9])()
#define gw_addsubq()	(*GWPROCPTRS[10])()
#define gw_copyzero()	(*GWPROCPTRS[11])()
#define gw_prothmod()	(*GWPROCPTRS[12])()
#define gw_addf()	(*GWPROCPTRS[13])()
#define gw_subf()	(*GWPROCPTRS[14])()
#define gw_addsubf()	(*GWPROCPTRS[15])()

#define fftinc(x)	(fft_count += x)
#define gwswap(s,d)	{gwnum t; t = s; s = d; d = t;}
#define gwcopy(s,d)	{SRCARG = s; DESTARG = d; gw_copy ();}
#define gwadd3quick(s1,s2,d){SRCARG = s1; SRC2ARG = s2; DESTARG=d; gw_addq ();}
#define gwsub3quick(s1,s2,d){SRCARG = s2; SRC2ARG = s1; DESTARG=d; gw_subq ();}
#define gwaddquick(s,d)	{SRCARG = s; SRC2ARG = d; DESTARG = d; gw_addq ();}
#define gwsubquick(s,d)	{SRCARG = s; SRC2ARG = d; DESTARG = d; gw_subq ();}
#define gwadd(s,d)	{SRCARG = s; SRC2ARG = d; DESTARG = d; gw_add ();}
#define gwsub(s,d)	{SRCARG = s; SRC2ARG = d; DESTARG = d; gw_sub ();}
#define gwadd3(s1,s2,d)	{SRCARG = s1; SRC2ARG = s2; DESTARG = d; gw_add ();}
#define gwsub3(s1,s2,d)	{SRCARG = s2; SRC2ARG = s1; DESTARG = d; gw_sub ();}
#define gwaddsub(a,b)	{SRCARG=a;SRC2ARG=b;DESTARG=a;DEST2ARG=b;gw_addsub();}
#define gwaddsub4(s1,s2,d1,d2){SRCARG=s1;SRC2ARG=s2;DESTARG=d1;DEST2ARG=d2;gw_addsub();}
#define gwaddsubquick(a,b){SRCARG=a;SRC2ARG=b;DESTARG=a;DEST2ARG=b;gw_addsubq();}
#define gwaddsub4quick(s1,s2,d1,d2){SRCARG=s1;SRC2ARG=s2;DESTARG=d1;DEST2ARG=d2;gw_addsubq();}
#define gwfft(s,d)	{SRCARG = s; DESTARG = d; gw_fft (); fftinc(1);}
#define gwsquare(s)	(DESTARG = s, gw_square(), fftinc(2))
#define gwfftmul(s,d)	{SRCARG = s; DESTARG = d; gw_mul (); fftinc(2);}
#define gwfftfftmul(s,s2,d){SRCARG=s;SRC2ARG=s2;DESTARG=d;gw_mulf();fftinc(1);}
#define gwmul(s,d)	{gwfft(s,s); gwfftmul(s,d);}
#define gwsafemul(s,d)	{gwnum qqq;qqq=gwalloc();gwfft(s,qqq); gwfftmul(qqq,d);gwfree(qqq);}
#define gwinfo(p,f,t)	{INFP=p;INFF=f;INFT=t;gwinfo1();}
#define gwsetmulbyconst(s) {SRCARG = (void*)(long)(s); eset_mul_const();}
#define gwsetnormroutine(z,e,c) {NORMRTN=GWPROCPTRS[16+4*((z)!=0)+2*(c)+(e)];if(z)gwsetzero(z);}
#define gwstartnextfft(f) {POSTFFT=f;}
#define gwprothmod(s)	{SRCARG=s;gw_prothmod();}
#define gwtouch(s)	gwcopy(s,s)
#define gwfftadd(s,d)	gwfftadd3(s,d,d)
#define gwfftsub(s,d)	gwfftsub3(d,s,d)
#define gwfftaddsub(a,b) gwfftaddsub4(a,b,a,b)
#define gwfftadd3(s1,s2,d){SRCARG = s1; SRC2ARG = s2; DESTARG=d; gw_addf ();}
#define gwfftsub3(s1,s2,d){SRCARG = s2; SRC2ARG = s1; DESTARG=d; gw_subf ();}
#define gwfftaddsub4(s1,s2,d1,d2){SRCARG=s1;SRC2ARG=s2;DESTARG=d1;DEST2ARG=d2;gw_addsubf();}
void gwcopyzero (gwnum s, gwnum d, unsigned long n);
void gwsetaddin (unsigned long, long);
void gwsetzero (unsigned long);

/* Mod k*2^n+/-1 routines.  Only works if bits-per-word is an integer */

void gwprothsetup (unsigned long k, unsigned long n, int inc);

/* Speed of other processors compared to a Pentium II of same clock speed */

#define REL_486_SPEED	8.4	/* 486 is over 8 times slower than PII */
#define REL_K6_SPEED	3.0	/* K6 is 3 times slower than PII */
#define REL_PENT_SPEED	1.2	/* Pentium is 20% slower than PII */
#define REL_K7_SPEED	0.7	/* Assume K7 is faster than a PII */
#define REL_P4_SPEED	0.7	/* Assume P4 is faster than a PII factoring*/

/* Other low-level math routines the caller can use for multi-precision */
/* arithmetic */

EXTERNC unsigned long CARRYH;	/* For multi-precision asm routines */
EXTERNC unsigned long CARRYL;
EXTERNC unsigned long RES;

EXTERNC void eaddhlp (void);
EXTERNC void esubhlp (void);
EXTERNC void emuladdhlp (void);
EXTERNC void emuladd2hlp (void);
EXTERNC void emulsubhlp (void);
#define addhlp(a)	SRCARG=(void*)a, eaddhlp()
#define subhlp(a)	SRCARG=(void*)a, esubhlp()
#define muladdhlp(a,b)	{SRCARG=(void*)a; SRC2ARG=(void*)b; emuladdhlp();}
#define muladd2hlp(a,b)	{SRCARG=(void*)a; SRC2ARG=(void*)b; emuladd2hlp();}
#define mulsubhlp(a,b)	{SRCARG=(void*)a; SRC2ARG=(void*)b; emulsubhlp();}

/* Specialized routines that let the giants code share the free */
/* memory pool used by gwnums. */

void gwfree_temporarily (gwnum);
void gwrealloc_temporarily (gwnum);

/* Other routines used internally */

unsigned long addr_offset (unsigned long, unsigned long);
EXTERNC void eset_mul_const (void);
