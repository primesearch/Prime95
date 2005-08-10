/*----------------------------------------------------------------------
| gwrename.h
|
| Handle the difference between the naming conventions in
| C compilers.  We need to do this for global variables that are
| referenced by the assembly routines.
| 
|  Copyright 2004-2005 Just For Fun Software, Inc.
|  All Rights Reserved.
+---------------------------------------------------------------------*/

#ifndef _GWRENAME_H
#define _GWRENAME_H

#define ERRCHK		_ERRCHK
#define MAXERR		_MAXERR
#define GWERROR		_GWERROR

#define KARG		_KARG
#define BARG		_BARG
#define PARG		_PARG
#define CARG		_CARG
#define FFTLEN		_FFTLEN
#define SRCARG		_SRCARG
#define SRC2ARG		_SRC2ARG
#define DESTARG		_DESTARG
#define DEST2ARG	_DEST2ARG
#define NUMARG		_NUMARG
#define NUM2ARG		_NUM2ARG
#define ZERO_PADDED_FFT	_ZERO_PADDED_FFT
#define RATIONAL_FFT	_RATIONAL_FFT
#define SPREAD_CARRY_OVER_4_WORDS _SPREAD_CARRY_OVER_4_WORDS
#define MAXDIFF		_MAXDIFF
#define INFT		_INFT
#define GWVERNUM	_GWVERNUM
#define GWPROCPTRS	_GWPROCPTRS

#define NORMRTN		_NORMRTN
#define BITS_PER_WORD	_BITS_PER_WORD
#define ADDIN_ROW	_ADDIN_ROW
#define ADDIN_OFFSET	_ADDIN_OFFSET
#define ADDIN_VALUE	_ADDIN_VALUE
#define BIGLIT_INCR2	_BIGLIT_INCR2
#define BIGLIT_INCR4	_BIGLIT_INCR4
#define COPYZERO	_COPYZERO
#define ASM_TIMERS	_ASM_TIMERS
#define POSTFFT		_POSTFFT

#define TOP_CARRY_NEEDS_ADJUSTING _TOP_CARRY_NEEDS_ADJUSTING
#define INVERSE_KARG	_INVERSE_KARG
#define KARG_HI		_KARG_HI
#define KARG_LO		_KARG_LO
#define CARRY_ADJUST1	_CARRY_ADJUST1
#define CARRY_ADJUST2	_CARRY_ADJUST2
#define CARRY_ADJUST3	_CARRY_ADJUST3
#define CARRY_ADJUST4	_CARRY_ADJUST4
#define CARRY_ADJUST5	_CARRY_ADJUST5
#define CARRY_ADJUST6	_CARRY_ADJUST6
#define HIGH_WORD1_OFFSET _HIGH_WORD1_OFFSET
#define HIGH_WORD2_OFFSET _HIGH_WORD2_OFFSET
#define HIGH_WORD3_OFFSET _HIGH_WORD3_OFFSET
#define HIGH_SCRATCH1_OFFSET _HIGH_SCRATCH1_OFFSET
#define HIGH_SCRATCH2_OFFSET _HIGH_SCRATCH2_OFFSET
#define HIGH_SCRATCH3_OFFSET _HIGH_SCRATCH3_OFFSET
#define ZPAD_TYPE	_ZPAD_TYPE
#define ZPAD_INVERSE_K6	_ZPAD_INVERSE_K6
#define ZPAD_INVERSE_K5	_ZPAD_INVERSE_K5
#define ZPAD_INVERSE_K4	_ZPAD_INVERSE_K4
#define ZPAD_INVERSE_K3	_ZPAD_INVERSE_K3
#define ZPAD_INVERSE_K2	_ZPAD_INVERSE_K2
#define ZPAD_INVERSE_K1	_ZPAD_INVERSE_K1
#define ZPAD_K6_HI	_ZPAD_K6_HI
#define ZPAD_K6_MID	_ZPAD_K6_MID
#define ZPAD_K6_LO	_ZPAD_K6_LO
#define ZPAD_SHIFT6	_ZPAD_SHIFT6
#define ZPAD_K5_HI	_ZPAD_K5_HI
#define ZPAD_K5_MID	_ZPAD_K5_MID
#define ZPAD_K5_LO	_ZPAD_K5_LO
#define ZPAD_SHIFT5	_ZPAD_SHIFT5
#define ZPAD_K4_HI	_ZPAD_K4_HI
#define ZPAD_K4_MID	_ZPAD_K4_MID
#define ZPAD_K4_LO	_ZPAD_K4_LO
#define ZPAD_SHIFT4	_ZPAD_SHIFT4
#define ZPAD_K3_HI	_ZPAD_K3_HI
#define ZPAD_K3_MID	_ZPAD_K3_MID
#define ZPAD_K3_LO	_ZPAD_K3_LO
#define ZPAD_SHIFT3	_ZPAD_SHIFT3
#define ZPAD_K2_HI	_ZPAD_K2_HI
#define ZPAD_K2_MID	_ZPAD_K2_MID
#define ZPAD_K2_LO	_ZPAD_K2_LO
#define ZPAD_SHIFT2	_ZPAD_SHIFT2
#define ZPAD_K1_HI	_ZPAD_K1_HI
#define ZPAD_K1_LO	_ZPAD_K1_LO
#define ZPAD_SHIFT1	_ZPAD_SHIFT1

#define UU		_UU
#define VV		_VV
#define Ulen		_Ulen
#define Vlen		_Vlen
#define EGCD_A		_EGCD_A
#define EGCD_B		_EGCD_B
#define EGCD_C		_EGCD_C
#define EGCD_D		_EGCD_D
#define EGCD_ODD	_EGCD_ODD
#define RES		_RES
#define CARRYL		_CARRYL
#define CARRYH		_CARRYH

/* Handle the difference in the way the two C compilers name routines */

#define gwinfo1		_gwinfo1
#define gwsetup2	_gwsetup2
#define eset_mul_const	_eset_mul_const
#define timeit		_timeit

#define eaddhlp		_eaddhlp
#define esubhlp		_esubhlp
#define emuladdhlp	_emuladdhlp
#define emulsubhlp	_emulsubhlp
#define emuladd2hlp	_emuladd2hlp
#define egcdhlp		_egcdhlp

#endif
