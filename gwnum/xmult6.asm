; Copyright 2001-2005 Just For Fun Software, Inc., all rights reserved
; Author:  George Woltman
; Email: woltman@alum.mit.edu
;
; This routine implements a discrete weighted transform to quickly multiply
; two numbers.
;
; This code uses Pentium 4's SSE2 instructions for very fast FFTs.
; FFT sizes of 40K and above are supported.
; This code does two passes, 10 levels on the second pass.
;
; You will not stand a chance of understanding any of this code without
; thoroughly familiarizing yourself with fast fourier transforms.  This
; code was adapted from an algorithm described in Richard Crandall's article
; on Discrete Weighted Transforms and Large-Integer Arithmetic.
;

	TITLE   setup

IFNDEF X86_64
	.686
	.XMM
	.MODEL	FLAT
ENDIF

INCLUDE	unravel.mac
INCLUDE extrn.mac
INCLUDE xfft6.mac
INCLUDE	xlucas.mac
INCLUDE xmult.mac
INCLUDE xpass1.mac
INCLUDE xpass1sc.mac

EXTRN	xgw_finish_fft:PROC
EXTRN	xgw_carries:PROC
EXTRN	xgw_finish_mult:PROC

EXTRNP	xpass2_10_levels
EXTRNP	xpass2_10_levels_p

_TEXT SEGMENT

;; Pass 2 does 10 FFT levels 2 sets of data (2 * 2^10 complex values =
;; 2^12 doubles = 32KB).

PREFETCHING = 1
blkdst = (4*(8192+128)+GAP2_10_4)
xpass2_levels = 10

;; All the FFT routines for each FFT length

PROCP	_xmm_gw_ffts6
	EXPANDING = 2
;	xfft	20K
;	xfft	24K
;	xfft	24Kp
;	xfft	28K
;	xfft	32K
;	xfft	32Kp
IFDEF trying_ten_levels_pass2
	xfft	40K
	xfft	48K
	xfft	48Kp
	xfft	56K
	xfft	64K
	xfft	64Kp
	xfft	80K
	xfft	96K
	xfft	96Kp
	xfft	112K
	xfft	128K
	xfft	128Kp
ENDIF
;	xfft	160K
;	xfft	192K
;	xfft	192Kp
;	xfft	224K
;	xfft	256K
;	xfft	256Kp

IFDEF AMD
allfft	xfftclm	160K, 4
	xfftclm	160K, 2
allfft	xfftclm	160K, 1
allfft	xfftclm	192K, 4
	xfftclm	192K, 2
allfft	xfftclm	192K, 1
allfft	xfftclm	192Kp, 4
	xfftclm	192Kp, 2
allfft	xfftclm	192Kp, 1
allfft	xfftclm	224K, 4
	xfftclm	224K, 2
allfft	xfftclm	224K, 1
allfft	xfftclm	256K, 4
	xfftclm	256K, 2
allfft	xfftclm	256K, 1
allfft	xfftclm	256Kp, 4
allfft	xfftclm	256Kp, 2
	xfftclm	256Kp, 1
allfft	xfftclm	320K, 4
allfft	xfftclm	320K, 2
	xfftclm	320K, 1
allfft	xfftclm	384K, 4
allfft	xfftclm	384K, 2
	xfftclm	384K, 1
allfft	xfftclm	384Kp, 4
allfft	xfftclm	384Kp, 2
	xfftclm	384Kp, 1
allfft	xfftclm	448K, 4
allfft	xfftclm	448K, 2
	xfftclm	448K, 1
allfft	xfftclm	512K, 4
allfft	xfftclm	512K, 2
	xfftclm	512K, 1
allfft	xfftclm	512Kp, 4
allfft	xfftclm	512Kp, 2
	xfftclm	512Kp, 1
allfft	xfftclm	640K, 4
allfft	xfftclm	640K, 2
allfft	xfftclm	640K, 1
allfft	xfftclm	640K, 0
allfft	xfftclm	768K, 4
allfft	xfftclm	768K, 2
allfft	xfftclm	768K, 1
allfft	xfftclm	768K, 0
allfft	xfftclm	768Kp, 4
allfft	xfftclm	768Kp, 2
allfft	xfftclm	896K, 4
allfft	xfftclm	896K, 2
allfft	xfftclm	896K, 1
allfft	xfftclm	896K, 0
allfft	xfftclm	1024K, 4
allfft	xfftclm	1024K, 2
allfft	xfftclm	1024K, 1
allfft	xfftclm	1024K, 0
allfft	xfftclm	1024Kp, 2
allfft	xfftclm	1280K, 2
allfft	xfftclm	1280K, 1
allfft	xfftclm	1280K, 0
allfft	xfftclm	1536K, 2
allfft	xfftclm	1536K, 1
allfft	xfftclm	1536K, 0
allfft	xfftclm	1536Kp, 1
allfft	xfftclm	1536Kp, 0
allfft	xfftclm	1792K, 2
allfft	xfftclm	1792K, 1
allfft	xfftclm	1792K, 0
allfft	xfftclm	2048K, 1
allfft	xfftclm	2048K, 0
allfft	xfftclm	2048Kp, 1
allfft	xfftclm	2048Kp, 0
allfft	xfftclm	2560K, 0
allfft	xfftclm	3072K, 0
allfft	xfftclm	3072Kp, 0
allfft	xfftclm	3584K, 0
allfft	xfftclm	4096K, 0
allfft	xfftclm	4096Kp, 0
ELSE
allfft	xfftclm	160K, 4
	xfftclm	160K, 2
allfft	xfftclm	160K, 1
	xfftclm	192K, 4
allfft	xfftclm	192K, 2
	xfftclm	192K, 1
allfft	xfftclm	192Kp, 4
allfft	xfftclm	192Kp, 2
	xfftclm	192Kp, 1
	xfftclm	224K, 4
	xfftclm	224K, 2
allfft	xfftclm	224K, 1
	xfftclm	256K, 4
	xfftclm	256K, 2
allfft	xfftclm	256K, 1
	xfftclm	256Kp, 4
	xfftclm	256Kp, 2
allfft	xfftclm	256Kp, 1
	xfftclm	320K, 4
allfft	xfftclm	320K, 2
	xfftclm	320K, 1
	xfftclm	384K, 4
	xfftclm	384K, 2
allfft	xfftclm	384K, 1
	xfftclm	384Kp, 4
	xfftclm	384Kp, 2
allfft	xfftclm	384Kp, 1
	xfftclm	448K, 4
	xfftclm	448K, 2
allfft	xfftclm	448K, 1
	xfftclm	512K, 4
	xfftclm	512K, 2
allfft	xfftclm	512K, 1
allfft	xfftclm	512K, 0
	xfftclm	512Kp, 4
	xfftclm	512Kp, 2
allfft	xfftclm	512Kp, 1
allfft	xfftclm	512Kp, 0
allfft	xfftclm	640K, 4
allfft	xfftclm	640K, 2
allfft	xfftclm	640K, 1
allfft	xfftclm	640K, 0
allfft	xfftclm	768K, 4
allfft	xfftclm	768K, 2
allfft	xfftclm	768K, 1
	xfftclm	768K, 0
allfft	xfftclm	768Kp, 4
allfft	xfftclm	768Kp, 2
allfft	xfftclm	768Kp, 1
allfft	xfftclm	768Kp, 0
allfft	xfftclm	896K, 4
allfft	xfftclm	896K, 2
allfft	xfftclm	896K, 1
allfft	xfftclm	896K, 0
allfft	xfftclm	1024K, 4
allfft	xfftclm	1024K, 2
allfft	xfftclm	1024K, 1
allfft	xfftclm	1024K, 0
allfft	xfftclm	1024Kp, 2
allfft	xfftclm	1024Kp, 1
allfft	xfftclm	1024Kp, 0
allfft	xfftclm	1280K, 2
allfft	xfftclm	1280K, 1
allfft	xfftclm	1280K, 0
allfft	xfftclm	1536K, 2
allfft	xfftclm	1536K, 1
allfft	xfftclm	1536K, 0
allfft	xfftclm	1536Kp, 1
allfft	xfftclm	1536Kp, 0
allfft	xfftclm	1792K, 2
allfft	xfftclm	1792K, 1
allfft	xfftclm	1792K, 0
allfft	xfftclm	2048K, 1
allfft	xfftclm	2048K, 0
allfft	xfftclm	2048Kp, 1
allfft	xfftclm	2048Kp, 0
allfft	xfftclm	2560K, 0
allfft	xfftclm	3072K, 0
allfft	xfftclm	3072Kp, 0
allfft	xfftclm	3584K, 0
allfft	xfftclm	4096K, 0
allfft	xfftclm	4096Kp, 0
ENDIF

ENDPP	_xmm_gw_ffts6

_TEXT	ENDS
END
