; Copyright 2001 Just For Fun Software, Inc., all rights reserved
; Author:  George Woltman
; Email: woltman@alum.mit.edu
;
; This routine implements a discrete weighted transform to quickly multiply
; two numbers.
;
; This code uses Pentium 4's SSE2 instructions for very fast FFTs.
; FFT sizes up to 8K are supported in a single pass.
;
; You will not stand a chance of understanding any of this code without
; thoroughly familiarizing yourself with fast fourier transforms.  This
; code was adapted from an algorithm described in Richard Crandall's article
; on Discrete Weighted Transforms and Large-Integer Arithmetic.
;

	TITLE   setup

	.686
	.XMM

_TEXT32 SEGMENT PARA USE32 PUBLIC 'DATA'

	ASSUME  CS: _TEXT32, DS: _TEXT32, SS: _TEXT32, ES: _TEXT32

INCLUDE extrn.mac
INCLUDE	unravel.mac
INCLUDE xfft1.mac
INCLUDE	xlucas.mac
INCLUDE xmult.mac

;; All the FFT routines for each FFT length

_gw_xffts1 PROC NEAR
	EXPANDING = 1
	xfft	32
	xfft	64
	xfft	80
	xfft	96
	xfft	112
	xfft	128
	xfft	160
	xfft	192
	xfft	224
	xfft	256
	xfft	320
	xfft	384
	xfft	448
	xfft	512
	xfft	640
	xfft	768
	xfft	896
	xfft	1024
	xfft	1280
	xfft	1536
	xfft	1792
	xfft	2048
	xfft	2560
	xfft	3072
	xfft	3584
	xfft	4096
	xfft	5120
	xfft	6144
	xfft	7168
	xfft	8192
_gw_xffts1 ENDP

_TEXT32	ENDS
END
