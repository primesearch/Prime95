; Copyright 1995-2000 Just For Fun Software, Inc., all rights reserved
; Author:  George Woltman
; Email: woltman@alum.mit.edu
;
; This routine implements a discrete weighted transform to quickly multiply
; two numbers.
;
; This code handles FFTs that use the almost flat memory model and/or
; the complex normalization code.  FFT sizes from 160 to 8192 doubles
; are supported.  This code does one or two passes, if two passes then
; six levels are done on the second pass.
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
INCLUDE	lucas.mac
INCLUDE pfa.mac
INCLUDE mult.mac
INCLUDE fft1b.mac
INCLUDE fft2a.mac
INCLUDE memory.mac

IFDEF PPRO
INCLUDE	lucasp.mac
ENDIF

EXTRNP	pass2_six_levels_type_123
EXTRNP	pass2_six_levels_type_123p
EXTRNP	pass2_six_levels_type_4
EXTRNP	pass2_six_levels_type_4p

	flat_distances

;; All the FFT routines for each FFT length

PROCP	_gw_ffts2
	EXPANDING = 0
	fft	160
	fft	192
;	fft	224
	fft	256
	fft	256, p
	fft	320
	fft	384
	fft	448
	fft	512
	fft	512, p
	fft	640
	fft	768
	fft	896
	fft	1024
	fft	1024, p
	EXPANDING = 2
	fft	1280
	fft	1536
	fft	1792
	fft	2048
	fft	2048, p
	fft	2560
	fft	3072
	fft	3584
	fft	4096
	fft	4096, p
	fft	5120
	fft	6144
	fft	7168
	fft	8192
	fft	8192, p
ENDPP	_gw_ffts2

_TEXT32	ENDS
END
