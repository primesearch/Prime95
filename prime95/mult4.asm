; Copyright 1995-2000 Just For Fun Software, Inc., all rights reserved
; Author:  George Woltman
; Email: woltman@alum.mit.edu
;
; This routine implements a discrete weighted transform to quickly multiply
; two numbers.
;
; This code handles FFTs that use the very convoluted memory model.
; FFT sizes greater than 64K doubles are supported.
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
INCLUDE fft3.mac
INCLUDE memory.mac

IFDEF PPRO
INCLUDE	lucasp.mac
ENDIF

EXTRNP	pass1_6_levels_fft
EXTRNP	pass1_6_levels_unfft
EXTRNP	pass1_6a_levels_fft
EXTRNP	pass1_6a_levels_unfft
EXTRNP	pass1_6_levels_fftp
EXTRNP	pass1_6_levels_unfftp
EXTRNP	pass1_7_levels_fft
EXTRNP	pass1_7_levels_unfft
EXTRNP	pass1_7a_levels_fft
EXTRNP	pass1_7a_levels_unfft
EXTRNP	pass1_7_levels_fftp
EXTRNP	pass1_7_levels_unfftp
EXTRNP	pass2_8_levels_type_123
EXTRNP	pass2_8_levels_type_123p
EXTRNP	pass2_8_levels_type_4
EXTRNP	pass2_8_levels_type_4p

	very_convoluted_distances

;; All the FFT routines for each FFT length

PROCP	_gw_ffts4
	EXPANDING = 2
	fft	80K
	fft	96K
	fft	112K
	fft	128K
	fft	128K, p
	fft	160K
	fft	192K
	fft	224K
	fft	256K
	fft	256K, p
	fft	320K
	fft	384K
	fft	448K
	fft	512K
	fft	512K, p
	fft	640K
	fft	768K
	fft	896K
	fft	1024K
	fft	1024K, p
	fft	1280K
	fft	1536K
	fft	1792K
	fft	2048K
	fft	2048K, p
	fft	2560K
	fft	3072K
	fft	3584K
	fft	4096K
	fft	4096K, p
ENDPP	_gw_ffts4

_TEXT32	ENDS
END
