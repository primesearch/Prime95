; Copyright 1995-2000 Just For Fun Software, Inc., all rights reserved
; Author:  George Woltman
; Email: woltman@alum.mit.edu
;
; This routine implements a discrete weighted transform to quickly multiply
; two numbers.
;
; This code handles FFTs that use the convoluted memory model.
; FFT sizes between than 10K and 64K doubles are supported.
; This code does two passes, eight levels on the second pass.
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
INCLUDE fft2b.mac
INCLUDE memory.mac

IFDEF PPRO
INCLUDE	lucasp.mac
ENDIF

EXTRNP	pass2_eight_levels_type_123
EXTRNP	pass2_eight_levels_type_123p
EXTRNP	pass2_eight_levels_type_4
EXTRNP	pass2_eight_levels_type_4p

	convoluted_distances

;; All the FFT routines for each FFT length

PROCP	_gw_ffts3
	EXPANDING = 2
	fft	10K
	fft	12K
	fft	14K
	fft	16K
	fft	16K, p
	fft	20K
	fft	24K
	fft	28K
	fft	32K
	fft	32K, p
	fft	40K
	fft	48K
	fft	56K
	fft	64K
	fft	64K, p
ENDPP	_gw_ffts3

_TEXT32	ENDS
END
