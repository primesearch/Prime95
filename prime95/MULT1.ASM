; Copyright 1995-2000 Just For Fun Software, Inc., all rights reserved
; Author:  George Woltman
; Email: woltman@alum.mit.edu
;
; This routine implements a discrete weighted transform to quickly multiply
; two numbers.
;
; This code handles FFTs that use a simple linear memory model and
; the simplified normalization code.  FFT sizes from 32 doubles to
; 128 doubles are supported.
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
INCLUDE fft1a.mac
INCLUDE memory.mac

IFDEF PPRO
INCLUDE	lucasp.mac
ENDIF

	dist1 =	8

;; All the FFT routines for each FFT length

PROCP	_gw_ffts1
	EXPANDING = 1
	fft	32
	fft	32, p
	fft	40
	fft	48
;	fft	56
	fft	64
	fft	64, p
	fft	80
	fft	96
;	fft	112
	fft	128
	fft	128, p
ENDPP	_gw_ffts1

_TEXT32	ENDS
END
