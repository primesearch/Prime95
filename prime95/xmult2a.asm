; Copyright 2001 Just For Fun Software, Inc., all rights reserved
; Author:  George Woltman
; Email: woltman@alum.mit.edu
;
; This routine implements part of a discrete weighted transform to
; quickly multiply two numbers.
;
; This code handles the last 8 levels of two pass FFTs that use the
; SSE2 instructions.
;

	TITLE   setup

	.686
	.XMM

_TEXT32 SEGMENT PARA USE32 PUBLIC 'DATA'

	ASSUME  CS: _TEXT32, DS: _TEXT32, SS: _TEXT32, ES: _TEXT32

INCLUDE extrn.mac
INCLUDE	unravel.mac
INCLUDE xmult.mac
INCLUDE	xlucas.mac
INCLUDE xpass2.mac

;; Routines to do the last 8 levels in a two-pass FFT

PUBLIC	xpass2_8_levels

;; Distance between two pass 2 data blocks.  Pass 2 does 8 FFT levels
;; 2 sets of data (2 * 2^8 complex values = 2^10 doubles = 8KB).

blkdst	EQU	(8192+128)


;; Do the last 8 levels of a two pass FFT

xpass2_8_levels PROC NEAR
	xpass2_8_levels_real blkdst
	mov	ecx, count1		; Number of complex iterations
	mov	edx, pass2_premults	; Address of the group multipliers
p2lp:	xpass2_8_levels_complex blkdst
	sub	ecx, 1
	JNZ_X	p2lp
	mov	esi, _DESTARG		; Restore source pointer
	ret
xpass2_8_levels ENDP

_TEXT32	ENDS
END
