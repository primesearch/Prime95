; Copyright 2001-2003 Just For Fun Software, Inc., all rights reserved
; Author:  George Woltman
; Email: woltman@alum.mit.edu
;
; This routine implements part of a discrete weighted transform to
; quickly multiply two numbers.
;
; This code handles the last 11 levels of two pass FFTs that use the
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

PREFETCHING = 1

;; Routines to do the last 11 levels in a two-pass FFT

PUBLIC	xpass2_11_levels

;; Do the last 11 levels of a two pass FFT

xpass2_11_levels PROC NEAR
start_timer 2
	xpass2_11_levels_real 0
	add	esi, pass1blkdst
end_timer 2
start_timer 3
	mov	ecx, count1		; Number of complex iterations
	mov	edx, pass2_premults	; Address of the group multipliers
p2lp:	xpass2_11_levels_complex 0
	add	esi, pass1blkdst
	sub	ecx, 1
	JNZ_X	p2lp
	mov	esi, _DESTARG		; Restore source pointer
end_timer 3
	ret
xpass2_11_levels ENDP

_TEXT32	ENDS
END
