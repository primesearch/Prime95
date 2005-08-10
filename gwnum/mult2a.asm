; Copyright 1995-2005 Just For Fun Software, Inc., all rights reserved
; Author:  George Woltman
; Email: woltman@alum.mit.edu
;
; This routine implements the 8 level 2nd pass for FFTs.
;

	TITLE   setup

	.686
	.XMM
	.MODEL	FLAT

INCLUDE	unravel.mac
INCLUDE extrn.mac
INCLUDE	lucas.mac
INCLUDE mult.mac
INCLUDE pass2.mac
INCLUDE memory.mac

_TEXT SEGMENT

	flat_distances

;; Routines to do the last 8 levels in a two-pass FFT

PUBLICP	pass2_8_levels
PUBLICP	pass2_8_levels_p

;; Distance between two pass 2 data blocks.  Pass 2 does 8 FFT levels
;; 2^8 complex values = 2^9 doubles = 4KB.

blkdst	EQU	(4096+64+64)

;; Do the last 8 levels of a two pass FFT

PROCP	pass2_8_levels
	pass2_eight_levels_real
LABELP	pass2_8_levels_p
	mov	ecx, count1		; Number of complex iterations
	mov	edx, pass2_premults	; Address of the group multipliers
p2lp:	pass2_eight_levels_complex
	dec	ecx
	jnz	p2lp
	mov	esi, _DESTARG		; Restore source pointer
	ret
ENDPP	pass2_8_levels

_TEXT	ENDS
END
