; Copyright 1995-2005 Just For Fun Software, Inc., all rights reserved
; Author:  George Woltman
; Email: woltman@alum.mit.edu
;
; This routine implements the 10 level 2nd pass for FFTs.
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

;; Routines to do the last 10 levels in a two-pass FFT

PUBLICP	pass2_10_levels
PUBLICP	pass2_10_levels_p

;; Distance between two pass 2 data blocks.  Pass 2 does 10 FFT levels
;; 2^10 complex values = 2^11 doubles = 16KB.

blkdst	=	(4*(4096+64)+64)

;; Do the last 10 levels of a two pass FFT

PROCP	pass2_10_levels
	pass2_ten_levels_real
LABELP	pass2_10_levels_p
	mov	ecx, count1		; Number of complex iterations
	mov	edx, pass2_premults	; Address of the group multipliers
p10lp:	pass2_ten_levels_complex
	dec	ecx
	jnz	p10lp
	mov	esi, _DESTARG		; Restore source pointer
	ret
ENDPP	pass2_10_levels



;; Routines to do the last 12 levels in a two-pass FFT

PUBLICP	pass2_12_levels
PUBLICP	pass2_12_levels_p

;; Distance between two pass 2 data blocks.  Pass 2 does 12 FFT levels
;; 2^12 complex values = 2^13 doubles = 64KB.

blkdst	=	(16*(4096+64)+64)

;; Do the last 12 levels of a two pass FFT

PROCP	pass2_12_levels
	pass2_twelve_levels_real
LABELP	pass2_12_levels_p
	mov	ecx, count1		; Number of complex iterations
	mov	edx, pass2_premults	; Address of the group multipliers
p12lp:	pass2_twelve_levels_complex
	dec	ecx
	jnz	p12lp
	mov	esi, _DESTARG		; Restore source pointer
	ret
ENDPP	pass2_12_levels


_TEXT	ENDS
END
