; Copyright 2001-2005 Just For Fun Software, Inc., all rights reserved
; Author:  George Woltman
; Email: woltman@alum.mit.edu
;
; This routine implements part of a discrete weighted transform to
; quickly multiply two numbers.
;
; This code handles the last 11-13 levels of two pass FFTs that use the
; SSE2 instructions.
;

	TITLE   setup

IFNDEF X86_64
	.686
	.XMM
	.MODEL	FLAT
ENDIF

INCLUDE	unravel.mac
INCLUDE extrn.mac
INCLUDE xmult.mac
INCLUDE	xlucas.mac
INCLUDE xpass2.mac

_TEXT SEGMENT

PREFETCHING = 1

;; Routines to do the last 10 levels in a two-pass FFT

PUBLICP	xpass2_10_levels
PUBLICP	xpass2_10_levels_p

PROCP	xpass2_10_levels
	start_timer 2
	xpass2_10_levels_real 0
	add	rsi, pass1blkdst
	end_timer 2
LABELP	xpass2_10_levels_p
	mov	ecx, count1		; Number of complex iterations
	mov	rdx, pass2_premults	; Address of the group multipliers
p10lp:	xpass2_10_levels_complex 0
	add	rsi, pass1blkdst
	sub	ecx, 1
	jnz	p10lp
	mov	rsi, _DESTARG		; Restore source pointer
	ret
ENDPP	xpass2_10_levels


;; Routines to do the last 11 levels in a two-pass FFT

PUBLICP	xpass2_11_levels
PUBLICP	xpass2_11_levels_p

PROCP	xpass2_11_levels
	start_timer 2
	xpass2_11_levels_real 0
	add	rsi, pass1blkdst
	end_timer 2
LABELP	xpass2_11_levels_p
	mov	ecx, count1		; Number of complex iterations
	mov	rdx, pass2_premults	; Address of the group multipliers
p11lp:	xpass2_11_levels_complex 0
	add	rsi, pass1blkdst
	sub	ecx, 1
	jnz	p11lp
	mov	rsi, _DESTARG		; Restore source pointer
	ret
ENDPP	xpass2_11_levels


;; Routines to do the last 12 levels in a two-pass FFT

PUBLICP	xpass2_12_levels
PUBLICP	xpass2_12_levels_p

PROCP	xpass2_12_levels
	start_timer 2
	xpass2_12_levels_real 0
	add	rsi, pass1blkdst
	end_timer 2
LABELP	xpass2_12_levels_p
	mov	ecx, count1		; Number of complex iterations
	mov	rdx, pass2_premults	; Address of the group multipliers
p12lp:	xpass2_12_levels_complex 0
	add	rsi, pass1blkdst
	sub	ecx, 1
	jnz	p12lp
	mov	rsi, _DESTARG		; Restore source pointer
	ret
ENDPP	xpass2_12_levels


;; Routines to do the last 13 levels in a two-pass FFT

PUBLICP	xpass2_13_levels
PUBLICP	xpass2_13_levels_p

PROCP	xpass2_13_levels
	start_timer 2
	xpass2_13_levels_real 0
	add	rsi, pass1blkdst
	end_timer 2
LABELP	xpass2_13_levels_p
	mov	ecx, count1		; Number of complex iterations
	mov	rdx, pass2_premults	; Address of the group multipliers
p13lp:	xpass2_13_levels_complex 0
	add	rsi, pass1blkdst
	sub	ecx, 1
	jnz	p13lp
	mov	rsi, _DESTARG		; Restore source pointer
	ret
ENDPP	xpass2_13_levels


PREFETCHING = 0

;; Routines to do the last 11 levels in a two-pass FFT without prefetching

PUBLICP	xpass2_11_levels_np
PUBLICP	xpass2_11_levels_p_np

PROCP	xpass2_11_levels_np
	start_timer 2
	xpass2_11_levels_real 0
	add	rsi, pass1blkdst
	end_timer 2
LABELP	xpass2_11_levels_p_np
	mov	ecx, count1		; Number of complex iterations
	mov	rdx, pass2_premults	; Address of the group multipliers
p11np:	xpass2_11_levels_complex 0
	add	rsi, pass1blkdst
	sub	ecx, 1
	jnz	p11np
	mov	rsi, _DESTARG		; Restore source pointer
	ret
ENDPP	xpass2_11_levels_np


;; Routines to do the last 12 levels in a two-pass FFT without prefetching

PUBLICP	xpass2_12_levels_np
PUBLICP	xpass2_12_levels_p_np

PROCP	xpass2_12_levels_np
	start_timer 2
	xpass2_12_levels_real 0
	add	rsi, pass1blkdst
	end_timer 2
LABELP	xpass2_12_levels_p_np
	mov	ecx, count1		; Number of complex iterations
	mov	rdx, pass2_premults	; Address of the group multipliers
p12np:	xpass2_12_levels_complex 0
	add	rsi, pass1blkdst
	sub	ecx, 1
	jnz	p12np
	mov	rsi, _DESTARG		; Restore source pointer
	ret
ENDPP	xpass2_12_levels_np


;; Routines to do the last 13 levels in a two-pass FFT without prefetching

PUBLICP	xpass2_13_levels_np
PUBLICP	xpass2_13_levels_p_np

PROCP	xpass2_13_levels_np
	start_timer 2
	xpass2_13_levels_real 0
	add	rsi, pass1blkdst
	end_timer 2
LABELP	xpass2_13_levels_p_np
	mov	ecx, count1		; Number of complex iterations
	mov	rdx, pass2_premults	; Address of the group multipliers
p13np:	xpass2_13_levels_complex 0
	add	rsi, pass1blkdst
	sub	ecx, 1
	jnz	p13np
	mov	rsi, _DESTARG		; Restore source pointer
	ret
ENDPP	xpass2_13_levels_np


_TEXT	ENDS
END
