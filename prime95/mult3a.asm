; Copyright 1995-2001 Just For Fun Software, Inc., all rights reserved
; Author:  George Woltman
; Email: woltman@alum.mit.edu
;
; This routine implements part of a discrete weighted transform to
; quickly multiply two numbers.
;
; This code handles the last 8 levels of FFTs that use the
; convoluted memory model.
;

	TITLE   setup

	.686
	.XMM

_TEXT32 SEGMENT PARA USE32 PUBLIC 'DATA'

	ASSUME  CS: _TEXT32, DS: _TEXT32, SS: _TEXT32, ES: _TEXT32

INCLUDE extrn.mac
INCLUDE	unravel.mac
INCLUDE	lucas.mac
INCLUDE mult.mac
INCLUDE pass2.mac
INCLUDE memory.mac

IFDEF PPRO
INCLUDE	lucasp.mac
ENDIF

	convoluted_distances

;; Routines to do the last eight levels in a two-pass FFT

	PUBLICP	pass2_eight_levels_type_123
	PUBLICP	pass2_eight_levels_type_123p
	PUBLICP	pass2_eight_levels_type_4
	PUBLICP	pass2_eight_levels_type_4p

PROCP	pass2_procs3

;; Branch to the proper forward FFT code - mod 2^N+1 arithmetic

LABELP	pass2_eight_levels_type_123p
	cmp	ffttype, 2
	JE_X	pass2_eight_levels_type_2p
	JG_X	pass2_eight_levels_type_3p
	JMP_X	pass2_eight_levels_type_1p

;; Branch to the proper forward FFT code

LABELP	pass2_eight_levels_type_123
	cmp	ffttype, 2
	JE_X	pass2_eight_levels_type_2
	JG_X	pass2_eight_levels_type_3
	;; Fall through

;; Do the last eight levels of a forward FFT

pass2_eight_levels_type_1:
	pass2_eight_levels_real 1
pass2_eight_levels_type_1p:
	mov	ecx, count1		;; Number of complex iterations
	mov	edx, pass2_premults	;; Address of the group multipliers
p281lp:	pass2_eight_levels_complex 1, 3
	dec	cl
	jnz	short p281lq
	lea	esi, [esi-64*dist128+dist8192]
	mov	cl, 16
p281lq:	dec	ch
	JNZ_X	p281lp
	pop	esi			;; Pop values and return to C code
	pop	edi
	pop	ebx
	pop	ebp
	ret

;; Do the last eight levels of a squaring (forward FFT, square, inverse FFT)

pass2_eight_levels_type_2:
start_timer 0
	pass2_eight_levels_real 2
pass2_eight_levels_type_2p:
	mov	ecx, count1		;; Number of complex iterations
	mov	edx, pass2_premults	;; Address of the group multipliers
p282lp:	pass2_eight_levels_complex 2, 3
	dec	cl
	jnz	short p282lq
	lea	esi, [esi-64*dist128+dist8192]
	mov	cl, 16
p282lq:	dec	ch
	JNZ_X	p282lp
end_timer 0
	ret

;; Do the last eight levels of a multiply (forward FFT, multiply, inverse FFT)

pass2_eight_levels_type_3:
	pass2_eight_levels_real 3
pass2_eight_levels_type_3p:
	mov	ecx, count1		;; Number of complex iterations
	mov	edx, pass2_premults	;; Address of the group multipliers
p283lp:	pass2_eight_levels_complex 3, 3
	dec	cl
	jnz	short p283lq
	lea	esi, [esi-64*dist128+dist8192]
	mov	cl, 16
p283lq:	dec	ch
	JNZ_X	p283lp
	ret

;; Do the last eight levels of a FFT multiply (multiply, inverse FFT)

LABELP	pass2_eight_levels_type_4
	pass2_eight_levels_real 4
LABELP	pass2_eight_levels_type_4p
	mov	ecx, count1		;; Number of complex iterations
	mov	edx, pass2_premults	;; Address of the group multipliers
p284lp:	pass2_eight_levels_complex 4, 3
	dec	cl
	jnz	short p284lq
	lea	esi, [esi-64*dist128+dist8192]
	mov	cl, 16
p284lq:	dec	ch
	JNZ_X	p284lp
	ret
ENDPP	pass2_procs3

_TEXT32	ENDS
END
