; Copyright 1995-1999 Just For Fun Software, Inc., all rights reserved
; Author:  George Woltman
; Email: woltman@alum.mit.edu
;
; This code handles the last 8 levels of three pass FFTs that use the very
; convoluted memory model.
;

	TITLE   setup

	.386

_TEXT32 SEGMENT PARA USE32 PUBLIC 'DATA'

	ASSUME  CS: _TEXT32, DS: _TEXT32, SS: _TEXT32, ES: _TEXT32

INCLUDE extrn.mac
INCLUDE	unravel.mac
INCLUDE	lucas1.mac
INCLUDE mult.mac
INCLUDE pass1.mac
INCLUDE pass2.mac
INCLUDE memory.mac

IFDEF PPRO
INCLUDE	lucas1p.mac
ENDIF

	very_convoluted_distances

	PUBLICP	pass2_8_levels_type_123
	PUBLICP	pass2_8_levels_type_123p
	PUBLICP	pass2_8_levels_type_4
	PUBLICP	pass2_8_levels_type_4p

pass2_procs PROC NEAR

;; Branch to the proper forward FFT code - mod 2^N+1 arithmetic

LABELP	pass2_8_levels_type_123p
	cmp	ffttype, 2
	JE_X	pass2_8_levels_type_2p
	JG_X	pass2_8_levels_type_3p
	JMP_X	pass2_8_levels_type_1p

;; Branch to the proper forward FFT code

LABELP	pass2_8_levels_type_123
	cmp	ffttype, 2
	JE_X	pass2_8_levels_type_2
	JG_X	pass2_8_levels_type_3
	;; Fall through

;; Do the last eight levels of a forward FFT

pass2_8_levels_type_1:
	pass2_eight_levels_real 1
pass2_8_levels_type_1p:
	mov	ecx, count1		;; Number of complex iterations
	mov	edx, pass2_premults	;; Address of the group multipliers
p281lp:	pass2_eight_levels_complex 1
	add	cl, 256/2
	JNC_X	p281lp
	lea	esi, [esi-8*dist128+dist1K]
	dec	cl
	JNZ_X	p281lp
	lea	esi, [esi-16*dist1K+dist16K]
	sub	ecx, 65536-16
	js	short p281dn
	add	ch, 256/8
	JNC_X	p281lp
	lea	esi, [esi-8*dist16K+dist128K]
	dec	ch
	JNZ_X	p281lp
	mov	ch, 8
	lea	esi, [esi-8*dist128K+dist1M]
	JMP_X	p281lp
p281dn:	pop	esi			;; Pop values and return to C code
	pop	edi
	pop	ebx
	pop	ebp
	ret

;; Do the last eight levels of a squaring (forward FFT, square, inverse FFT)

pass2_8_levels_type_2:
	pass2_eight_levels_real 2
pass2_8_levels_type_2p:
	mov	ecx, count1		;; Number of complex iterations
	mov	edx, pass2_premults	;; Address of the group multipliers
p282lp:	pass2_eight_levels_complex 2
	add	cl, 256/2
	JNC_X	p282lp
	lea	esi, [esi-8*dist128+dist1K]
	dec	cl
	JNZ_X	p282lp
	lea	esi, [esi-16*dist1K+dist16K]
	sub	ecx, 65536-16
	js	short p282dn
	add	ch, 256/8
	JNC_X	p282lp
	lea	esi, [esi-8*dist16K+dist128K]
	dec	ch
	JNZ_X	p282lp
	mov	ch, 8
	lea	esi, [esi-8*dist128K+dist1M]
	JMP_X	p282lp
p282dn:	ret

;; Do the last eight levels of a multiply (forward FFT, multiply, inverse FFT)

pass2_8_levels_type_3:
	pass2_eight_levels_real 3
pass2_8_levels_type_3p:
	mov	ecx, count1		;; Number of complex iterations
	mov	edx, pass2_premults	;; Address of the group multipliers
p283lp:	pass2_eight_levels_complex 3
	add	cl, 256/2
	JNC_X	p283lp
	lea	esi, [esi-8*dist128+dist1K]
	dec	cl
	JNZ_X	p283lp
	lea	esi, [esi-16*dist1K+dist16K]
	sub	ecx, 65536-16
	js	short p283dn
	add	ch, 256/8
	JNC_X	p283lp
	lea	esi, [esi-8*dist16K+dist128K]
	dec	ch
	JNZ_X	p283lp
	mov	ch, 8
	lea	esi, [esi-8*dist128K+dist1M]
	JMP_X	p283lp
p283dn:	ret

;; Do the last eight levels of a FFT multiply (multiply, inverse FFT)

LABELP	pass2_8_levels_type_4
	pass2_eight_levels_real 4
LABELP	pass2_8_levels_type_4p
	mov	ecx, count1		;; Number of complex iterations
	mov	edx, pass2_premults	;; Address of the group multipliers
p284lp:	pass2_eight_levels_complex 4
	add	cl, 256/2
	JNC_X	p284lp
	lea	esi, [esi-8*dist128+dist1K]
	dec	cl
	JNZ_X	p284lp
	lea	esi, [esi-16*dist1K+dist16K]
	sub	ecx, 65536-16
	js	short p284dn
	add	ch, 256/8
	JNC_X	p284lp
	lea	esi, [esi-8*dist16K+dist128K]
	dec	ch
	JNZ_X	p284lp
	mov	ch, 8
	lea	esi, [esi-8*dist128K+dist1M]
	JMP_X	p284lp
p284dn:	ret
pass2_procs ENDP

_TEXT32	ENDS
END
