; Copyright 1995-2000 Just For Fun Software, Inc., all rights reserved
; Author:  George Woltman
; Email: woltman@alum.mit.edu
;
; This routine implements the 6 level 2nd pass for FFTs from 1024 doubles to
; 8192 doubles.
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

	flat_distances

;; Routines to do the last six levels in a two-pass FFT

	PUBLICP	pass2_six_levels_type_123
	PUBLICP	pass2_six_levels_type_123p
	PUBLICP	pass2_six_levels_type_4
	PUBLICP	pass2_six_levels_type_4p

PROCP	pass2_procs2

;; Branch to the proper forward FFT code - mod 2^N+1 arithmetic

LABELP	pass2_six_levels_type_123p
	cmp	ffttype, 2
	JE_X	pass2_six_levels_type_2p
	JG_X	pass2_six_levels_type_3p
	JMP_X	pass2_six_levels_type_1p

;; Branch to the proper forward FFT code

LABELP	pass2_six_levels_type_123
	cmp	ffttype, 2
	JE_X	pass2_six_levels_type_2
	JG_X	pass2_six_levels_type_3
	;; Fall through

;; Do the last six levels of a forward FFT

pass2_six_levels_type_1:
	pass2_six_levels_real 1
pass2_six_levels_type_1p:
	mov	ecx, count1		;; Number of complex iterations
	mov	edx, pass2_premults	;; Address of the group multipliers
p261lp:	pass2_six_levels_complex 1
	dec	cl
	jnz	short p261lq
	lea	esi, [esi-512*dist1+dist512]
	mov	cl, 4
p261lq:	dec	ch
	JNZ_X	p261lp
	pop	esi			;; Pop values and return to C code
	pop	edi
	pop	ebx
	pop	ebp
	ret

;; Do the last six levels of a squaring (forward FFT, square, inverse FFT)

pass2_six_levels_type_2:
	pass2_six_levels_real 2
pass2_six_levels_type_2p:
	mov	ecx, count1		;; Number of complex iterations
	mov	edx, pass2_premults	;; Address of the group multipliers
p262lp:	pass2_six_levels_complex 2
	dec	cl
	jnz	short p262lq
	lea	esi, [esi-512*dist1+dist512]
	mov	cl, 4
p262lq:	dec	ch
	JNZ_X	p262lp
	ret

;; Do the last six levels of a multiply (forward FFT, multiply, inverse FFT)

pass2_six_levels_type_3:
	pass2_six_levels_real 3
pass2_six_levels_type_3p:
	mov	ecx, count1		;; Number of complex iterations
	mov	edx, pass2_premults	;; Address of the group multipliers
p263lp:	pass2_six_levels_complex 3
	dec	cl
	jnz	short p263lq
	lea	esi, [esi-512*dist1+dist512]
	mov	cl, 4
p263lq:	dec	ch
	JNZ_X	p263lp
	ret

;; Do the last six levels of a FFT multiply (multiply, inverse FFT)

LABELP	pass2_six_levels_type_4
	pass2_six_levels_real 4
LABELP	pass2_six_levels_type_4p
	mov	ecx, count1		;; Number of complex iterations
	mov	edx, pass2_premults	;; Address of the group multipliers
p264lp:	pass2_six_levels_complex 4
	dec	cl
	jnz	short p264lq
	lea	esi, [esi-512*dist1+dist512]
	mov	cl, 4
p264lq:	dec	ch
	JNZ_X	p264lp
	ret
ENDPP	pass2_procs2

_TEXT32	ENDS
END
