; Copyright 1995-1999 Just For Fun Software, Inc., all rights reserved
; Author:  George Woltman
; Email: woltman@alum.mit.edu
;
; This routine implements a discrete weighted transform to quickly multiply
; two numbers.
;
; This code handles FFTs that use the very convoluted memory model.
; FFT sizes greater than 64K doubles are supported.
;
; You will not stand a chance of understanding any of this code without
; thoroughly familiarizing yourself with fast fourier transforms.  This
; code was adapted from an algorithm described in Richard Crandall's article
; on Discrete Weighted Transforms and Large-Integer Arithmetic.
;

	TITLE   setup

	.386

_TEXT32 SEGMENT PARA USE32 PUBLIC 'DATA'

	ASSUME  CS: _TEXT32, DS: _TEXT32, SS: _TEXT32, ES: _TEXT32

INCLUDE extrn.mac
INCLUDE	unravel.mac
INCLUDE mult.mac
INCLUDE memory.mac

	very_convoluted_distances

INCLUDE normal.mac
PURGE simple_normalize
PURGE simple_normalize_012
PURGE simple_norm_op
PURGE simple_norm_addsub
PURGE simple_error_check

;;
;; Add two numbers without carry propogation.  Caller can use this for
;; consecutive add or subtract operations.  However, the last operation
;; before a multiply must use the routine that will normalize data.
;;

	PUBLIC	gwaddq4
gwaddq4	PROC NEAR
	mov	ecx, _SRCARG		; U - Address of first number
	mov	edx, _SRC2ARG		; V - Address of second number
	push	esi			; U - Save registers
	push	edi			; V - Save registers
	mov	eax, [ecx-4]		; U - Load needs-normalize counter
	mov	edi, [edx-4]		; V - Load needs-normalize counter
	mov	esi, _DESTARG		; U - Address of destination
	sub	edx, ecx		; V - Compute diff between sources
	lea	eax, [eax+edi+1]	; U - Set new needs-normalize counter
	jmp	short uadd		; Do the unnormalized add
gwaddq4	ENDP


;;
;; Add two numbers.  This will be done without carry propogation
;; (for extra speed) if it can be done safely.
;;

	PUBLIC	gwadd4
gwadd4	PROC NEAR

; See if add can be done without normalization

	mov	ecx, _SRCARG		; U - Address of first number
	mov	edx, _SRC2ARG		; V - Address of second number
	push	esi			; U - Save registers
	push	edi			; V - Save registers
	mov	eax, [ecx-4]		; U - Load needs-normalize counter
	mov	edi, [edx-4]		; V - Load needs-normalize counter
	mov	esi, _DESTARG		; U - Address of destination
	sub	edx, ecx		; V - Compute diff between sources
	lea	eax, [eax+edi+1]	; U - Set new needs-normalize counter
	mov	edi, extra_bits		; V - Load compare value
	cmp	eax, edi		; U - Is normalization needed?
	JG_X	nadd			; V - yes

; Do an unnormalized add
; ecx = src arg #1
; edx = distance to src arg #2
; esi = dest arg
; eax = new needs-normalize counter

uadd:	mov	[esi-4], eax		; U - Store needs-normalize counter
	mov	eax, normcount1		; V - Load loop counter
uaddlp:	fld	QWORD PTR [ecx]		; Load first number
	fadd	QWORD PTR [ecx][edx]	; Add in second number
	fld	QWORD PTR [ecx+dist1]	; Load first number
	fadd	QWORD PTR [ecx+dist1][edx]; Add in second number
	fld	QWORD PTR [ecx+2*dist1]	; Load first number
	fadd	QWORD PTR [ecx+2*dist1][edx]; Add in second number
	fxch	st(1)
	fld	QWORD PTR [ecx+3*dist1]	; Load first number
	fadd	QWORD PTR [ecx+3*dist1][edx]; Add in second number
	fxch	st(3)
	fstp	QWORD PTR [esi]		; Save result
	fstp	QWORD PTR [esi+dist1]	; Save result
	fstp	QWORD PTR [esi+2*dist1]	; Save result
	fstp	QWORD PTR [esi+3*dist1]	; Save result
	fld	QWORD PTR [ecx+4*dist1]	; Load first number
	fadd	QWORD PTR [ecx+4*dist1][edx]; Add in second number
	fld	QWORD PTR [ecx+5*dist1]	; Load first number
	fadd	QWORD PTR [ecx+5*dist1][edx]; Add in second number
	fld	QWORD PTR [ecx+6*dist1]	; Load first number
	fadd	QWORD PTR [ecx+6*dist1][edx]; Add in second number
	fxch	st(1)
	fld	QWORD PTR [ecx+7*dist1]	; Load first number
	fadd	QWORD PTR [ecx+7*dist1][edx]; Add in second number
	fxch	st(3)
	fstp	QWORD PTR [esi+4*dist1]	; Save result
	fstp	QWORD PTR [esi+5*dist1]	; Save result
	fstp	QWORD PTR [esi+6*dist1]	; Save result
	fstp	QWORD PTR [esi+7*dist1]	; Save result
	lea	ecx, [ecx+dist8]	; U - Next source
	lea	esi, [esi+dist8]	; V - Next source
	add	eax, 65536/16*65536	; U - Decrement loop counter
	JNC_X	uaddlp			; V - Loop til done
	lea	ecx, [ecx-16*dist8+dist128]; U - Next source
	lea	esi, [esi-16*dist8+dist128]; V - Next source
	test	esi, 31			; U - Test loop counter
	JNZ_X	uaddlp			; V - Loop til done
	sub	eax, 32768-1		; U - Decrement loop counter
	js	short uadddn		; V - Loop til done
	test	eax, 2-1		;UV - Test next loop counter
	JNZ_X	uaddlp			;*V - Loop til done
	test	eax, 16*2-1		;UV - Test next loop counter
	lea	ecx, [ecx-8*dist128+dist1K]; U - Next source
	lea	esi, [esi-8*dist128+dist1K]; V - Next source
	JNZ_X	uaddlp			; V - Loop til done
	test	eax, 8*16*2-1		;UV - Test next loop counter
	lea	ecx, [ecx-16*dist1K+dist16K]; U - Next source
	lea	esi, [esi-16*dist1K+dist16K]; V - Next source
	JNZ_X	uaddlp			; V - Loop til done
	test	eax, 8*8*16*2-1		;UV - Test next loop counter
	lea	ecx, [ecx-8*dist16K+dist128K]; U - Next source
	lea	esi, [esi-8*dist16K+dist128K]; V - Next source
	JNZ_X	uaddlp			; V - Loop til done
	lea	ecx, [ecx-8*dist128K+dist1M]; U - Next source
	lea	esi, [esi-8*dist128K+dist1M]; U - Next source
	JMP_X	uaddlp			; V - Loop til done
uadddn:	pop	edi
	pop	esi
	ret

; Do a normalized add
; ecx = src arg #1
; edx = distance to src arg #2
; esi = dest arg

nadd:	fninit
	push	ebp			; U - Save registers
	push	ebx			; V - Save registers
	fld	BIGVAL			; UV - Start process with no carry
	mov	edi, norm_grp_mults	; U - Address of group multipliers
	mov	ebp, scaled_numlit	; V - Computes big vs little words
	mov	eax, -1			; U - First word is a big word (-1)
	mov	ebx, normcount1		; V - Load loop counter
naddlp2:fld	QWORD PTR [edi]		; UV - Two-to-phi
	fld	QWORD PTR [edi+8]	; UV - Two-to-minus-phi group mult
	fmul	ttmp_ff			; UV - Apply two-to-minus-phi fudge
	fld	QWORD PTR [edi]		; UV - Two-to-phi group multiplier
	lea	edi, [edi+NMD]		; U - Next group multiplier
	mov	normgrpptr, edi		;*U - Save group multiplier address
	mov	edi, norm_col_mults	; V - Address of the multipliers
	fld	st(1)			; UV - Two-to-minus-phi
naddlp:	norm_op	fadd dist8		; Add and normalize 8 values
	add	ebx, 65536/16*65536	; U - Decrement loop counter
	JNC_X	naddlp			; V - Loop til done
	lea	ecx, [ecx-16*dist8+dist128]; U - Next source
	lea	esi, [esi-16*dist8+dist128]; V - Next dest
	mov	edi, normgrpptr		; U - Restore group multiplier address
	fcompp				; UV - Pop group multipliers
	fcompp				; UV - Pop group multipliers
	test	esi, 31			; U - Test loop counter
	JNZ_X	naddlp2			; V - Loop til done
	sub	ebx, 32768-1		; U - Decrement loop counter
	js	short nadddn		; V - Loop til done
	sub	ebp, scaling_ff		;UU - Make edi accurate again
	test	ebx, 2-1		;UV - Test next loop counter
	JNZ_X	naddlp2			;*V - Loop til done
	test	ebx, 16*2-1		;UV - Test next loop counter
	lea	ecx, [ecx-8*dist128+dist1K]; U - Next source
	lea	esi, [esi-8*dist128+dist1K]; V - Next source
	JNZ_X	naddlp2			; V - Loop til done
	test	ebx, 8*16*2-1		;UV - Test next loop counter
	lea	ecx, [ecx-16*dist1K+dist16K]; U - Next source
	lea	esi, [esi-16*dist1K+dist16K]; V - Next source
	JNZ_X	naddlp2			; V - Loop til done
	sub	ebp, scaling_ff2	;UU - Make edi accurate again
	test	ebx, 8*8*16*2-1		;UV - Test next loop counter
	lea	ecx, [ecx-8*dist16K+dist128K]; U - Next source
	lea	esi, [esi-8*dist16K+dist128K]; V - Next source
	JNZ_X	naddlp2			; V - Loop til done
	lea	ecx, [ecx-8*dist128K+dist1M]; U - Next source
	lea	esi, [esi-8*dist128K+dist1M]; V - Next source
	JMP_X	naddlp2			; V - Loop til done
nadddn:	fsub	BIGVAL			; Compute carry
	mov	ecx, _DESTARG		; Address of result
	cmp	_PLUS1, 0		; Change sign of carry
	jz	short nadd1		; if we are doing
	fchs				; mod 2^N+1 operations
nadd1:	fadd	QWORD PTR [ecx]		; Add in carry
	mov	DWORD PTR [ecx-4], 0	; Clear needs-normalize flag
	pop	ebx			; Restore registers
	pop	ebp
	pop	edi
	pop	esi
	fstp	QWORD PTR [ecx]		; Save new value
	ret
gwadd4	ENDP


;;
;; Subtract two numbers without carry propogation.  Caller can use this for
;; consecutive add or subtract operations.  However, the last operation
;; before a multiply must use the routine that will normalize data.
;;

	PUBLIC	gwsubq4
gwsubq4	PROC NEAR
	mov	ecx, _SRCARG		; U - Address of first number
	mov	edx, _SRC2ARG		; V - Address of second number
	push	esi			; U - Save registers
	push	edi			; V - Save registers
	mov	eax, [ecx-4]		; U - Load needs-normalize counter
	mov	edi, [edx-4]		; V - Load needs-normalize counter
	mov	esi, _DESTARG		; U - Address of destination
	sub	edx, ecx		; V - Compute diff between sources
	lea	eax, [eax+edi+1]	; U - Set new needs-normalize counter
	jmp	short usub		; Do the unnormalized sub
gwsubq4	ENDP

;;
;; Subtract two numbers.  This will be done without carry propogation
;; (for extra speed) if it can be done safely.
;;

	PUBLIC	gwsub4
gwsub4	PROC NEAR

; See if subtract can be done without normalization

	mov	ecx, _SRCARG		; U - Address of first number
	mov	edx, _SRC2ARG		; V - Address of second number
	push	esi			; U - Save registers
	push	edi			; V - Save registers
	mov	eax, [ecx-4]		; U - Load needs-normalize counter
	mov	edi, [edx-4]		; V - Load needs-normalize counter
	mov	esi, _DESTARG		; U - Address of destination
	sub	edx, ecx		; V - Compute diff between sources
	lea	eax, [eax+edi+1]	; U - Set new needs-normalize counter
	mov	edi, extra_bits		; V - Load compare value
	cmp	eax, edi		; U - Is normalization needed?
	JG_X	nsub			; V - yes

; Do an unnormalized subtract
; ecx = src arg #1
; edx = distance to src arg #2
; esi = dest arg
; eax = new needs-normalize counter

usub:	mov	[esi-4], eax		; U - Store needs-normalize counter
	mov	eax, normcount1		; V - Load loop counter
usublp:	fld	QWORD PTR [ecx][edx]	; Load second number
	fsub	QWORD PTR [ecx]		; Subtract first number
	fld	QWORD PTR [ecx+dist1][edx]; Load second number
	fsub	QWORD PTR [ecx+dist1]	; Subtract first number
	fld	QWORD PTR [ecx+2*dist1][edx]; Load second number
	fsub	QWORD PTR [ecx+2*dist1]	; Subtract first number
	fxch	st(1)
	fld	QWORD PTR [ecx+3*dist1][edx]; Load second number
	fsub	QWORD PTR [ecx+3*dist1]	; Subtract first number
	fxch	st(3)
	fstp	QWORD PTR [esi]		; Save result
	fstp	QWORD PTR [esi+dist1]	; Save result
	fstp	QWORD PTR [esi+2*dist1]	; Save result
	fstp	QWORD PTR [esi+3*dist1]	; Save result
	fld	QWORD PTR [ecx+4*dist1][edx]; Load second number
	fsub	QWORD PTR [ecx+4*dist1]	; Subtract first number
	fld	QWORD PTR [ecx+5*dist1][edx]; Load second number
	fsub	QWORD PTR [ecx+5*dist1]	; Subtract first number
	fld	QWORD PTR [ecx+6*dist1][edx]; Load second number
	fsub	QWORD PTR [ecx+6*dist1]	; Subtract first number
	fxch	st(1)
	fld	QWORD PTR [ecx+7*dist1][edx]; Load second number
	fsub	QWORD PTR [ecx+7*dist1]	; Subtract first number
	fxch	st(3)
	fstp	QWORD PTR [esi+4*dist1]	; Save result
	fstp	QWORD PTR [esi+5*dist1]	; Save result
	fstp	QWORD PTR [esi+6*dist1]	; Save result
	fstp	QWORD PTR [esi+7*dist1]	; Save result
	lea	ecx, [ecx+dist8]	; U - Next source
	lea	esi, [esi+dist8]	; V - Next source
	add	eax, 65536/16*65536	; U - Decrement loop counter
	JNC_X	usublp			; V - Loop til done
	lea	ecx, [ecx-16*dist8+dist128]; U - Next source
	lea	esi, [esi-16*dist8+dist128]; V - Next source
	test	esi, 31			; U - Test loop counter
	JNZ_X	usublp			; V - Loop til done
	sub	eax, 32768-1		; U - Decrement loop counter
	js	short usubdn		; V - Loop til done
	test	eax, 2-1		;UV - Test next loop counter
	JNZ_X	usublp			;*V - Loop til done
	test	eax, 16*2-1		;UV - Test next loop counter
	lea	ecx, [ecx-8*dist128+dist1K]; U - Next source
	lea	esi, [esi-8*dist128+dist1K]; V - Next source
	JNZ_X	usublp			; V - Loop til done
	test	eax, 8*16*2-1		;UV - Test next loop counter
	lea	ecx, [ecx-16*dist1K+dist16K]; U - Next source
	lea	esi, [esi-16*dist1K+dist16K]; V - Next source
	JNZ_X	usublp			; V - Loop til done
	test	eax, 8*8*16*2-1		;UV - Test next loop counter
	lea	ecx, [ecx-8*dist16K+dist128K]; U - Next source
	lea	esi, [esi-8*dist16K+dist128K]; V - Next source
	JNZ_X	usublp			; V - Loop til done
	lea	ecx, [ecx-8*dist128K+dist1M]; U - Next source
	lea	esi, [esi-8*dist128K+dist1M]; U - Next source
	JMP_X	usublp			; V - Loop til done
usubdn:	pop	edi
	pop	esi
	ret

; Do a normalized subtract
; ecx = src arg #1
; edx = distance to src arg #2
; esi = dest arg

nsub:	fninit
	push	ebp			; U - Save registers
	push	ebx			; V - Save registers
	fld	BIGVAL			; UV - Start process with no carry
	mov	edi, norm_grp_mults	; U - Address of group multipliers
	mov	ebp, scaled_numlit	; V - Computes big vs little words
	mov	eax, -1			; U - First word is a big word (-1)
	mov	ebx, normcount1		; V - Load loop counter
nsublp2:fld	QWORD PTR [edi]		; UV - Two-to-phi
	fld	QWORD PTR [edi+8]	; UV - Two-to-minus-phi group mult
	fmul	ttmp_ff			; UV - Apply two-to-minus-phi fudge
	fld	QWORD PTR [edi]		; UV - Two-to-phi group multiplier
	lea	edi, [edi+NMD]		; U - Next group multiplier
	mov	normgrpptr, edi		;*U - Save group multiplier address
	mov	edi, norm_col_mults	; V - Address of the multipliers
	fld	st(1)			; UV - Two-to-minus-phi
nsublp:	norm_op	fsub dist8		; Subtract and normalize 8 values
	add	ebx, 65536/16*65536	; U - Decrement loop counter
	JNC_X	nsublp			; V - Loop til done
	lea	ecx, [ecx-16*dist8+dist128]; U - Next source
	lea	esi, [esi-16*dist8+dist128]; V - Next dest
	mov	edi, normgrpptr		; U - Restore group multiplier address
	fcompp				; UV - Pop group multipliers
	fcompp				; UV - Pop group multipliers
	test	esi, 31			; U - Test loop counter
	JNZ_X	nsublp2			; V - Loop til done
	sub	ebx, 32768-1		; U - Decrement loop counter
	js	short nsubdn		; V - Loop til done
	sub	ebp, scaling_ff		;UU - Make edi accurate again
	test	ebx, 2-1		;UV - Test next loop counter
	JNZ_X	nsublp2			;*V - Loop til done
	test	ebx, 16*2-1		;UV - Test next loop counter
	lea	ecx, [ecx-8*dist128+dist1K]; U - Next source
	lea	esi, [esi-8*dist128+dist1K]; V - Next source
	JNZ_X	nsublp2			; V - Loop til done
	test	ebx, 8*16*2-1		;UV - Test next loop counter
	lea	ecx, [ecx-16*dist1K+dist16K]; U - Next source
	lea	esi, [esi-16*dist1K+dist16K]; V - Next source
	JNZ_X	nsublp2			; V - Loop til done
	sub	ebp, scaling_ff2	;UU - Make edi accurate again
	test	ebx, 8*8*16*2-1		;UV - Test next loop counter
	lea	ecx, [ecx-8*dist16K+dist128K]; U - Next source
	lea	esi, [esi-8*dist16K+dist128K]; V - Next source
	JNZ_X	nsublp2			; V - Loop til done
	lea	ecx, [ecx-8*dist128K+dist1M]; U - Next source
	lea	esi, [esi-8*dist128K+dist1M]; V - Next source
	JMP_X	nsublp2			; V - Loop til done
nsubdn:	fsub	BIGVAL			; Compute carry
	mov	ecx, _DESTARG		; Address of result
	cmp	_PLUS1, 0		; Change sign of carry
	jz	short nsub1		; if we are doing
	fchs				; mod 2^N+1 operations
nsub1:	fadd	QWORD PTR [ecx]		; Add in carry
	mov	DWORD PTR [ecx-4], 0	; Clear needs-normalize flag
	pop	ebx			; Restore registers
	pop	ebp
	pop	edi
	pop	esi
	fstp	QWORD PTR [ecx]		; Save new value
	ret
gwsub4	ENDP

;;
;; Add and subtract two numbers without carry propogation.
;;

	PUBLIC	gwaddsubq4
gwaddsubq4 PROC NEAR
	mov	ecx, _SRCARG		; U - Address of first number
	mov	edx, _SRC2ARG		; V - Address of second number
	push	esi			; U - Save registers
	push	edi			; V - Save registers
	mov	eax, [ecx-4]		; U - Load needs-normalize counter
	mov	edi, [edx-4]		; V - Load needs-normalize counter
	push	ebp			; U - Save registers
	mov	esi, _DESTARG		; V - Address of destination
	mov	ebp, _DEST2ARG  	; U - Address of destination #2
	sub	edx, ecx		; V - Compute diff between sources
	lea	eax, [eax+edi+1]	; U - Set new needs-normalize counter
	jmp	short uaddsub		; Do the unnormalized add
gwaddsubq4 ENDP

;;
;; Add and subtract two numbers.  This will be done without carry propogation
;; (for extra speed) if it can be done safely.
;;

	PUBLIC	gwaddsub4
gwaddsub4 PROC NEAR

; See if add can be done without normalization

	mov	ecx, _SRCARG		; U - Address of first number
	mov	edx, _SRC2ARG		; V - Address of second number
	push	esi			; U - Save registers
	push	edi			; V - Save registers
	mov	eax, [ecx-4]		; U - Load needs-normalize counter
	mov	edi, [edx-4]		; V - Load needs-normalize counter
	push	ebp			; U - Save registers
	mov	esi, _DESTARG		; V - Address of destination
	mov	ebp, _DEST2ARG  	; U - Address of destination #2
	sub	edx, ecx		; V - Compute diff between sources
	lea	eax, [eax+edi+1]	; U - Set new needs-normalize counter
	mov	edi, extra_bits		; V - Load compare value
	cmp	eax, edi		; U - Is normalization needed?
	JG_X	naddsub			; V - yes

; Do an unnormalized add and subtract
; ecx = src arg #1
; edx = distance to src arg #2
; esi = dest arg #1
; ebp = dest arg #2
; eax = new needs-normalize counter

uaddsub:mov	[esi-4], eax		; U - Store needs-normalize counter
	mov	[ebp-4], eax		; V - Store needs-normalize counter
	mov	eax, normcount1		; U - Load loop counter
	sub	ebp, esi		; V - Compute diff between destinations
uaddsublp:
	fld	QWORD PTR [ecx]		; Load first number
	fld	st(0)			; Dup first number
	fadd	QWORD PTR [ecx][edx]	; Add in second number
	fxch	st(1)			; S0,A0
	fsub	QWORD PTR [ecx][edx]	; Subtract out second number
	fld	QWORD PTR [ecx+dist1]	; Load first number
	fld	st(0)			; Dup first number
	fadd	QWORD PTR [ecx+dist1][edx]; Add in second number
	fxch	st(1)			; S1,A1,S0,A0
	fsub	QWORD PTR [ecx+dist1][edx]; Subtract out second number
	fld	QWORD PTR [ecx+2*dist1]	; Load first number
	fld	st(0)			; Dup first number
	fadd	QWORD PTR [ecx+2*dist1][edx]; Add in second number
	fxch	st(1)			; S2,A2,S1,A1,S0,A0
	fsub	QWORD PTR [ecx+2*dist1][edx]; Subtract out second number
	fld	QWORD PTR [ecx+3*dist1]	; Load first number
	fld	st(0)			; Dup first number
	fadd	QWORD PTR [ecx+3*dist1][edx]; Add in second number
	fxch	st(7)			; A0,S3,S2,A2,S1,A1,S0,A3
	fstp	QWORD PTR [esi]		; Save result
	fsub	QWORD PTR [ecx+3*dist1][edx]; Subtract out second number
	fxch	st(5)			; S0,S2,A2,S1,A1,S3,A3
	fstp	QWORD PTR [esi][ebp]	; Save result
	fstp	QWORD PTR [esi+2*dist1][ebp]; Save result
	fstp	QWORD PTR [esi+2*dist1]	; Save result
	fstp	QWORD PTR [esi+dist1][ebp]; Save result
	fstp	QWORD PTR [esi+dist1]	; Save result
	fstp	QWORD PTR [esi+3*dist1][ebp]; Save result
	fstp	QWORD PTR [esi+3*dist1]	; Save result
	fld	QWORD PTR [ecx+4*dist1]	; Load first number
	fld	st(0)			; Dup first number
	fadd	QWORD PTR [ecx+4*dist1][edx]; Add in second number
	fxch	st(1)			; S0,A0
	fsub	QWORD PTR [ecx+4*dist1][edx]; Subtract out second number
	fld	QWORD PTR [ecx+5*dist1]	; Load first number
	fld	st(0)			; Dup first number
	fadd	QWORD PTR [ecx+5*dist1][edx]; Add in second number
	fxch	st(1)			; S1,A1,S0,A0
	fsub	QWORD PTR [ecx+5*dist1][edx]; Subtract out second number
	fld	QWORD PTR [ecx+6*dist1]	; Load first number
	fld	st(0)			; Dup first number
	fadd	QWORD PTR [ecx+6*dist1][edx]; Add in second number
	fxch	st(1)			; S2,A2,S1,A1,S0,A0
	fsub	QWORD PTR [ecx+6*dist1][edx]; Subtract out second number
	fld	QWORD PTR [ecx+7*dist1]	; Load first number
	fld	st(0)			; Dup first number
	fadd	QWORD PTR [ecx+7*dist1][edx]; Add in second number
	fxch	st(7)			; A0,S3,S2,A2,S1,A1,S0,A3
	fstp	QWORD PTR [esi+4*dist1]	; Save result
	fsub	QWORD PTR [ecx+7*dist1][edx]; Subtract out second number
	fxch	st(5)			; S0,S2,A2,S1,A1,S3,A3
	fstp	QWORD PTR [esi+4*dist1][ebp]; Save result
	fstp	QWORD PTR [esi+6*dist1][ebp]; Save result
	fstp	QWORD PTR [esi+6*dist1]	; Save result
	fstp	QWORD PTR [esi+5*dist1][ebp]; Save result
	fstp	QWORD PTR [esi+5*dist1]	; Save result
	fstp	QWORD PTR [esi+7*dist1][ebp]; Save result
	fstp	QWORD PTR [esi+7*dist1]	; Save result
	lea	ecx, [ecx+dist8]	; U - Next source
	lea	esi, [esi+dist8]	; V - Next source
	add	eax, 65536/16*65536	; U - Decrement loop counter
	JNC_X	uaddsublp		; V - Loop til done
	lea	ecx, [ecx-16*dist8+dist128]; U - Next source
	lea	esi, [esi-16*dist8+dist128]; V - Next source
	test	esi, 31			; U - Test loop counter
	JNZ_X	uaddsublp		; V - Loop til done
	sub	eax, 32768-1		; U - Decrement loop counter
	js	short uaddsubdn		; V - Loop til done
	test	eax, 2-1		;UV - Test next loop counter
	JNZ_X	uaddsublp		;*V - Loop til done
	test	eax, 16*2-1		;UV - Test next loop counter
	lea	ecx, [ecx-8*dist128+dist1K]; U - Next source
	lea	esi, [esi-8*dist128+dist1K]; V - Next source
	JNZ_X	uaddsublp		; V - Loop til done
	test	eax, 8*16*2-1		;UV - Test next loop counter
	lea	ecx, [ecx-16*dist1K+dist16K]; U - Next source
	lea	esi, [esi-16*dist1K+dist16K]; V - Next source
	JNZ_X	uaddsublp		; V - Loop til done
	test	eax, 8*8*16*2-1		;UV - Test next loop counter
	lea	ecx, [ecx-8*dist16K+dist128K]; U - Next source
	lea	esi, [esi-8*dist16K+dist128K]; V - Next source
	JNZ_X	uaddsublp		; V - Loop til done
	lea	ecx, [ecx-8*dist128K+dist1M]; U - Next source
	lea	esi, [esi-8*dist128K+dist1M]; U - Next source
	JMP_X	uaddsublp		; V - Loop til done
uaddsubdn:
	pop	ebp
	pop	edi
	pop	esi
	ret

; Do a normalized add and subtract
; ecx = src arg #1
; edx = distance to src arg #2
; esi = dest arg #1
; ebp = dest arg #2

naddsub:fninit
	push	ebx			; U - Save registers
	sub	ebp, esi		; V - Compute diff between destinations
	fld	BIGVAL			; UV - Start process with no carry
	fld	BIGVAL			; UV - Start process with no carry
	mov	addsubtemp, edx		; U - Save edx
	mov	edx, normcount1		; V - Load loop counter
	mov	edi, norm_grp_mults	; U - Address of group multipliers
	mov	ebx, scaled_numlit	; V - Computes big vs little words
	mov	eax, -1			; U - First word is a big word (-1)
naddsublp2:
	fld	QWORD PTR [edi]		; UV - Two-to-phi
	fld	QWORD PTR [edi+8]	; UV - Two-to-minus-phi group mult
	fmul	ttmp_ff			; UV - Apply two-to-minus-phi fudge
	fld	QWORD PTR [edi]		; UV - Two-to-phi group multiplier
	lea	edi, [edi+NMD]		; U - Next group multiplier
	mov	normgrpptr, edi		;*U - Save group multiplier address
	mov	edi, norm_col_mults	; V - Address of the multipliers
	fld	st(1)			; UV - Two-to-minus-phi
naddsublp:
	push	edx			; U - Save loop counter
	mov	edx, addsubtemp		; V - Restore source #2 distance
	norm_addsub dist8		; Add and normalize 8 values
	add	edx, 65536/16*65536	; U - Decrement loop counter
	JNC_X	naddsublp		; V - Loop til done
	lea	ecx, [ecx-16*dist8+dist128]; U - Next source
	lea	esi, [esi-16*dist8+dist128]; V - Next dest
	mov	edi, normgrpptr		; U - Restore group multiplier address
	fcompp				; UV - Pop group multipliers
	fcompp				; UV - Pop group multipliers
	test	esi, 31			; U - Test loop counter
	JNZ_X	naddsublp2		; V - Loop til done
	sub	edx, 32768-1		; U - Decrement loop counter
	js	short naddsubdn		; V - Loop til done
	sub	ebx, scaling_ff		;UU - Make ebx accurate again
	test	edx, 2-1		;UV - Test next loop counter
	JNZ_X	naddsublp2		;*V - Loop til done
	test	edx, 16*2-1		;UV - Test next loop counter
	lea	ecx, [ecx-8*dist128+dist1K]; U - Next source
	lea	esi, [esi-8*dist128+dist1K]; V - Next source
	JNZ_X	naddsublp2		; V - Loop til done
	test	edx, 8*16*2-1		;UV - Test next loop counter
	lea	ecx, [ecx-16*dist1K+dist16K]; U - Next source
	lea	esi, [esi-16*dist1K+dist16K]; V - Next source
	JNZ_X	naddsublp2		; V - Loop til done
	sub	ebx, scaling_ff2	;UU - Make edi accurate again
	test	edx, 8*8*16*2-1		;UV - Test next loop counter
	lea	ecx, [ecx-8*dist16K+dist128K]; U - Next source
	lea	esi, [esi-8*dist16K+dist128K]; V - Next source
	JNZ_X	naddsublp2		; V - Loop til done
	lea	ecx, [ecx-8*dist128K+dist1M]; U - Next source
	lea	esi, [esi-8*dist128K+dist1M]; V - Next source
	JMP_X	naddsublp2		; V - Loop til done
naddsubdn:
	fsub	BIGVAL			; Compute carry
	fxch	st(1)
	fsub	BIGVAL			; Compute carry
	mov	ecx, _DESTARG		; Address of result
	mov	eax, _DEST2ARG		; Address of result
	cmp	_PLUS1, 0		; Change sign of carry
	jz	short naddsub1		; if we are doing
	fchs				; mod 2^N+1 operations
	fxch	st(1)
	fchs				; UV - Yes, flip carry's sign bit
	fxch	st(1)
naddsub1:
	fadd	QWORD PTR [ecx]		; Add in carry
	fxch	st(1)
	fadd	QWORD PTR [eax]		; Add in carry
	mov	DWORD PTR [ecx-4], 0	; Clear needs-normalize flag
	mov	DWORD PTR [eax-4], 0	; Clear needs-normalize flag
	pop	ebx			; Restore registers
	pop	ebp
	pop	edi
	pop	esi
	fstp	QWORD PTR [eax]		; Save new value
	fstp	QWORD PTR [ecx]		; Save new value
	ret
gwaddsub4 ENDP

_TEXT32	ENDS
END
