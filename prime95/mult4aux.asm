; Copyright 1995-2001 Just For Fun Software, Inc., all rights reserved
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

	.686
	.XMM

_TEXT32 SEGMENT PARA USE32 PUBLIC 'DATA'

	ASSUME  CS: _TEXT32, DS: _TEXT32, SS: _TEXT32, ES: _TEXT32

INCLUDE extrn.mac
INCLUDE	unravel.mac
INCLUDE mult.mac
INCLUDE memory.mac

	very_convoluted_distances

INCLUDE normal.mac

;;
;; Add two numbers without carry propogation.  Caller can use this for
;; consecutive add or subtract operations.  However, the last operation
;; before a multiply must use the routine that will normalize data.
;;

IFNDEF PPRO
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
naddlp:	norm_op_2d fadd, dist8		; Add and normalize 8 values
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
nsublp:	norm_op_2d fsub, dist8		; Subtract and normalize 8 values
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
	mov	loopcount1, edx		; U - Save edx
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
	mov	edx, loopcount1		; V - Restore source #2 distance
	norm_addsub_2d dist8		; Add and normalize 8 values
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

;;
;; Copy one number and zero some low order words.
;;

	PUBLIC	gwcopyzero4
gwcopyzero4 PROC NEAR
	push	esi			; U - Save registers
	mov	esi, _SRCARG		; V - Address of first number
	push	edi			; U - Save registers
	mov	edi, _DESTARG		; V - Address of destination
	mov	eax, [esi-4]		; U - Load needs-normalize counter
	mov	ecx, _SRC2ARG		; V - Address of first number
	mov	[edi-4], eax		; U - Store needs-normalize counter
	mov	eax, normcount1		; V - Load loop counter
zlp:	dec	ecx			; U - Decrement count of words to zero
	js	short c8		; V - If negative we are now copying
	fldz
	dec	ecx			; U - Decrement count of words to zero
	js	short c7		; V - If negative we are now copying
	fldz
	dec	ecx			; U - Decrement count of words to zero
	js	short c6		; V - If negative we are now copying
	fldz
	dec	ecx			; U - Decrement count of words to zero
	js	short c5		; V - If negative we are now copying
	fldz
	dec	ecx			; U - Decrement count of words to zero
	js	short c4		; V - If negative we are now copying
	fldz
	dec	ecx			; U - Decrement count of words to zero
	js	short c3		; V - If negative we are now copying
	fldz
	dec	ecx			; U - Decrement count of words to zero
	js	short c2		; V - If negative we are now copying
	fldz
	dec	ecx			; U - Decrement count of words to zero
	js	short c1		; V - If negative we are now copying
	fldz
	jmp	short c0
c8:	fld	QWORD PTR [esi]		; Load number
c7:	fld	QWORD PTR [esi+dist1]	; Load number
c6:	fld	QWORD PTR [esi+2*dist1]	; Load number
c5:	fld	QWORD PTR [esi+3*dist1]	; Load number
c4:	fld	QWORD PTR [esi+4*dist1]	; Load number
c3:	fld	QWORD PTR [esi+5*dist1]	; Load number
c2:	fld	QWORD PTR [esi+6*dist1]	; Load number
c1:	fld	QWORD PTR [esi+7*dist1]	; Load number
c0:	fstp	QWORD PTR [edi+7*dist1]	; Save result
	fstp	QWORD PTR [edi+6*dist1]	; Save result
	fstp	QWORD PTR [edi+5*dist1]	; Save result
	fstp	QWORD PTR [edi+4*dist1]	; Save result
	fstp	QWORD PTR [edi+3*dist1]	; Save result
	fstp	QWORD PTR [edi+2*dist1]	; Save result
	fstp	QWORD PTR [edi+dist1]	; Save result
	fstp	QWORD PTR [edi]		; Save result
	lea	esi, [esi+dist8]	; U - Bump source pointer
	lea	edi, [edi+dist8]	; V - Bump dest pointer
	add	eax, 65536/16*65536	; U - Decrement loop counter
	JNC_X	zlp			; V - Loop til done
	lea	esi, [esi-16*dist8+dist128]; U - Next source
	lea	edi, [edi-16*dist8+dist128]; V - Next dest
	test	esi, 31			; U - Test loop counter
	JNZ_X	zlp			; V - Loop til done
	sub	eax, 32768-1		; U - Decrement loop counter
	js	short zlpdn		; V - Loop til done
	test	eax, 2-1		;UV - Test next loop counter
	JNZ_X	zlp			;*V - Loop til done
	test	eax, 16*2-1		;UV - Test next loop counter
	lea	esi, [esi-8*dist128+dist1K]; U - Next source
	lea	edi, [edi-8*dist128+dist1K]; V - Next source
	JNZ_X	zlp			; V - Loop til done
	test	eax, 8*16*2-1		;UV - Test next loop counter
	lea	esi, [esi-16*dist1K+dist16K]; U - Next source
	lea	edi, [edi-16*dist1K+dist16K]; V - Next source
	JNZ_X	zlp			; V - Loop til done
	test	eax, 8*8*16*2-1		;UV - Test next loop counter
	lea	esi, [esi-8*dist16K+dist128K]; U - Next source
	lea	edi, [edi-8*dist16K+dist128K]; V - Next source
	JNZ_X	zlp			; V - Loop til done
	lea	esi, [esi-8*dist128K+dist1M]; U - Next source
	lea	edi, [edi-8*dist128K+dist1M]; V - Next source
	JMP_X	zlp			; V - Loop til done
zlpdn:	pop	edi
	pop	esi
	ret
gwcopyzero4 ENDP

;;
;; Do a mod k*2^n+/-1
;;

	PUBLIC	gwprothmod4
gwprothmod4 PROC NEAR
	push	esi			; U - Save registers
	push	edi			; V - Save registers
	mov	esi, _SRCARG		; U - Address of first number
	prothmod_upper_prep_0d		; Prepare for prothmod
	mov	ecx, _FFTLEN		; U - Load loop counter
pmulp:	prothmod_upper_0d pmudn		; Divide upper values by k
	lea	esi, [esi-dist8]	; U - Next source
	lea	edi, [edi+8*8]		; V - Next scratch words
	add	ecx, 65536/16*65536	; Decrement loop counter
	JNC_X	pmulp			; Loop til done
	lea	esi, [esi+16*dist8-dist128];; Next source
	sub	ecx, 128		; Decrement loop counter
	test	ecx, 1023		; Test next loop counter
	JNZ_X	pmulp			; Loop til done
	lea	esi, [esi+8*dist128-dist1K]; Next source
	test	ecx, 16383		; Test next loop counter
	JNZ_X	pmulp			; Loop til done
	lea	esi, [esi+16*dist1K-dist16K]; Next source
	test	ecx, 131071		; Test next loop counter
	JNZ_X	pmulp			; Loop til done
	lea	esi, [esi+8*dist16K-dist128K]; Next source
	test	ecx, 1048575		; Test next loop counter
	JNZ_X	pmulp			; Loop til done
	lea	esi, [esi+8*dist128K-dist1M]; Next source
	JMP_X	pmulp			; Loop til done
pmudn:	prothmod_lower_prep_0d		; Prepare for prothmod
	mov	ecx, normcount1		; Load loop counter
pmllp:	prothmod_lower_0d		; Divide upper values by k
	dec	edx
	jz	short pmldn
	lea	esi, [esi+dist8]	; U - Next source ptr
	lea	edi, [edi-8*8]		; V - Next scratch words
	add	ecx, 65536/16*65536	;; Decrement loop counter
	JNC_X	pmllp			;; Loop til done
	lea	esi, [esi-16*dist8+dist128];; Next source
	test	esi, 31			;; Test loop counter
	JNZ_X	pmllp			;; Loop til done
	sub	ecx, 32768-1		;; Decrement loop counter
	test	ecx, 2-1		;; Test next loop counter
	JNZ_X	pmllp			;; Loop til done
	test	ecx, 16*2-1		;; Test next loop counter
	lea	esi, [esi-8*dist128+dist1K];; Next source
	JNZ_X	pmllp			;; Loop til done
	test	ecx, 8*16*2-1		;; Test next loop counter
	lea	esi, [esi-16*dist1K+dist16K];; Next source
	JNZ_X	pmllp			;; Loop til done
	test	ecx, 8*8*16*2-1		;; Test next loop counter
	lea	esi, [esi-8*dist16K+dist128K];; Next source
	JNZ_X	pmllp			;; Loop til done
	lea	esi, [esi-8*dist128K+dist1M];; Next source
	JMP_X	pmllp			;; Loop til done
pmldn:	prothmod_final_0d		; Finish off the mod
	pop	edi			; U - Restore registers
	pop	esi			; V - Restore registers
	ret
gwprothmod4 ENDP
ENDIF

;;
;; Routines to do the normalization after a multiply
;;

PROCP	_comp8norm4

; Macro to loop through all the FFT values and apply the proper normalization
; routine.  Used whenever we are using an irrational-base FFT.

inorm	MACRO	lab, norm_macro
	LOCAL	ilp2, ilp3, ilp4, pf1
	PUBLICP	lab
	LABELP	lab
	mov	esi, _DESTARG		;; Address of multiplied number
start_timer 30
	mov	edi, _ADDIN_OFFSET	;; Get address to add value into
	fld	QWORD PTR [esi][edi]	;; Get the value
	fadd	_ADDIN_VALUE		;; Add in the requested value
	fstp	QWORD PTR [esi][edi]	;; Save the new value
	fldz				;; Init SUMOUT
	fsub	_ADDIN_VALUE		;; Do not include addin in sumout
	fld	BIGVAL			;; Start process with no carry
	mov	ebx, norm_grp_mults	;; Address of the group multipliers
	mov	edi, scaled_numlit	;; Used to compute big vs little words
	mov	ebp, scaled_numbig
	mov	eax, -1			;; First word is a big word (-1)
	mov	ecx, normcount1		;; Load loop counter
ilp4:
	IFDEF PFETCH
	lea	edx, [esi+4*dist128]	;; Assume next block is 4*dist128 away
	test	cl, 2-1			;; Test loop counter
	jz	short ilp3		;; No TLB loads needed
	mov	al, 16*2-1		;; Mask for loop counter
	and	al, cl
	lea	edx, [esi+dist1K-4*dist128];; Assume next block is dist1K away
	cmp	al, 16*2-1		;; Test loop counter again
	jne	short pf1		;; Skip if assumption correct
	lea	edx, [esi+dist16K-15*dist1K-4*dist128];; Next is dist16K away
pf1:	mov	al, [edx+0*dist8]	;; Read from this page (loads the TLB)
	mov	al, [edx+1*dist8]
	mov	al, [edx+2*dist8]
	mov	al, [edx+3*dist8]
	mov	al, [edx+4*dist8]
	mov	al, [edx+5*dist8]
	mov	al, [edx+6*dist8]
	mov	al, [edx+7*dist8]
	mov	al, [edx+8*dist8]
	mov	al, [edx+9*dist8]
	mov	al, [edx+10*dist8]
	mov	al, [edx+11*dist8]
	mov	al, [edx+12*dist8]
	mov	al, [edx+13*dist8]
	mov	al, [edx+14*dist8]
	mov	al, [edx+15*dist8]
	sar	eax, 8			;; Restore eax to zero or minus one
	ENDIF
ilp3:	fld	QWORD PTR [ebx]		;; Two-to-phi
	fld	QWORD PTR [ebx+8]	;; Two-to-minus-phi group multiplier
	fld	QWORD PTR [ebx]		;; Two-to-phi group multiplier
	fld	QWORD PTR [ebx+8]	;; Two-to-minus-phi
	IFDEF PFETCH
	prefetcht0 [ebx+16]		;; Prefetch next group multipliers
	ENDIF
	mov	norm_ptr1, ebx		;; Save group multiplier pointer
	mov	ebx, norm_col_mults	;; Restart the column multipliers
ilp2:	norm_macro			;; Normalize 8 values
	lea	esi, [esi+dist8]	;; Next source
	lea	ebx, [ebx+8*NMD]	;; Next set of 8 column multipliers
	IFDEF PFETCH
	prefetcht0 [edx]		;; Prefetch from next 4KB block
	prefetcht0 [edx+dist1]
	lea	edx, [edx+dist8]
	ENDIF
	add	ecx, 65536/16*65536	;; Decrement loop counter
	JNC_X	ilp2			;; Loop til done
	lea	esi, [esi-16*dist8+dist128];; Next source
	mov	ebx, norm_ptr1		;; Load group multiplier pointer
	lea	ebx, [ebx+NMD]		;; Next group multiplier
	fcompp				;; Pop group multipliers
	fcompp				;; Pop group multipliers
	IFDEF PFETCH
	lea	edx, [edx-16*dist8+2*dist1]
	ENDIF
	test	esi, 31			;; Test loop counter
	JNZ_X	ilp3			;; Loop til done
	sub	ecx, 32768-1		;; Decrement loop counter
	JS_X	idn			;; Loop til done
	sub	edi, scaling_ff		;; Make edi accurate again
	test	ecx, 2-1		;; Test next loop counter
	JNZ_X	ilp4			;; Loop til done
	test	ecx, 16*2-1		;; Test next loop counter
	lea	esi, [esi-8*dist128+dist1K];; Next source
	JNZ_X	ilp4			;; Loop til done
	test	ecx, 8*16*2-1		;; Test next loop counter
	lea	esi, [esi-16*dist1K+dist16K];; Next source
	JNZ_X	ilp4			;; Loop til done
	sub	edi, scaling_ff2	;; Make edi accurate again
	test	ecx, 8*8*16*2-1		;; Test next loop counter
	lea	esi, [esi-8*dist16K+dist128K];; Next source
	JNZ_X	ilp4			;; Loop til done
	lea	esi, [esi-8*dist128K+dist1M];; Next source
	JMP_X	ilp4			;; Loop til done
ENDM

; Macro to loop through all the FFT values and apply the proper normalization
; routine.  Used whenever we are using a rational-base FFT.

rnorm	MACRO	lab, norm_macro
	LOCAL	rlp, rlp4, pf1
	PUBLICP	lab
	LABELP	lab
	mov	esi, _DESTARG		;; Address of squared number
	mov	edi, _ADDIN_OFFSET	;; Get address to add value into
	fld	QWORD PTR [esi][edi]	;; Get the value
	fadd	_ADDIN_VALUE		;; Add in the requested value
	fstp	QWORD PTR [esi][edi]	;; Save the new value
	fld	ttmp_ff_inv		;; Preload 2 / FFTLEN
	fldz				;; Init SUMOUT
	fsub	_ADDIN_VALUE		;; Do not include addin in sumout
	fld	BIGVAL			;; Load initial carry
	mov	edx, _FFTZERO		;; Count of words to NOT zero
	mov	ecx, normcount1		;; Load loop counter
rlp4:
	IFDEF PFETCH
	lea	ebx, [esi+4*dist128]	;; Assume next block is 4*dist128 away
	test	cl, 2-1			;; Test loop counter
	jz	short rlp		;; No TLB loads needed
	mov	al, 16*2-1		;; Mask for loop counter
	and	al, cl
	lea	ebx, [esi+dist1K-4*dist128];; Assume next block is dist1K away
	cmp	al, 16*2-1		;; Test loop counter again
	jne	short pf1		;; Skip if assumption correct
	lea	ebx, [esi+dist16K-15*dist1K-4*dist128];; Next is dist16K away
pf1:	mov	al, [ebx+0*dist8]	;; Read from this page (loads the TLB)
	mov	al, [ebx+1*dist8]
	mov	al, [ebx+2*dist8]
	mov	al, [ebx+3*dist8]
	mov	al, [ebx+4*dist8]
	mov	al, [ebx+5*dist8]
	mov	al, [ebx+6*dist8]
	mov	al, [ebx+7*dist8]
	mov	al, [ebx+8*dist8]
	mov	al, [ebx+9*dist8]
	mov	al, [ebx+10*dist8]
	mov	al, [ebx+11*dist8]
	mov	al, [ebx+12*dist8]
	mov	al, [ebx+13*dist8]
	mov	al, [ebx+14*dist8]
	mov	al, [ebx+15*dist8]
	ENDIF
rlp:	norm_macro			;; Normalize 8 values
	lea	esi, [esi+dist8]	;; Next source
	IFDEF PFETCH
	prefetcht0 [ebx]		;; Prefetch from next 4KB block
	prefetcht0 [ebx+dist1]
	lea	ebx, [ebx+dist8]
	ENDIF
	add	ecx, 65536/16*65536	;; Decrement loop counter
	JNC_X	rlp			;; Loop til done
	lea	esi, [esi-16*dist8+dist128];; Next source
	IFDEF PFETCH
	lea	ebx, [ebx-16*dist8+2*dist1]
	ENDIF
	test	esi, 31			;; Test loop counter
	JNZ_X	rlp			;; Loop til done
	sub	ecx, 32768-1		;; Decrement loop counter
	JS_X	rdn			;; Loop til done
	test	ecx, 2-1		;; Test next loop counter
	JNZ_X	rlp4			;; Loop til done
	test	ecx, 16*2-1		;; Test next loop counter
	lea	esi, [esi-8*dist128+dist1K];; Next source
	JNZ_X	rlp4			;; Loop til done
	test	ecx, 8*16*2-1		;; Test next loop counter
	lea	esi, [esi-16*dist1K+dist16K];; Next source
	JNZ_X	rlp4			;; Loop til done
	test	ecx, 8*8*16*2-1		;; Test next loop counter
	lea	esi, [esi-8*dist16K+dist128K];; Next source
	JNZ_X	rlp4			;; Loop til done
	lea	esi, [esi-8*dist128K+dist1M];; Next source
	JMP_X	rlp4			;; Loop til done
ENDM

; The 16 different normalization routines

	rnorm	r4, norm_0d
	rnorm	r4e, norm_0d_e
	rnorm	r4c, norm_0d_c
	rnorm	r4ec, norm_0d_e_c
	rnorm	r4z, norm_0d_z
	rnorm	r4ze, norm_0d_z_e
	inorm	i4, norm_2d
	inorm	i4e, norm_2d_e
	inorm	i4c, norm_2d_c
	inorm	i4ec, norm_2d_e_c
	inorm	i4z, norm_2d_z
	inorm	i4ze, norm_2d_z_e

; Finish off the normalization process by add any carry to first values.
; Handle both the with and without two-to-phi array cases.

rdn:	mov	esi, _DESTARG		; Address of squared number
	norm012_0d
	JMP_X	cmnend

idn:	mov	esi, _DESTARG		; Address of squared number
	mov	edi, scaled_numlit	; To compute big vs little words
	mov	ebx, norm_col_mults	; Restart the column multipliers
	norm012_2d
end_timer 30

; Clear needs-normalize counter

cmnend:	mov	DWORD PTR [esi-4], 0

; Normalize SUMOUT value by multiplying by 1 / (fftlen/2).

	fmul	ttmp_ff_inv
	fst	QWORD PTR [esi-24]	; Save sum of FFT outputs

; Test if the sum of the output values is an error (such as infinity or NaN)

	fxam				; Test the sum of FFT outputs
	fnstsw	ax
	and	eax, 0100h		; Isolate the C0 bit (nan or infinity)
	jz	short noerr1		; If zero, no error
	or	_GWERROR, 1		; Set error flag
	fcomp	st(0)			; Pop the bad value
	jmp	short exit		; Skip second error check

; Check that the sum of the input numbers squared is approximately
; equal to the sum of inverse fft results.

noerr1:	fsub	QWORD PTR [esi-16]	; Compare to product of sum of inputs
	fabs
	fcomp	_MAXDIFF		; Compare diff to maximum allowable
	fnstsw	ax
	and	eax, 4100h		; Isolate the C3 and C0 bits
	jnz	short exit		; If non-zero, no error
	or	_GWERROR, 2		; Set error flag

; Return

exit:	pop	esi
	pop	edi
	pop	ebx
	pop	ebp
	ret
ENDPP	_comp8norm4

_TEXT32	ENDS
END
