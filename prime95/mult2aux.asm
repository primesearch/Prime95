; Copyright 1995-2000 Just For Fun Software, Inc., all rights reserved
; Author:  George Woltman
; Email: woltman@alum.mit.edu
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
INCLUDE normal.mac

	flat_distances

;;
;; Add two numbers without carry propogation.  Caller can use this for
;; consecutive add or subtract operations.  However, the last operation
;; before a multiply must use the routine that will normalize data.
;;

IFNDEF PPRO
	PUBLIC	gwaddq2
gwaddq2	PROC NEAR
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
gwaddq2	ENDP


;;
;; Add two numbers.  This will be done without carry propogation
;; (for extra speed) if it can be done safely.
;;

	PUBLIC	gwadd2
gwadd2	PROC NEAR

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
	jg	short nadd		; V - yes

; Do an unnormalized add
; ecx = src arg #1
; edx = distance to src arg #2
; esi = dest arg
; eax = new needs-normalize counter

uadd:	mov	[esi-4], eax		; U - Store needs-normalize counter
	mov	eax, addcount1		; V - Load loop counter
uaddlp:	fld	QWORD PTR [ecx]		; Load first number
	fadd	QWORD PTR [ecx][edx]	; Add in second number
	fld	QWORD PTR [ecx+dist1]	; Load first number
	fadd	QWORD PTR [ecx+dist1][edx] ; Add in second number
	fld	QWORD PTR [ecx+2*dist1]	; Load first number
	fadd	QWORD PTR [ecx+2*dist1][edx] ; Add in second number
	fxch	st(1)
	fld	QWORD PTR [ecx+3*dist1]	; Load first number
	fadd	QWORD PTR [ecx+3*dist1][edx] ; Add in second number
	fxch	st(3)
	fstp	QWORD PTR [esi]		; Save result
	fstp	QWORD PTR [esi+dist1]	; Save result
	fstp	QWORD PTR [esi+2*dist1]	; Save result
	fstp	QWORD PTR [esi+3*dist1]	; Save result
	lea	ecx, [ecx+4*dist1]	; U - Bump source pointer
	lea	esi, [esi+4*dist1]	; V - Bump dest pointer
	dec	al			; U - Check inner loop counter
	jnz	short uaddlp		; V - Loop if necessary
	lea	ecx, [ecx-512*dist1+dist512] ; U - Bump source pointer
	lea	esi, [esi-512*dist1+dist512] ; V - Bump dest pointer
	sub	eax, 256-128		; U - Check outer/restore inner count
	jns	short uaddlp		; V - Loop if necessary
	mov	eax, addcount2		; U - Load a new loop counter
	jc	short uaddlp		; V - Do the second set of counts
	pop	edi
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
naddlp:	norm_op_2d fadd, 8*dist1	; Add and normalize 8 values
	dec	bl			; U - Decrement loop counter
	JNZ_X	naddlp			; V - Loop til done
	mov	edi, normgrpptr		; U - Restore group multiplier address
	dec	bh			; V - Decrement loop counter
	fcompp				; UV - Pop group multipliers
	fcompp				; UV - Pop group multipliers
	mov	bl, 4			; U - Reload inner loop counter
	JNZ_X	naddlp2			; V - Loop til done
	sub	ebp, scaling_ff		;*UU - Make ebp accurate again
	lea	ecx, [ecx-512*dist1+dist512]; V - Next source
	lea	esi, [esi-512*dist1+dist512]; V - Next dest
	sub	ebx, 65536-16*256	; U - Dec outer/restore mid loop count
	JNS_X	naddlp2			; V - Loop til done
	mov	ebx, normcount2		; U - Load a new loop counter
	JC_X	naddlp2			; V - Do a second set of counts
	fsub	BIGVAL			; Compute carry
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
gwadd2	ENDP


;;
;; Subtract two numbers without carry propogation.  Caller can use this for
;; consecutive add or subtract operations.  However, the last operation
;; before a multiply must use the routine that will normalize data.
;;

	PUBLIC	gwsubq2
gwsubq2	PROC NEAR
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
gwsubq2	ENDP

;;
;; Subtract two numbers.  This will be done without carry propogation
;; (for extra speed) if it can be done safely.
;;

	PUBLIC	gwsub2
gwsub2	PROC NEAR

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
	jg	short nsub		; V - yes

; Do an unnormalized subtract
; ecx = src arg #1
; edx = distance to src arg #2
; esi = dest arg
; eax = new needs-normalize counter

usub:	mov	[esi-4], eax		; U - Store needs-normalize counter
	mov	eax, addcount1		; V - Load loop counter
usublp:	fld	QWORD PTR [ecx][edx]	; Load second number
	fsub	QWORD PTR [ecx]		; Subtract first number
	fld	QWORD PTR [ecx+dist1][edx] ; Load second number
	fsub	QWORD PTR [ecx+dist1]	; Subtract first number
	fld	QWORD PTR [ecx+2*dist1][edx] ; Load second number
	fsub	QWORD PTR [ecx+2*dist1]	; Subtract first number
	fxch	st(1)
	fld	QWORD PTR [ecx+3*dist1][edx] ; Load second number
	fsub	QWORD PTR [ecx+3*dist1]	; Subtract first number
	fxch	st(3)
	fstp	QWORD PTR [esi]		; Save result
	fstp	QWORD PTR [esi+dist1]	; Save result
	fstp	QWORD PTR [esi+2*dist1]	; Save result
	fstp	QWORD PTR [esi+3*dist1]	; Save result
	lea	ecx, [ecx+4*dist1]	; U - Bump source pointer
	lea	esi, [esi+4*dist1]	; V - Bump dest pointer
	dec	al			; U - Check inner loop counter
	jnz	short usublp		; V - Loop if necessary
	lea	ecx, [ecx-512*dist1+dist512] ; U - Bump source pointer
	lea	esi, [esi-512*dist1+dist512] ; V - Bump dest pointer
	sub	eax, 256-128		; U - Check outer/restore inner count
	jns	short usublp		; V - Loop if necessary
	mov	eax, addcount2		; U - Load a new loop counter
	jc	short usublp		; V - Do the second set of counts
	pop	edi
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
nsublp:	norm_op_2d fsub, 8*dist1	; Subtract and normalize 8 values
	dec	bl			; U - Decrement loop counter
	JNZ_X	nsublp			; V - Loop til done
	mov	edi, normgrpptr		; U - Restore group multiplier address
	dec	bh			; V - Decrement loop counter
	fcompp				; UV - Pop group multipliers
	fcompp				; UV - Pop group multipliers
	mov	bl, 4			; U - Reload inner loop counter
	JNZ_X	nsublp2			; V - Loop til done
	sub	ebp, scaling_ff		;*UU - Make ebp accurate again
	lea	ecx, [ecx-512*dist1+dist512]; V - Next source
	lea	esi, [esi-512*dist1+dist512]; V - Next dest
	sub	ebx, 65536-16*256	; U - Dec outer/restore mid loop count
	JNS_X	nsublp2			; V - Loop til done
	mov	ebx, normcount2		; U - Load a new loop counter
	JC_X	nsublp2			; V - Do a second set of counts
	fsub	BIGVAL			; Compute carry
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
gwsub2	ENDP

;;
;; Add and subtract two numbers without carry propogation.
;;

	PUBLIC	gwaddsubq2
gwaddsubq2 PROC NEAR
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
gwaddsubq2 ENDP

;;
;; Add and subtract two numbers.  This will be done without carry propogation
;; (for extra speed) if it can be done safely.
;;

	PUBLIC	gwaddsub2
gwaddsub2 PROC NEAR

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
	mov	eax, addcount1		; U - Load loop counter
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
	lea	ecx, [ecx+4*dist1]	; U - Bump source pointer
	lea	esi, [esi+4*dist1]	; V - Bump dest pointer
	dec	al			; U - Check inner loop counter
	jnz	short uaddsublp		; V - Loop if necessary
	lea	ecx, [ecx-512*dist1+dist512] ; U - Bump source pointer
	lea	esi, [esi-512*dist1+dist512] ; V - Bump dest pointer
	sub	eax, 256-128		; U - Check outer/restore inner count
	jns	short uaddsublp		; V - Loop if necessary
	mov	eax, addcount2		; U - Load a new loop counter
	jc	short uaddsublp		; V - Do the second set of counts
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
	norm_addsub_2d 8*dist1		; Add and normalize 8 values
	dec	dl			; U - Decrement loop counter
	JNZ_X	naddsublp		; V - Loop til done
	mov	edi, normgrpptr		; U - Restore group multiplier address
	dec	dh			; V - Decrement loop counter
	fcompp				; UV - Pop group multipliers
	fcompp				; UV - Pop group multipliers
	mov	dl, 4			; U - Reload inner loop counter
	JNZ_X	naddsublp2		; V - Loop til done
	sub	ebx, scaling_ff		;*UU - Make ebx accurate again
	lea	ecx, [ecx-512*dist1+dist512]; V - Next source
	lea	esi, [esi-512*dist1+dist512]; V - Next dest
	sub	edx, 65536-16*256	; U - Dec outer/restore mid loop count
	JNS_X	naddsublp2		; V - Loop til done
	mov	edx, normcount2		; U - Load a new loop counter
	JC_X	naddsublp2		; V - Do a second set of counts
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
gwaddsub2 ENDP

;;
;; Copy one number and zero some low order words.
;;

	PUBLIC	gwcopyzero2
gwcopyzero2 PROC NEAR
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
	lea	esi, [esi+8*dist1]	; U - Bump source pointer
	lea	edi, [edi+8*dist1]	; V - Bump dest pointer
	dec	al			; U - Decrement loop counter
	JNZ_X	zlp			; V - Loop til done
	dec	ah			; U - Decrement loop counter
	mov	al, 4			; V - Reload inner loop counter
	JNZ_X	zlp			; V - Loop til done
	lea	esi, [esi-512*dist1+dist512]; U - Next source
	lea	edi, [edi-512*dist1+dist512]; V - Next dest
	sub	eax, 65536-16*256	; U - Dec outer/restore mid loop count
	JNS_X	zlp			; V - Loop til done
	mov	eax, normcount2		; U - Load a new loop counter
	JC_X	zlp			; V - Do a second set of counts
	pop	edi
	pop	esi
	ret
gwcopyzero2 ENDP

;;
;; Do a mod k*2^n+/-1
;;

	PUBLIC	gwprothmod2
gwprothmod2 PROC NEAR
	push	esi			; U - Save registers
	push	edi			; V - Save registers
	push	ebx
	mov	esi, _SRCARG		; U - Address of first number
	prothmod_upper_prep_0d		; Prepare for prothmod
	mov	ebx, normcount2		; U - Load a new loop counter
pmulp:	prothmod_upper_0d pmudn		; Divide upper values by k
	lea	esi, [esi-8*dist1]	; U - Next source
	lea	edi, [edi+8*8]		; V - Next scratch words
	dec	bl			; U - Decrement loop counter
	JNZ_X	pmulp			; V - Loop til done
	dec	bh			; V - Decrement loop counter
	mov	bl, 4			; U - Reload inner loop counter
	JNZ_X	pmulp			; V - Loop til done
	lea	esi, [esi+512*dist1-dist512]; V - Next dest
	sub	ebx, 65536-16*256	; U - Dec outer/restore mid loop count
	JNS_X	pmulp			; V - Loop til done
	mov	ebx, normcount1		; U - Load a new loop counter
	JMP_X	pmulp			; V - Do a second set of counts
pmudn:	prothmod_lower_prep_0d		; Prepare for prothmod
	mov	ebx, normcount1		; U - Load a new loop counter
pmllp:	prothmod_lower_0d		; Divide upper values by k
	dec	edx
	jz	short pmldn
	lea	esi, [esi+8*dist1]	; U - Next source ptr
	lea	edi, [edi-8*8]		; V - Next scratch words
	dec	bl			; U - Decrement loop counter
	JNZ_X	pmllp			; V - Loop til done
	dec	bh			; V - Decrement loop counter
	mov	bl, 4			; U - Reload inner loop counter
	JNZ_X	pmllp			; V - Loop til done
	lea	esi, [esi-512*dist1+dist512]; V - Next dest
	sub	ebx, 65536-16*256	; U - Dec outer/restore mid loop count
	JNS_X	pmllp			; V - Loop til done
	mov	ebx, normcount2		; U - Load a new loop counter
	JMP_X	pmllp			; V - Do a second set of counts
pmldn:	prothmod_final_0d		; Finish off the mod
	pop	ebx
	pop	edi			; U - Restore registers
	pop	esi			; V - Restore registers
	ret
gwprothmod2 ENDP
ENDIF

;; Routines to do the normalization after a multiply

IFNDEF PPRO
_comp6norm PROC NEAR

; Macro to loop through all the FFT values and apply the proper normalization
; routine.  Used whenever we are using an irrational-base FFT.

inorm	MACRO	lab, norm_macro
	LOCAL	ilp2, ilp3
	PUBLIC	lab
lab:	mov	esi, _DESTARG		;; Address of multiplied number
	mov	edi, _ADDIN_OFFSET	;; Get address to add value into
	fld	QWORD PTR [esi][edi]	;; Get the value
	fadd	_ADDIN_VALUE		;; Add in the requested value
	fstp	QWORD PTR [esi][edi]	;; Save the new value
	fldz				;; Init SUMOUT
	fsub	_ADDIN_VALUE		;; Do not include addin in sumout
	fld	BIGVAL			;; Start process with no carry
	mov	edx, norm_grp_mults	;; Address of the group multipliers
	mov	edi, scaled_numlit	;; Used to compute big vs little words
	mov	ebp, scaled_numbig
	mov	eax, -1			;; First word is a big word (-1)
	mov	ecx, normcount1		;; Load loop counter
ilp3:	fld	QWORD PTR [edx]		;; Two-to-phi
	fld	QWORD PTR [edx+8]	;; Two-to-minus-phi group multiplier
	fld	QWORD PTR [edx]		;; Two-to-phi group multiplier
	fld	QWORD PTR [edx+8]	;; Two-to-minus-phi
	mov	ebx, norm_col_mults	;; Restart the column multipliers
ilp2:	norm_macro			;; Normalize 8 values
	lea	esi, [esi+8*dist1]	;; Next source
	lea	ebx, [ebx+8*NMD]	;; Next set of 8 column multipliers
	dec	cl			;; Test inner loop counter
	JNZ_X	ilp2			;; Loop til done
	lea	edx, [edx+NMD]		;; Next group multiplier
	dec	ch			;; Test middle loop counter
	fcompp				;; Pop group multipliers
	fcompp				;; Pop group multipliers
	mov	cl, 4			;; Reload inner load counter
	JNZ_X	ilp3			;; Loop til done
	lea	esi, [esi-512*dist1+dist512];; Next source
	mov	ebx, scaling_ff		;; Make edi accurate again
	sub	edi, ebx
	sub	ecx, 65536-16*256	;; Dec outer/restore mid loop counter
	JNS_X	ilp3			;; Loop til done
	mov	ecx, normcount2		;; Load a new loop counter
	JC_X	ilp3			;; Do a second set of counts
	JMP_X	idn
	ENDM

; Macro to loop through all the FFT values and apply the proper normalization
; routine.  Used whenever we are using a rational-base FFT.

rnorm	MACRO	lab, norm_macro
	LOCAL	rlp
	PUBLIC	lab
lab:	mov	esi, _DESTARG		;; Address of squared number
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
rlp:	norm_macro			;; Normalize 8 values
	lea	esi, [esi+8*dist1]	;; Next source
	dec	cl			;; Test inner loop counter
	JNZ_X	rlp			;; Loop til done
	dec	ch			;; Test middle loop counter
	mov	cl, 4			;; Reload inner load counter
	JNZ_X	rlp			;; Loop til done
	lea	esi, [esi-512*dist1+dist512];; Next source
	sub	ecx, 65536-16*256	;; Dec outer/restore mid loop counter
	JNS_X	rlp			;; Loop til done
	mov	ecx, normcount2		;; Load a new loop counter
	JC_X	rlp			;; Do a second set of counts
	JMP_X	rdn
	ENDM

; The 16 different normalization routines

	rnorm	r2, norm_0d
	rnorm	r2e, norm_0d_e
	rnorm	r2c, norm_0d_c
	rnorm	r2ec, norm_0d_e_c
	rnorm	r2z, norm_0d_z
	rnorm	r2ze, norm_0d_z_e
	inorm	i2, norm_2d
	inorm	i2e, norm_2d_e
	inorm	i2c, norm_2d_c
	inorm	i2ec, norm_2d_e_c
	inorm	i2z, norm_2d_z
	inorm	i2ze, norm_2d_z_e

; Finish off the normalization process by add any carry to first values.
; Handle both the with and without two-to-phi array cases.

rdn:	mov	esi, _DESTARG		; Address of squared number
	norm012_0d
	JMP_X	cmnend

idn:	mov	esi, _DESTARG		; Address of squared number
	mov	edi, scaled_numlit	; To compute big vs little words
	mov	ebx, norm_col_mults	; Restart the column multipliers
	norm012_2d

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
_comp6norm ENDP
ENDIF

_TEXT32	ENDS
END
