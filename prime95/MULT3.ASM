; Copyright 1995-2000 Just For Fun Software, Inc., all rights reserved
; Author:  George Woltman
; Email: woltman@alum.mit.edu
;
; This routine implements a discrete weighted transform to quickly multiply
; two numbers.
;
; This code handles FFTs that use the convoluted memory model.
; FFT sizes between than 10K and 64K doubles are supported.
; This code does two passes, eight levels on the second pass.
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
INCLUDE	lucas1.mac
INCLUDE pfa.mac
INCLUDE mult.mac
INCLUDE fft2b.mac
INCLUDE memory.mac

IFDEF PPRO
INCLUDE	lucas1p.mac
ENDIF

EXTRNP	pass2_eight_levels_type_123
EXTRNP	pass2_eight_levels_type_123p
EXTRNP	pass2_eight_levels_type_4
EXTRNP	pass2_eight_levels_type_4p

	convoluted_distances

;; All the FFT routines for each FFT length

_gw_ffts PROC NEAR
	EXPANDING = 2
	fft	10K
	fft	12K
	fft	14K
	fft	16K
	fft	16K p
	fft	20K
	fft	24K
	fft	28K
	fft	32K
	fft	32K p
	fft	40K
	fft	48K
	fft	56K
	fft	64K
	fft	64K p
_gw_ffts ENDP

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

IFNDEF PPRO
	PUBLIC	gwaddq3
gwaddq3	PROC NEAR
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
gwaddq3	ENDP


;;
;; Add two numbers.  This will be done without carry propogation
;; (for extra speed) if it can be done safely.
;;

	PUBLIC	gwadd3
gwadd3	PROC NEAR

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
	mov	eax, addcount1		; V - Load loop counter
uaddlp:	fld	QWORD PTR [ecx]		; Load first number
	fadd	QWORD PTR [ecx][edx]	; Add in second number
	fld	QWORD PTR [ecx+8]	; Load first number
	fadd	QWORD PTR [ecx+8][edx]	; Add in second number
	fld	QWORD PTR [ecx+16]	; Load first number
	fadd	QWORD PTR [ecx+16][edx] ; Add in second number
	fxch	st(1)
	fld	QWORD PTR [ecx+24]	; Load first number
	fadd	QWORD PTR [ecx+24][edx] ; Add in second number
	fxch	st(3)
	fstp	QWORD PTR [esi]		; Save result
	fstp	QWORD PTR [esi+8]	; Save result
	fstp	QWORD PTR [esi+16]	; Save result
	fstp	QWORD PTR [esi+24]	; Save result
	lea	ecx, [ecx+dist1]	; U - Bump source pointer
	lea	esi, [esi+dist1]	; V - Bump dest pointer
	add	al, 256/8		; U - Check inner loop counter
	jnc	short uaddlp		; V - Loop if necessary
	lea	ecx, [ecx-8*dist1+dist8]; U - Bump source pointer
	lea	esi, [esi-8*dist1+dist8]; V - Bump dest pointer
	dec	al			; U - Check middle loop counter
	jnz	short uaddlp		; V - Loop if necessary
	sub	eax, 65536-16		; U - Test next loop counter
	js	short uadddn		; V - Loop til done
	lea	ecx, [ecx-16*dist8+4*dist128]; U - Bump source pointer
	lea	esi, [esi-16*dist8+4*dist128]; V - Bump dest pointer
	add	ah, 256/16		; U - Check middle loop counter
	jnc	short uaddlp		; V - Loop if necessary
	lea	ecx, [ecx-64*dist128+dist8192]; U - Skip wasted cache lines
	lea	esi, [esi-64*dist128+dist8192]; V - Skip wasted cache lines
	jmp	short uaddlp		; V - Loop if necessary
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
	add	bl, 256/16		; U - Decrement loop counter
	JNC_X	naddlp			; V - Loop til done
	lea	ecx, [ecx-16*dist8+dist128]; U - Next source
	lea	esi, [esi-16*dist8+dist128]; V - Next dest
	mov	edi, normgrpptr		; U - Restore group multiplier address
	dec	bl			; V - Decrement loop counter
	fcompp				; UV - Pop group multipliers
	fcompp				; UV - Pop group multipliers
	JNZ_X	naddlp2			; V - Loop til done
	sub	ebp, scaling_ff		;*UU - Make ebp accurate again
	sub	ebx, 65536-4		; U - Test next loop counter
	js	short nadddn		; V - Loop til done
	add	bh, 256/16		; U - Dec next loop counter
	JNC_X	naddlp2			; V - Loop til done
	lea	ecx, [ecx-64*dist128+dist8192]; V - Next source
	lea	esi, [esi-64*dist128+dist8192]; V - Next dest
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
gwadd3	ENDP


;;
;; Subtract two numbers without carry propogation.  Caller can use this for
;; consecutive add or subtract operations.  However, the last operation
;; before a multiply must use the routine that will normalize data.
;;

	PUBLIC	gwsubq3
gwsubq3	PROC NEAR
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
gwsubq3	ENDP

;;
;; Subtract two numbers.  This will be done without carry propogation
;; (for extra speed) if it can be done safely.
;;

	PUBLIC	gwsub3
gwsub3	PROC NEAR

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
	mov	eax, addcount1		; V - Load loop counter
usublp:	fld	QWORD PTR [ecx][edx]	; Load second number
	fsub	QWORD PTR [ecx]		; Subtract first number
	fld	QWORD PTR [ecx+8][edx]	; Load second number
	fsub	QWORD PTR [ecx+8]	; Subtract first number
	fld	QWORD PTR [ecx+16][edx] ; Load second number
	fsub	QWORD PTR [ecx+16]	; Subtract first number
	fxch	st(1)
	fld	QWORD PTR [ecx+24][edx] ; Load second number
	fsub	QWORD PTR [ecx+24]	; Subtract first number
	fxch	st(3)
	fstp	QWORD PTR [esi]		; Save result
	fstp	QWORD PTR [esi+8]	; Save result
	fstp	QWORD PTR [esi+16]	; Save result
	fstp	QWORD PTR [esi+24]	; Save result
	lea	ecx, [ecx+dist1]	; U - Bump source pointer
	lea	esi, [esi+dist1]	; V - Bump dest pointer
	add	al, 256/8		; U - Check inner loop counter
	jnc	short usublp		; V - Loop if necessary
	lea	ecx, [ecx-8*dist1+dist8]; U - Bump source pointer
	lea	esi, [esi-8*dist1+dist8]; V - Bump dest pointer
	dec	al			; U - Check middle loop counter
	jnz	short usublp		; V - Loop if necessary
	sub	eax, 65536-16		; U - Test next loop counter
	js	short usubdn		; V - Loop til done
	lea	ecx, [ecx-16*dist8+4*dist128]; U - Bump source pointer
	lea	esi, [esi-16*dist8+4*dist128]; V - Bump dest pointer
	add	ah, 256/16		; U - Check middle loop counter
	jnc	short usublp		; V - Loop if necessary
	lea	ecx, [ecx-64*dist128+dist8192]; U - Skip wasted cache lines
	lea	esi, [esi-64*dist128+dist8192]; V - Skip wasted cache lines
	jmp	short usublp		; V - Loop if necessary
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
	add	bl, 256/16		; U - Decrement loop counter
	JNC_X	nsublp			; V - Loop til done
	lea	ecx, [ecx-16*dist8+dist128]; U - Next source
	lea	esi, [esi-16*dist8+dist128]; V - Next dest
	mov	edi, normgrpptr		; U - Restore group multiplier address
	dec	bl			; V - Decrement loop counter
	fcompp				; UV - Pop group multipliers
	fcompp				; UV - Pop group multipliers
	JNZ_X	nsublp2			; V - Loop til done
	sub	ebp, scaling_ff		;*UU - Make ebp accurate again
	sub	ebx, 65536-4		; U - Test next loop counter
	js	short nsubdn		; V - Loop til done
	add	bh, 256/16		; U - Dec next loop counter
	JNC_X	nsublp2			; V - Loop til done
	lea	ecx, [ecx-64*dist128+dist8192]; V - Next source
	lea	esi, [esi-64*dist128+dist8192]; V - Next dest
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
gwsub3	ENDP

;;
;; Add and subtract two numbers without carry propogation.
;;

	PUBLIC	gwaddsubq3
gwaddsubq3 PROC NEAR
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
gwaddsubq3 ENDP


;;
;; Add and subtract two numbers.  This will be done without carry propogation
;; (for extra speed) if it can be done safely.
;;

	PUBLIC	gwaddsub3
gwaddsub3 PROC NEAR

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
	fld	QWORD PTR [ecx+8]	; Load first number
	fld	st(0)			; Dup first number
	fadd	QWORD PTR [ecx+8][edx]	; Add in second number
	fxch	st(1)			; S1,A1,S0,A0
	fsub	QWORD PTR [ecx+8][edx]	; Subtract out second number
	fld	QWORD PTR [ecx+16]	; Load first number
	fld	st(0)			; Dup first number
	fadd	QWORD PTR [ecx+16][edx]	; Add in second number
	fxch	st(1)			; S2,A2,S1,A1,S0,A0
	fsub	QWORD PTR [ecx+16][edx]	; Subtract out second number
	fld	QWORD PTR [ecx+24]	; Load first number
	fld	st(0)			; Dup first number
	fadd	QWORD PTR [ecx+24][edx]	; Add in second number
	fxch	st(7)			; A0,S3,S2,A2,S1,A1,S0,A3
	fstp	QWORD PTR [esi]		; Save result
	fsub	QWORD PTR [ecx+24][edx]	; Subtract out second number
	fxch	st(5)			; S0,S2,A2,S1,A1,S3,A3
	fstp	QWORD PTR [esi][ebp]	; Save result
	fstp	QWORD PTR [esi+16][ebp]	; Save result
	fstp	QWORD PTR [esi+16]	; Save result
	fstp	QWORD PTR [esi+8][ebp]	; Save result
	fstp	QWORD PTR [esi+8]	; Save result
	fstp	QWORD PTR [esi+24][ebp]	; Save result
	fstp	QWORD PTR [esi+24]	; Save result
	lea	ecx, [ecx+dist1]	; U - Bump source pointer
	lea	esi, [esi+dist1]	; V - Bump dest pointer
	add	al, 256/8		; U - Check inner loop counter
	jnc	short uaddsublp		; V - Loop if necessary
	lea	ecx, [ecx-8*dist1+dist8]; U - Bump source pointer
	lea	esi, [esi-8*dist1+dist8]; V - Bump dest pointer
	dec	al			; U - Check middle loop counter
	jnz	short uaddsublp		; V - Loop if necessary
	sub	eax, 65536-16		; U - Test next loop counter
	js	short uaddsubdn		; V - Loop til done
	lea	ecx, [ecx-16*dist8+4*dist128]; U - Bump source pointer
	lea	esi, [esi-16*dist8+4*dist128]; V - Bump dest pointer
	add	ah, 256/16		; U - Check middle loop counter
	JNC_X	uaddsublp		; V - Loop if necessary
	lea	ecx, [ecx-64*dist128+dist8192]; U - Skip wasted cache lines
	lea	esi, [esi-64*dist128+dist8192]; V - Skip wasted cache lines
	JMP_X	uaddsublp		; V - Loop if necessary
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
	add	dl, 256/16		; U - Decrement loop counter
	JNC_X	naddsublp		; V - Loop til done
	lea	ecx, [ecx-16*dist8+dist128]; U - Next source
	lea	esi, [esi-16*dist8+dist128]; V - Next dest
	mov	edi, normgrpptr		; U - Restore group multiplier address
	dec	dl			; V - Decrement loop counter
	fcompp				; UV - Pop group multipliers
	fcompp				; UV - Pop group multipliers
	JNZ_X	naddsublp2		; V - Loop til done
	sub	ebx, scaling_ff		;*UU - Make ebx accurate again
	sub	edx, 65536-4		; U - Test next loop counter
	js	short naddsubdn		; V - Loop til done
	add	dh, 256/16		; U - Dec next loop counter
	JNC_X	naddsublp2		; V - Loop til done
	lea	ecx, [ecx-64*dist128+dist8192]; V - Next source
	lea	esi, [esi-64*dist128+dist8192]; V - Next dest
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
gwaddsub3 ENDP
ENDIF

;; Routines to do the normalization after a multiply

_comp8norm PROC NEAR

; Normalize the results and multiply them by 3.  Used in P-1 factoring.
; As a safety measure, make sure the top bits of ERRCHK is AB93.  This
; makes it highly unlikely that ERRCHK will be corrupted and LL tests
; will all of a sudden start doing a multiply by 3.

nt3:	mov	eax, _ERRCHK		; Load variable
	and	eax, 0FFFF0000h		; Mask out bits
	cmp	eax, 0AB930000h		; See if the upper bits are correct
	JNE_X	noechk			; This should never happen
	fldz				; Init SUMOUT
	fld	BIGVAL			; Load initial carry
	mov	edx, norm_grp_mults	; Address of the group multipliers
	mov	esi, _DESTARG		; Address of squared number
	mov	edi, scaled_numlit	; Used to compute big vs little words
	mov	ebp, scaled_numbig
	mov	eax, -1			; First word is a big word (-1)
	mov	ecx, normcount1		; Load loop counter
plp3:	fld	QWORD PTR [edx]		; Two-to-phi
	fld	QWORD PTR [edx+8]	; Two-to-minus-phi group multiplier
	fld	QWORD PTR [edx]		; Two-to-phi group multiplier
	fld	QWORD PTR [edx+8]	; Two-to-minus-phi
	mov	ebx, norm_col_mults	; Restart the column multipliers
plp2:	normalize_times_3		; Normalize 8 values
	lea	esi, [esi+dist8]	; Next source
	lea	ebx, [ebx+8*NMD]	; Next set of 8 column multipliers
	add	cl, 256/16		; U - Decrement loop counter
	JNC_X	plp2			; V - Loop til done
	lea	esi, [esi-16*dist8+dist128]; U - Next source
	lea	edx, [edx+NMD]		; V - Next group multiplier
	fcompp				; UV - Pop group multipliers
	fcompp				; UV - Pop group multipliers
	dec	cl			; U - Decrement loop counter
	JNZ_X	plp3			; V - Loop til done
	sub	edi, scaling_ff		;*UU - Make ebp accurate again
	sub	ecx, 65536-4		; U - Test counter / restore inner
	JS_X	ndn			; V - Loop til done
	add	ch, 256/16		; U - Test next loop counter
	JNC_X	plp3			; V - Loop til done
	lea	esi, [esi-64*dist128+dist8192]; U - Next source
	JMP_X	plp3			; V - Loop til done

; Error check - find largest round off error

echk:	fld	_MAXERR			; Maximum error
	mov	edx, norm_grp_mults	; Address of the group multipliers
	mov	esi, _DESTARG		; Address of squared number
	mov	edi, scaled_numlit	; Used to compute big vs little words
	mov	ebp, scaled_numbig
	mov	eax, -1			; First word is a big word (-1)
	mov	ecx, normcount1		; N iterations of
elp3:	fld	QWORD PTR [edx+8]	; Two-to-minus-phi group multiplier
	fld	QWORD PTR [edx+8]	; Two-to-minus-phi
	mov	ebx, norm_col_mults	; Restart the column multipliers
elp2:	error_check			; Error check eight values
	lea	esi, [esi+dist8]	; Next source
	lea	ebx, [ebx+8*NMD]	; Next set of 8 column multipliers
	add	cl, 256/16		; U - Decrement loop counter
	JNC_X	elp2			; V - Loop til done
	lea	esi, [esi-16*dist8+dist128]; U - Next source
	lea	edx, [edx+NMD]		; V - Next group multiplier
	fcompp				; UV - Pop group multipliers
	dec	cl			; U - Decrement loop counter
	JNZ_X	elp3			; V - Loop til done
	sub	edi, scaling_ff		;*UU - Make ebp accurate again
	sub	ecx, 65536-4		; U - Test counter / restore inner
	js	short edn		; V - Loop til done
	add	ch, 256/16		; U - Test next loop counter
	JNC_X	elp3			; V - Loop til done
	lea	esi, [esi-64*dist128+dist8192]; U - Next source
	JMP_X	elp3			; V - Loop til done
edn:	fstp	_MAXERR			; Store maximum error

; Normalize the results

noechk:	fldz				; Init SUMOUT
	fld	BIGVAL			; Load initial carry
	mov	edx, norm_grp_mults	; Address of the group multipliers
	mov	esi, _DESTARG		; Address of squared number
	mov	edi, scaled_numlit	; Used to compute big vs little words
	mov	ebp, scaled_numbig
	mov	eax, -1			; First word is a big word (-1)
	mov	ecx, normcount1		; Load loop counter
nlp3:	fld	QWORD PTR [edx]		; Two-to-phi
	fld	QWORD PTR [edx+8]	; Two-to-minus-phi group multiplier
	fld	QWORD PTR [edx]		; Two-to-phi group multiplier
	fld	QWORD PTR [edx+8]	; Two-to-minus-phi
	mov	ebx, norm_col_mults	; Restart the column multipliers
nlp2:	normalize			; Normalize 8 values
	lea	esi, [esi+dist8]	; Next source
	lea	ebx, [ebx+8*NMD]	; Next set of 8 column multipliers
	add	cl, 256/16		; U - Decrement loop counter
	JNC_X	nlp2			; V - Loop til done
	lea	esi, [esi-16*dist8+dist128]; U - Next source
	lea	edx, [edx+NMD]		; V - Next group multiplier
	fcompp				; UV - Pop group multipliers
	fcompp				; UV - Pop group multipliers
	dec	cl			; U - Decrement loop counter
	JNZ_X	nlp3			; V - Loop til done
	sub	edi, scaling_ff		;*UU - Make ebp accurate again
	sub	ecx, 65536-4		; U - Test counter / restore inner
	js	short ndn		; V - Loop til done
	add	ch, 256/16		; U - Test next loop counter
	JNC_X	nlp3			; V - Loop til done
	lea	esi, [esi-64*dist128+dist8192]; U - Next source
	JMP_X	nlp3			; V - Loop til done
					; Add any carry to first values
ndn:	mov	esi, _DESTARG		; Address of squared number
	mov	edi, scaled_numlit	; To compute big vs little words
	mov	ebx, norm_col_mults	; Restart the column multipliers
	cmp	_PLUS1, 0
	jz	short nlp1
	fsub	BIGVAL
	fchs
	fadd	BIGVAL
nlp1:	normalize_012
	mov	DWORD PTR [esi-4], 0	; Clear needs-normalize counter

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
_comp8norm ENDP

_TEXT32	ENDS
END
