; Copyright 1995-2000 Just For Fun Software, Inc., all rights reserved
; Author:  George Woltman
; Email: woltman@alum.mit.edu
;
; This routine implements a discrete weighted transform to quickly multiply
; two numbers.
;
; This code handles FFTs that use the almost flat memory model and/or
; the complex normalization code.  FFT sizes from 160 to 8192 doubles
; are supported.  This code does one or two passes, if two passes then
; six levels are done on the second pass.
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
INCLUDE fft1b.mac
INCLUDE fft2a.mac
INCLUDE memory.mac

IFDEF PPRO
INCLUDE	lucas1p.mac
ENDIF

EXTRNP	pass2_six_levels_type_123
EXTRNP	pass2_six_levels_type_123p
EXTRNP	pass2_six_levels_type_4
EXTRNP	pass2_six_levels_type_4p

	flat_distances

;; All the FFT routines for each FFT length

_gw_ffts PROC NEAR
	EXPANDING = 0
	fft	160
	fft	192
;	fft	224
	fft	256
	fft	256 p
	fft	320
	fft	384
	fft	448
	fft	512
	fft	512 p
	fft	640
	fft	768
	fft	896
	fft	1024
	fft	1024 p
	EXPANDING = 2
	fft	1280
	fft	1536
	fft	1792
	fft	2048
	fft	2048 p
	fft	2560
	fft	3072
	fft	3584
	fft	4096
	fft	4096 p
	fft	5120
	fft	6144
	fft	7168
	fft	8192
	fft	8192 p
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
naddlp:	norm_op	fadd 8*dist1		; Add and normalize 8 values
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
nsublp:	norm_op	fsub 8*dist1		; Subtract and normalize 8 values
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
	norm_addsub 8*dist1		; Add and normalize 8 values
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
ENDIF

;; Routines to do the normalization after a multiply

_comp6norm PROC NEAR

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
	lea	esi, [esi+8*dist1]	; Next source
	lea	ebx, [ebx+8*NMD]	; Next set of 8 column multipliers
	dec	cl			; Test inner loop counter
	JNZ_X	plp2			; Loop til done
	lea	edx, [edx+NMD]		; Next group multiplier
	dec	ch			; Test middle loop counter
	fcompp				; Pop group multipliers
	fcompp				; Pop group multipliers
	mov	cl, 4			; Reload inner load counter
	JNZ_X	plp3			; Loop til done
	lea	esi, [esi-512*dist1+dist512]; Next source
	mov	ebx, scaling_ff		; Make edi accurate again
	sub	edi, ebx
	sub	ecx, 65536-16*256	; Dec outer/restore mid loop counter
	JNS_X	plp3			; Loop til done
	mov	ecx, normcount2		; U - Load a new loop counter
	JC_X	plp3			; V - Do a second set of counts
	JMP_X	ndn			; Go add in final carry

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
	lea	esi, [esi+8*dist1]	; Next source
	lea	ebx, [ebx+8*NMD]	; Next set of 8 column multipliers
	dec	cl			; Test inner loop counter
	JNZ_X	elp2			; Loop til done
	lea	edx, [edx+NMD]		; Next group multiplier
	dec	ch			; Test middle loop counter
	fcompp				; Pop group multipliers
	mov	cl, 4			; Reload inner load counter
	JNZ_X	elp3			; Loop til done
	lea	esi, [esi-512*dist1+dist512]; Next source
	mov	ebx, scaling_ff		; Make edi accurate again
	sub	edi, ebx
	sub	ecx, 65536-16*256	; Dec outer/restore mid loop counter
	JNS_X	elp3			; Loop til done
	mov	ecx, normcount2		; U - Load a new loop counter
	JC_X	elp3			; V - Do a second set of counts
	fstp	_MAXERR			; Store maximum error

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
	lea	esi, [esi+8*dist1]	; Next source
	lea	ebx, [ebx+8*NMD]	; Next set of 8 column multipliers
	dec	cl			; Test inner loop counter
	JNZ_X	nlp2			; Loop til done
	lea	edx, [edx+NMD]		; Next group multiplier
	dec	ch			; Test middle loop counter
	fcompp				; Pop group multipliers
	fcompp				; Pop group multipliers
	mov	cl, 4			; Reload inner load counter
	JNZ_X	nlp3			; Loop til done
	lea	esi, [esi-512*dist1+dist512]; Next source
	mov	ebx, scaling_ff		; Make edi accurate again
	sub	edi, ebx
	sub	ecx, 65536-16*256	; Dec outer/restore mid loop counter
	JNS_X	nlp3			; Loop til done
	mov	ecx, normcount2		; U - Load a new loop counter
	JC_X	nlp3			; V - Do a second set of counts

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
_comp6norm ENDP

_TEXT32	ENDS
END
