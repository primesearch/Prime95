; Copyright 1995-2000 Just For Fun Software, Inc., all rights reserved
; Author:  George Woltman
; Email: woltman@alum.mit.edu
;
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

	dist1 =	8

;;
;; Add two numbers without carry propogation.  Caller can use this for
;; consecutive add or subtract operations.  However, the last operation
;; before a multiply must use the routine that will normalize data.
;;

IFNDEF PPRO
	PUBLIC	gwaddq1
gwaddq1	PROC NEAR
	mov	ecx, _SRCARG		; U - Address of first number
	mov	edx, _SRC2ARG		; V - Address of second number
	push	esi			; U - Save registers
	push	edi			; V - Save registers
	mov	eax, [ecx-4]		; U - Load needs-normalize counter
	mov	edi, [edx-4]		; V - Load needs-normalize counter
	mov	esi, _DESTARG		; U - Address of destination
					; V - Stall
	lea	eax, [eax+edi+1]	; U - Set new needs-normalize counter
	mov	edi, _FFTLEN		; V - Load loop counter
	mov	[esi-4], eax		; U - Store needs-normalize counter
	jmp	short uaddlp		; V - Do the unnormalized add
gwaddq1	ENDP


;;
;; Add two numbers.  This will be done without carry propogation
;; (for extra speed) if it can be done safely.
;;

	PUBLIC	gwadd1
gwadd1	PROC NEAR

; See if add can be done without normalization

	mov	ecx, _SRCARG		; U - Address of first number
	mov	edx, _SRC2ARG		; V - Address of second number
	push	esi			; U - Save registers
	push	edi			; V - Save registers
	mov	eax, [ecx-4]		; U - Load needs-normalize counter
	mov	edi, [edx-4]		; V - Load needs-normalize counter
	mov	esi, _DESTARG		; U - Address of destination
					; V - Stall
	lea	eax, [eax+edi+1]	; U - Set new needs-normalize counter
	mov	edi, extra_bits		; V - Load compare value
	cmp	eax, edi		; U - Is normalization needed?
	mov	edi, _FFTLEN		; V - Load loop counter
	mov	[esi-4], eax		; U - Store needs-normalize counter
	jg	short nadd		; V - Yes, do a normalized add

; Do an unnormalized add
; ecx = src arg #1
; edx = src arg #2
; esi = dest arg
; edi = loop counter

uaddlp:	fld	QWORD PTR [ecx][edi*8-32] ; Load first number
	fadd	QWORD PTR [edx][edi*8-32] ; Add in second number
	fld	QWORD PTR [ecx][edi*8-24] ; Load first number
	fadd	QWORD PTR [edx][edi*8-24] ; Add in second number
	fld	QWORD PTR [ecx][edi*8-16] ; Load first number
	fadd	QWORD PTR [edx][edi*8-16] ; Add in second number
	fxch	st(1)
	fld	QWORD PTR [ecx][edi*8-8]  ; Load first number
	fadd	QWORD PTR [edx][edi*8-8]  ; Add in second number
	fxch	st(3)
	fstp	QWORD PTR [esi][edi*8-32] ; Save result
	fstp	QWORD PTR [esi][edi*8-24] ; Save result
	fstp	QWORD PTR [esi][edi*8-16] ; Save result
	fstp	QWORD PTR [esi][edi*8-8]  ; Save result
	sub	edi, 4			; U - Check loop counter
	jnz	short uaddlp		; V - Loop if necessary
	pop	edi			; U - Restore registers
	pop	esi			; V - Restore registers
	ret

; Do a normalized add
; ecx = src arg #1
; edx = src arg #2
; esi = dest arg
; edi = loop counter

nadd:	fninit
	push	ebp			; U - Save registers
	push	ebx			; V - Save registers
	fld	BIGVAL			; UV - Start process with no carry
	mov	ebx, norm_col_mults	; U - Address of the multipliers
	mov	ebp, scaled_numlit	; V - Computes big vs little word flag
	mov	eax, -1			; U - First word is a big word (-1)
naddlp:	norm_op_1d fadd			; Add and normalize 8 values
	sub	edi, 8			; U - Decrement loop counter
	JNZ_X	naddlp			; V - Loop til done
	fsub	BIGVAL			; UV - Get carry without BIGVAL
	mov	ecx, _DESTARG		; U - Address of squared number
	cmp	_PLUS1, 0		; UV - Are we working mod 2^N+1
	jz	short nadd1		; V - No, leave carry alone
	fchs				; UV - Yes, flip carry's sign bit
nadd1:	fadd	QWORD PTR [ecx]		; UV - Add in carry
	mov	DWORD PTR [ecx-4], 0	; UV - Clear needs-normalize flag
	pop	ebx			; U - Restore registers
	pop	ebp			; V - Restore registers
	pop	edi			; U - Restore registers
	pop	esi			; V - Restore registers
	fstp	QWORD PTR [ecx]		; UV - Save new value
	ret
gwadd1	ENDP


;;
;; Subtract two numbers without carry propogation.  Caller can use this for
;; consecutive add or subtract operations.  However, the last operation
;; before a multiply must use the routine that will normalize data.
;;

	PUBLIC	gwsubq1
gwsubq1	PROC NEAR
	mov	ecx, _SRCARG		; U - Address of first number
	mov	edx, _SRC2ARG		; V - Address of second number
	push	esi			; U - Save registers
	push	edi			; V - Save registers
	mov	eax, [ecx-4]		; U - Load needs-normalize counter
	mov	edi, [edx-4]		; V - Load needs-normalize counter
	mov	esi, _DESTARG		; U - Address of destination
					; V - Stall
	lea	eax, [eax+edi+1]	; U - Set new needs-normalize counter
	mov	edi, _FFTLEN		; V - Load loop counter
	mov	[esi-4], eax		; U - Store needs-normalize counter
	jmp	short usublp		; V - Do the unnormalized sub
gwsubq1	ENDP

;;
;; Subtract two numbers.  This will be done without carry propogation
;; (for extra speed) if it can be done safely.
;;

	PUBLIC	gwsub1
gwsub1	PROC NEAR

; See if subtract can be done without normalization

	mov	ecx, _SRCARG		; U - Address of first number
	mov	edx, _SRC2ARG		; V - Address of second number
	push	esi			; U - Save registers
	push	edi			; V - Save registers
	mov	eax, [ecx-4]		; U - Load needs-normalize counter
	mov	edi, [edx-4]		; V - Load needs-normalize counter
	mov	esi, _DESTARG		; U - Address of destination
					; V - Stall
	lea	eax, [eax+edi+1]	; U - Set new needs-normalize counter
	mov	edi, extra_bits		; V - Load compare value
	cmp	eax, edi		; U - Is normalization needed?
	mov	edi, _FFTLEN		; V - Load loop counter
	mov	[esi-4], eax		; U - Store needs-normalize counter
	jg	short nsub		; V - Yes, do a normalized subtract

; Do an unnormalized subtract
; ecx = src arg #1
; edx = src arg #2
; esi = dest arg
; edi = loop counter

usublp:	fld	QWORD PTR [edx][edi*8-32] ; Load second number
	fsub	QWORD PTR [ecx][edi*8-32] ; Subtract first number
	fld	QWORD PTR [edx][edi*8-24] ; Load second number
	fsub	QWORD PTR [ecx][edi*8-24] ; Subtract first number
	fld	QWORD PTR [edx][edi*8-16] ; Load second number
	fsub	QWORD PTR [ecx][edi*8-16] ; Subtract first number
	fxch	st(1)
	fld	QWORD PTR [edx][edi*8-8]  ; Load second number
	fsub	QWORD PTR [ecx][edi*8-8]  ; Subtract first number
	fxch	st(3)
	fstp	QWORD PTR [esi][edi*8-32] ; Save result
	fstp	QWORD PTR [esi][edi*8-24] ; Save result
	fstp	QWORD PTR [esi][edi*8-16] ; Save result
	fstp	QWORD PTR [esi][edi*8-8]  ; Save result
	sub	edi, 4			; U - Check loop counter
	jnz	short usublp		; V - Loop if necessary
	pop	edi			; U - Restore registers
	pop	esi			; V - Restore registers
	ret

; Do a normalized subtract
; ecx = src arg #1
; edx = src arg #2
; esi = dest arg
; edi = loop counter

nsub:	fninit
	push	ebp			; U - Save registers
	push	ebx			; V - Save registers
	fld	BIGVAL			; UV - Start process with no carry
	mov	ebx, norm_col_mults	; U - Address of the multipliers
	mov	ebp, scaled_numlit	; V - Computes big vs little word flag
	mov	eax, -1			; U - First word is a big word (-1)
nsublp:	norm_op_1d fsub			; Subtract and normalize 8 values
	sub	edi, 8			; Decrement loop counter
	JNZ_X	nsublp			; Loop til done
	fsub	BIGVAL
	mov	ecx, _DESTARG		; Address of squared number
	cmp	_PLUS1, 0
	jz	short nsub1
	fchs
nsub1:	fadd	QWORD PTR [ecx]		; Add in carry
	mov	DWORD PTR [ecx-4], 0	; Clear needs-normalize flag
	pop	ebx
	pop	ebp
	pop	edi
	pop	esi
	fstp	QWORD PTR [ecx]		; Save new value
	ret
gwsub1	ENDP


;;
;; Add and subtract two numbers without carry propogation.
;;

	PUBLIC	gwaddsubq1
gwaddsubq1 PROC NEAR
	mov	ecx, _SRCARG		; U - Address of first number
	mov	edx, _SRC2ARG		; V - Address of second number
	push	esi			; U - Save registers
	push	edi			; V - Save registers
	mov	eax, [ecx-4]		; U - Load needs-normalize counter
	mov	edi, [edx-4]		; V - Load needs-normalize counter
	push	ebp			; U - Save registers
	mov	esi, _DESTARG		; V - Address of destination #1
	mov	ebp, _DEST2ARG  	; U - Address of destination #2
	lea	eax, [eax+edi+1]	; V - Set new needs-normalize counter
	mov	edi, _FFTLEN		; U - Load loop counter
	mov	[esi-4], eax		; V - Store needs-normalize counter
	mov	[ebp-4], eax		; U - Store needs-normalize counter
	jmp	short uaddsublp		; V - Do the unnormalized add
gwaddsubq1 ENDP

;;
;; Add and subtract two numbers.  This will be done without carry propogation
;; (for extra speed) if it can be done safely.
;;

	PUBLIC	gwaddsub1
gwaddsub1 PROC NEAR

; See if add can be done without normalization

	mov	ecx, _SRCARG		; U - Address of first number
	mov	edx, _SRC2ARG		; V - Address of second number
	push	esi			; U - Save registers
	push	edi			; V - Save registers
	mov	eax, [ecx-4]		; U - Load needs-normalize counter
	mov	edi, [edx-4]		; V - Load needs-normalize counter
	push	ebp			; U - Save registers
	mov	esi, _DESTARG		; V - Address of destination
	lea	eax, [eax+edi+1]	; U - Set new needs-normalize counter
	mov	edi, extra_bits		; V - Load compare value
	mov	ebp, _DEST2ARG  	; U - Address of destination #2
	cmp	eax, edi		; V - Is normalization needed?
	mov	edi, _FFTLEN		; U - Load loop counter
	mov	[esi-4], eax		; V - Store needs-normalize counter
	mov	[ebp-4], eax		; U - Store needs-normalize counter
	jg	short naddsub		; V - Yes, do a normalized add

; Do an unnormalized add
; ecx = src arg #1
; edx = src arg #2
; esi = dest arg #1
; ebp = dest arg #2
; edi = loop counter

uaddsublp:
	fld	QWORD PTR [ecx][edi*8-32] ; Load first number
	fld	st(0)			  ; Dup first number
	fadd	QWORD PTR [edx][edi*8-32] ; Add in second number
	fxch	st(1)			  ; S0,A0
	fsub	QWORD PTR [edx][edi*8-32] ; Subtract out second number
	fld	QWORD PTR [ecx][edi*8-24] ; Load first number
	fld	st(0)			  ; Dup first number
	fadd	QWORD PTR [edx][edi*8-24] ; Add in second number
	fxch	st(1)			  ; S1,A1,S0,A0
	fsub	QWORD PTR [edx][edi*8-24] ; Subtract out second number
	fld	QWORD PTR [ecx][edi*8-16] ; Load first number
	fld	st(0)			  ; Dup first number
	fadd	QWORD PTR [edx][edi*8-16] ; Add in second number
	fxch	st(1)			  ; S2,A2,S1,A1,S0,A0
	fsub	QWORD PTR [edx][edi*8-16] ; Subtract out second number
	fld	QWORD PTR [ecx][edi*8-8]  ; Load first number
	fld	st(0)			  ; Dup first number
	fadd	QWORD PTR [edx][edi*8-8]  ; Add in second number
	fxch	st(7)			  ; A0,S3,S2,A2,S1,A1,S0,A3
	fstp	QWORD PTR [esi][edi*8-32] ; Save result
	fsub	QWORD PTR [edx][edi*8-8]  ; Subtract out second number
	fxch	st(5)			  ; S0,S2,A2,S1,A1,S3,A3
	fstp	QWORD PTR [ebp][edi*8-32] ; Save result
	fstp	QWORD PTR [ebp][edi*8-16] ; Save result
	fstp	QWORD PTR [esi][edi*8-16] ; Save result
	fstp	QWORD PTR [ebp][edi*8-24] ; Save result
	fstp	QWORD PTR [esi][edi*8-24] ; Save result
	fstp	QWORD PTR [ebp][edi*8-8]  ; Save result
	fstp	QWORD PTR [esi][edi*8-8]  ; Save result
	sub	edi, 4			; U - Check loop counter
	jnz	short uaddsublp		; V - Loop if necessary
	pop	ebp			; U - Restore registers
	pop	edi			; U - Restore registers
	pop	esi			; V - Restore registers
	ret

; Do a normalized add and subtract
; ecx = src arg #1
; edx = src arg #2
; esi = dest arg #1
; ebp = dest arg #2
; edi = loop counter

naddsub:fninit
	mov	loopcount1, edi		; U - Loop counter
	push	ebx			; V - Save registers
	fld	BIGVAL			; UV - Start process with no carry
	fld	BIGVAL			; UV - Start process with no carry
	mov	ebx, norm_col_mults	; U - Address of the multipliers
	mov	edi, scaled_numlit	; V - Computes big vs little word flag
	mov	eax, -1			; U - First word is a big word (-1)
naddsublp:
	norm_addsub_1d			; Add/sub and normalize 8 values
	sub	loopcount1, 8		; UU - Decrement loop counter
	JNZ_X	naddsublp		; V - Loop til done
	fsub	BIGVAL			; UV - Get carry without BIGVAL
	fxch	st(1)
	fsub	BIGVAL			; UV - Get carry without BIGVAL
	mov	ecx, _DESTARG		; U - Address of squared number
	mov	eax, _DEST2ARG		; V - Address of squared number
	cmp	_PLUS1, 0		; UV - Are we working mod 2^N+1
	jz	short naddsub1		; V - No, leave carry alone
	fchs				; UV - Yes, flip carry's sign bit
	fxch	st(1)
	fchs				; UV - Yes, flip carry's sign bit
	fxch	st(1)
naddsub1:
	fadd	QWORD PTR [ecx]		; UV - Add in carry
	fxch	st(1)
	fadd	QWORD PTR [eax]		; UV - Add in carry
	mov	DWORD PTR [ecx-4], 0	; UV - Clear needs-normalize flag
	mov	DWORD PTR [eax-4], 0	; UV - Clear needs-normalize flag
	pop	ebx			; U - Restore registers
	pop	ebp			; V - Restore registers
	pop	edi			; U - Restore registers
	pop	esi			; V - Restore registers
	fstp	QWORD PTR [eax]		; UV - Save new value
	fstp	QWORD PTR [ecx]		; UV - Save new value
	ret
gwaddsub1 ENDP

;;
;; Copy one number and zero some low order words.
;;

	PUBLIC	gwcopyzero1
gwcopyzero1 PROC NEAR
	push	esi			; U - Save registers
	push	edi			; V - Save registers
	mov	esi, _SRCARG		; U - Address of first number
	mov	edi, _DESTARG		; V - Address of destination
	mov	edx, _FFTLEN		; U - Total number of values
	mov	ecx, _SRC2ARG		; V - Values to zero
	sub	edx, ecx		; U - Number of values to copy
	mov	eax, [esi-4]		; V - Load needs-normalize counter
	mov	[edi-4], eax		; U - Store needs-normalize counter
	lea	esi, [esi+8*ecx]	; V - First value to copy
	sub	eax, eax		; U - Create the zero to store
	add	ecx, ecx		; V - Convert num values to num words
	rep	stosd			; UV - Zero words
	lea	ecx, [edx+edx]		; U - Convert num values to num words
	rep	movsd			; UV - Copy words
	pop	edi			; U - Restore registers
	pop	esi			; V - Restore registers
	ret
gwcopyzero1 ENDP

;;
;; Do a mod k*2^n+/-1
;;

	PUBLIC	gwprothmod1
gwprothmod1 PROC NEAR
	push	esi			; U - Save registers
	push	edi			; V - Save registers
	mov	esi, _SRCARG		; U - Address of first number
	prothmod_upper_prep_0d		; Prepare for prothmod
pmulp:	prothmod_upper_0d pmudn		; Divide upper values by k
	lea	esi, [esi-8*dist1]	; Next source word
	lea	edi, [edi+8*8]		; Next scratch words
	JMP_X	pmulp			; Loop til done
pmudn:	prothmod_lower_prep_0d		; Prepare for prothmod
pmllp:	prothmod_lower_0d		; Divide upper values by k
	lea	esi, [esi+8*dist1]	; Next source word
	lea	edi, [edi-8*8]		; Next scratch words
	dec	edx
	jz	short pmldn
	JMP_X	pmllp			; Loop til done
pmldn:	prothmod_final_0d		; Finish off the mod
	pop	edi			; U - Restore registers
	pop	esi			; V - Restore registers
	ret
gwprothmod1 ENDP
ENDIF

;;
;; Routines to do the normalization after a multiply
;;

IFNDEF PPRO
_simpnorm PROC NEAR

; Macro to loop through all the FFT values and apply the proper normalization
; routine.  Used whenever we are using an irrational-base FFT.

inorm	MACRO	lab, norm_macro
	LOCAL	ilp
	PUBLIC	lab
lab:	mov	esi, _DESTARG		;; Address of multiplied number
	mov	edi, _ADDIN_OFFSET	;; Get address to add value into
	fld	QWORD PTR [esi][edi]	;; Get the value
	fadd	_ADDIN_VALUE		;; Add in the requested value
	fstp	QWORD PTR [esi][edi]	;; Save the new value
	fldz				;; Init SUMOUT
	fsub	_ADDIN_VALUE		;; Do not include addin in sumout
	fld	BIGVAL			;; Start process with no carry
	mov	ebx, norm_col_mults	;; Address of the multipliers
	mov	edi, scaled_numlit	;; Used to compute big vs little words
	mov	ebp, scaled_numbig
	mov	eax, -1			;; First word is a big word (-1)
	mov	ecx, _FFTLEN		;; Load loop counter
	mov	edx, _FFTZERO		;; Count of words to NOT zero
ilp:	norm_macro			;; Normalize 8 values
	lea	esi, [esi+8*dist1]	;; Next source
	lea	ebx, [ebx+8*NMD]	;; Next set of 8 multipliers
	sub	ecx, 8			;; Test loop counter
	JNZ_X	ilp			;; Loop til done
	JMP_X	idn
	ENDM

; Macro to loop through all the FFT values and apply the proper normalization
; routine.  Used whenever we are using a rational-base FFT.

rnorm	MACRO	lab, norm_macro
	LOCAL	rlp
	PUBLIC	lab
lab:	mov	esi, _DESTARG		;; Address of multiplied number
	mov	edi, _ADDIN_OFFSET	;; Get address to add value into
	fld	QWORD PTR [esi][edi]	;; Get the value
	fadd	_ADDIN_VALUE		;; Add in the requested value
	fstp	QWORD PTR [esi][edi]	;; Save the new value
	fld	ttmp_ff_inv		;; Preload 2 / FFTLEN
	fldz				;; Init SUMOUT
	fsub	_ADDIN_VALUE		;; Do not include addin in sumout
	fld	BIGVAL			;; Start process with no carry
	mov	ecx, _FFTLEN		;; Load loop counter
	mov	edx, _FFTZERO		;; Count of words to NOT zero
rlp:	norm_macro			;; Normalize 8 values
	lea	esi, [esi+8*dist1]	;; Next source
	sub	ecx, 8			;; Test loop counter
	JNZ_X	rlp			;; Loop til done
	JMP_X	rdn
	ENDM

; The 16 different normalization routines

	rnorm	r1, norm_0d
	rnorm	r1e, norm_0d_e
	rnorm	r1c, norm_0d_c
	rnorm	r1ec, norm_0d_e_c
	rnorm	r1z, norm_0d_z
	rnorm	r1ze, norm_0d_z_e
	inorm	i1, norm_1d
	inorm	i1e, norm_1d_e
	inorm	i1c, norm_1d_c
	inorm	i1ec, norm_1d_e_c
	inorm	i1z, norm_1d_z
	inorm	i1ze, norm_1d_z_e

; Finish off the normalization process by add any carry to first values.
; Handle both the with and without two-to-phi array cases.

rdn:	mov	esi, _DESTARG		; Address of squared number
	norm012_0d
	jmp	short cmnend

idn:	mov	esi, _DESTARG		; Address of squared number
	mov	ebx, norm_col_mults	; Restart the column multipliers
	norm012_1d

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
_simpnorm ENDP
ENDIF

_TEXT32	ENDS
END
