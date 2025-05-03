; Copyright 1995-2024 Mersenne Research, Inc.  All rights reserved
; Author:  George Woltman
; Email: woltman@alum.mit.edu
;

	TITLE   setup

	.686
	.MODEL	FLAT

INCLUDE	unravel.mac
INCLUDE extrn.mac
INCLUDE mult.mac
INCLUDE memory.mac
INCLUDE normal.mac

_TEXT SEGMENT

	flat_distances

;;
;; Add two numbers without carry propagation.  Caller can use this for
;; consecutive add or subtract operations.  However, the last operation
;; before a multiply must use the routine that will normalize data.
;;

PROCFL	gwaddq1
	ad_prolog 0,0,rsi,rdi
	mov	ecx, SRCARG		; Address of first number
	mov	edx, SRC2ARG		; Address of second number
	mov	esi, DESTARG		; Address of destination
	mov	edi, FFTLEN		; Load loop counter
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
	sub	edi, 4			; Check loop counter
	jnz	short uaddlp		; Loop if necessary
	ad_epilog 0,0,rsi,rdi
gwaddq1	ENDP

;;
;; Add two numbers with carry propagation
;;

PROCFL	gwadd1
	ad_prolog 0,0,rbx,rbp,rsi,rdi
	mov	ecx, SRCARG		; Address of first number
	mov	edx, SRC2ARG		; Address of second number
	mov	esi, DESTARG		; Address of destination
	mov	ebp, FFTLEN		; Load loop counter
	fld	BIGVAL			; Start process with no carry
	fld	BIGVAL
	mov	ebx, norm_col_mults	; Address of the multipliers
	mov	edi, norm_biglit_array	; Computes big vs little word flag
	sub	eax, eax		; Clear big/lit flags
naddlp:	norm_op_1d fadd			; Add and normalize 8 values
	sub	ebp, 8			; Decrement loop counter
	jnz	naddlp			; Loop til done
	mov	esi, DESTARG		; Address of squared number
	mov	ebx, norm_col_mults	; Address of the multipliers
	norm_op_1d_cleanup
	ad_epilog 0,0,rbx,rbp,rsi,rdi
gwadd1	ENDP


;;
;; Subtract two numbers without carry propagation.  Caller can use this for
;; consecutive add or subtract operations.  However, the last operation
;; before a multiply must use the routine that will normalize data.
;;

PROCFL	gwsubq1
	ad_prolog 0,0,rsi,rdi
	mov	ecx, SRCARG		; Address of first number
	mov	edx, SRC2ARG		; Address of second number
	mov	esi, DESTARG		; Address of destination
	mov	edi, FFTLEN		; Load loop counter
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
	sub	edi, 4			; Check loop counter
	jnz	short usublp		; Loop if necessary
	ad_epilog 0,0,rsi,rdi
gwsubq1	ENDP

;;
;; Subtract two numbers with carry propagation
;;

PROCFL	gwsub1
	ad_prolog 0,0,rbx,rbp,rsi,rdi
	mov	ecx, SRCARG		; Address of first number
	mov	edx, SRC2ARG		; Address of second number
	mov	esi, DESTARG		; Address of destination
	mov	ebp, FFTLEN		; Load loop counter
	fld	BIGVAL			; Start process with no carry
	fld	BIGVAL
	mov	ebx, norm_col_mults	; Address of the multipliers
	mov	edi, norm_biglit_array	; Computes big vs little word flag
	sub	eax, eax		; Clear big/lit flags
nsublp:	norm_op_1d fsub			; Subtract and normalize 8 values
	sub	ebp, 8			; Decrement loop counter
	jnz	nsublp			; Loop til done
	mov	esi, DESTARG		; Address of squared number
	mov	ebx, norm_col_mults	; Address of the multipliers
	norm_op_1d_cleanup
	ad_epilog 0,0,rbx,rbp,rsi,rdi
gwsub1	ENDP


;;
;; Add and subtract two numbers without carry propagation.
;;

PROCFL	gwaddsubq1
	ad_prolog 0,0,rbp,rsi,rdi
	mov	ecx, SRCARG		; Address of first number
	mov	edx, SRC2ARG		; Address of second number
	mov	esi, DESTARG		; Address of destination #1
	mov	ebp, DEST2ARG	  	; Address of destination #2
	mov	edi, FFTLEN		; Load loop counter
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
	sub	edi, 4			; Check loop counter
	jnz	short uaddsublp		; Loop if necessary
	ad_epilog 0,0,rbp,rsi,rdi
gwaddsubq1 ENDP

;;
;; Add and subtract two numbers with carry propagation
;;

loopcount1	EQU	DPTR [rsp+first_local]

PROCFL	gwaddsub1
	ad_prolog 4,0,rbx,rbp,rsi,rdi
	mov	ecx, SRCARG		; Address of first number
	mov	edx, SRC2ARG		; Address of second number
	mov	esi, DESTARG		; Address of destination
	mov	ebp, DEST2ARG	  	; Address of destination #2
	mov	eax, FFTLEN		; Load loop counter
	mov	loopcount1, eax		; Save loop counter
	fld	BIGVAL			; Start process with no carry
	fld	BIGVAL
	fld	BIGVAL
	fld	BIGVAL
	mov	ebx, norm_col_mults	; Address of the multipliers
	mov	edi, norm_biglit_array	; Computes big vs little word flag
	sub	eax, eax		; Clear big/lit flag
naddsublp:
	norm_addsub_1d			; Add/sub and normalize 4 values
	sub	loopcount1, 4		; Decrement loop counter
	jnz	naddsublp		; Loop til done
	mov	esi, DESTARG		; Address of squared number
	mov	ebp, DEST2ARG		; Address of squared number
	mov	ebx, norm_col_mults	; Address of the multipliers
	norm_addsub_1d_cleanup
	ad_epilog 4,0,rbx,rbp,rsi,rdi
gwaddsub1 ENDP

;;
;; Mul by a small value with carry propagation
;;

PROCFL	gwmuls1
	ad_prolog 0,0,rbx,rbp,rsi,rdi
	mov	esi, DESTARG		; Address of destination
	fld	DBLARG			; Load small value
	fmul	NORM012_FF		; Mul by two-to-minus-phi fudge
	fstp	TMP5			; Save multiplier
	mov	ebp, FFTLEN		; Load loop counter
	fld	BIGVAL			; Start process with no carry
	fld	BIGVAL
	mov	ebx, norm_col_mults	; Address of the multipliers
	mov	edi, norm_biglit_array	; Computes big vs little word flag
	sub	eax, eax		; Clear big/lit flags
nmullp:	norm_smallmul_1d		; Mul and normalize 8 values
	sub	ebp, 8			; Decrement loop counter
	jnz	nmullp			; Loop til done
	mov	esi, DESTARG		; Address of squared number
	mov	ebx, norm_col_mults	; Address of the multipliers
	norm_op_1d_cleanup
	ad_epilog 0,0,rbx,rbp,rsi,rdi
gwmuls1	ENDP


;;
;; Routines to do the normalization after a multiply
;;

; Macro to loop through all the FFT values and apply the proper normalization
; routine.  Used whenever we are using an irrational-base FFT.

inorm	MACRO	lab, ttp, echk, const
	LOCAL	ilp
	PROCFL	lab
	int_prolog 0,0,4
	mov	esi, DESTARG		;; Address of multiplied number
	fld	MAXERR			;; Load MAXERR
	fldz				;; Init SUMOUT
	mov	edi, ADDIN_OFFSET	;; Get address to add value into
	fld	QWORD PTR [esi][edi]	;; Get the value
	fadd	ADDIN_VALUE		;; Add in the requested value
	fstp	QWORD PTR [esi][edi]	;; Save the new value
	fsub	ADDIN_VALUE		;; Do not include addin in sumout
	fld	BIGVAL			;; Init carry #1
	fld	BIGVAL			;; Init carry #2
	mov	ebx, norm_col_mults	;; Address of the multipliers
	mov	edi, norm_biglit_array	;; Big/lit array pointer
	sub	eax, eax		;; Clear big/lit flag
	mov	ecx, FFTLEN		;; Load loop counter
ilp:	norm_1d ttp, echk, const	;; Normalize 8 values
	lea	esi, [esi+64]		;; Next source
ttp	lea	ebx, [ebx+128]		;; Next set of 8 multipliers
ttp	lea	edi, [edi+4]		;; Next big/lit array ptr
	sub	ecx, 8			;; Test loop counter
	jnz	ilp			;; Loop til done
	mov	esi, DESTARG		;; Address of multiplied number
	mov	edi, ADDIN_OFFSET	;; Get address to add value into
	fld	QWORD PTR [esi][edi]	;; Get the value
	fadd	POSTADDIN_VALUE		;; Add in the requested value
	fstp	QWORD PTR [esi][edi]	;; Save the new value
	jmp	idn			;; Go to normal end code
&lab	ENDP
	ENDM

; Do a zero-padded normalization.

zpnorm	MACRO	lab, ttp, echk, const
	LOCAL	ilp
	PROCFL	lab
	int_prolog 0,0,4

	c_call	ZPAD_SUB7		;; Subtract 7 ZPAD words from lowest FFT words

	mov	esi, DESTARG		;; Address of multiplied number
	fld	MAXERR			;; Load MAXERR
	fldz				;; Init SUMOUT
	fld	BIGVAL			;; Init traditional carry
	fldz				;; Init high FFT data carry
	mov	ebx, norm_col_mults	;; Address of the multipliers
	mov	edi, norm_biglit_array	;; Big/lit array pointer
	sub	eax, eax		;; Clear big/lit flag
	mov	ecx, FFTLEN		;; Load loop counter
ilp:	norm_1d_zpad ttp, echk, const	;; Normalize 4 values
	lea	esi, [esi+4*8]		;; Next source
ttp	lea	ebx, [ebx+4*16]		;; Next set of 4 multipliers
ttp	lea	edi, [edi+2]		;; Next big/lit array ptr
	sub	ecx, 4			;; Test loop counter
	jnz	ilp			;; Loop til done
const	jmp	zpcdn			;; Go to zero padded FFT end code
no const jmp	zpdn			;; Go to zero padded FFT end code
&lab	ENDP
	ENDM

; The 16 different normalization routines

	inorm	r1, noexec, noexec, noexec
	inorm	r1e, noexec, exec, noexec
	inorm	r1c, noexec, noexec, exec
	inorm	r1ec, noexec, exec, exec
	inorm	i1, exec, noexec, noexec
	inorm	i1e, exec, exec, noexec
	inorm	i1c, exec, noexec, exec
	inorm	i1ec, exec, exec, exec

	zpnorm	r1zp, noexec, noexec, noexec
	zpnorm	r1zpe, noexec, exec, noexec
	zpnorm	r1zpc, noexec, noexec, exec
	zpnorm	r1zpec, noexec, exec, exec
	zpnorm	i1zp, exec, noexec, noexec
	zpnorm	i1zpe, exec, exec, noexec
	zpnorm	i1zpc, exec, noexec, exec
	zpnorm	i1zpec, exec, exec, exec

; Common code to finish off the one-pass FFTs normalization.  The
; Windows 64-bit ABI frowns on us jumping from one procedure into another.
; However, my reading of the spec is that as long as the two procedures have
; identical prologs then stack unwinding for exception handling will work OK.
; Of course, this x87 is not linked in to a Windows 64-bit executable, but
; we make the code identical to xmult1ax.asm for consistency sake.

PROCF	__common_norm1_end_code

	;; Dummy prolog to match normalization code
	int_prolog 0,0,4

; Finish off the normalization process by adding any carry to first values.
; FPU stack is now: carry #2, carry #1, sumout, maxerr

zpcdn:	mov	esi, DESTARG		; Address of squared number
	mov	edi, norm_biglit_array	; Address of the big/little flags array
	mov	ebx, norm_col_mults	; Restart the column multipliers
	sub	eax, eax
	norm012_1d_zpad exec		; Add in carries
	jmp	cmnend			; All done, go cleanup

zpdn:	mov	esi, DESTARG		; Address of squared number
	mov	edi, norm_biglit_array	; Address of the big/little flags array
	mov	ebx, norm_col_mults	; Restart the column multipliers
	sub	eax, eax
	norm012_1d_zpad noexec		; Add in carries
	jmp	cmnend			; All done, go cleanup

idn:	mov	esi, DESTARG		; Address of squared number
	norm_top_carry_1d		; Adjust top carry when k > 1
	sub	eax, eax		; Clear big/lit flag
	mov	edi, norm_biglit_array	; Big/lit array pointer
	mov	ebx, norm_col_mults	; Restart the column multipliers
	norm012_1d noexec

; Return unnormalized SUMOUT -- no one uses it
; FPU stack is now: sumout, maxerr

cmnend:	fstp	QWORD PTR [esi-24]	; Save sum of FFT outputs
	fstp	MAXERR			; Save new maximum roundoff error

; Return

	int_epilog 0,0,4
__common_norm1_end_code ENDP

_TEXT	ENDS
END
