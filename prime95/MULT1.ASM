; Copyright 1995-2000 Just For Fun Software, Inc., all rights reserved
; Author:  George Woltman
; Email: woltman@alum.mit.edu
;
; This routine implements a discrete weighted transform to quickly multiply
; two numbers.
;
; This code handles FFTs that use a simple linear memory model and
; the simplified normalization code.  FFT sizes from 32 doubles to
; 128 doubles are supported.
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
INCLUDE fft1a.mac
INCLUDE memory.mac
INCLUDE normal.mac

IFDEF PPRO
INCLUDE	lucas1p.mac
ENDIF

	dist1 =	8

;; All the FFT routines for each FFT length

_gw_ffts PROC NEAR
	EXPANDING = 1
	fft	32
	fft	32, p
	fft	40
	fft	48
;	fft	56
	fft	64
	fft	64, p
	fft	80
	fft	96
;	fft	112
	fft	128
	fft	128, p
_gw_ffts ENDP


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
naddlp:	simple_norm_op fadd		; Add and normalize 8 values
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
nsublp:	simple_norm_op fsub		; Subtract and normalize 8 values
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
	mov	addsubtemp, edi		; U - Loop counter
	push	ebx			; V - Save registers
	fld	BIGVAL			; UV - Start process with no carry
	fld	BIGVAL			; UV - Start process with no carry
	mov	ebx, norm_col_mults	; U - Address of the multipliers
	mov	edi, scaled_numlit	; V - Computes big vs little word flag
	mov	eax, -1			; U - First word is a big word (-1)
naddsublp:
	simple_norm_addsub		; Add/sub and normalize 8 values
	sub	addsubtemp, 8		; UU - Decrement loop counter
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
ENDIF

;; Routines to do the normalization after a multiply

_simpnorm PROC NEAR

; Normalize the results and multiply them by 3.  Used in P-1 factoring.
; As a safety measure, make sure the top bits of ERRCHK is AB93.  This
; makes it highly unlikely that ERRCHK will be corrupted and LL tests
; will all of a sudden start doing a multiply by 3.

nt3:	mov	eax, _ERRCHK		; Load variable
	and	eax, 0FFFF0000h		; Mask out bits
	cmp	eax, 0AB930000h		; See if the upper bits are correct
	JNE_X	noechk			; This should never happen
	fldz				; Init SUMOUT
	fld	BIGVAL			; Start process with no carry
	mov	ebx, norm_col_mults	; Address of the multipliers
	mov	esi, _DESTARG		; Address of multiplied number
	mov	edi, scaled_numlit	; Used to compute big vs little words
	mov	ebp, scaled_numbig
	mov	eax, -1			; First word is a big word (-1)
	mov	ecx, _FFTLEN		; Load loop counter
plp:	simple_normalize_times_3	; Normalize 8 values
	lea	esi, [esi+8*dist1]	; Next source
	lea	ebx, [ebx+8*NMD]	; Next set of 8 multipliers
	sub	ecx, 8			; Test loop counter
	JNZ_X	plp			; Loop til done
	JMP_X	ndn			; Go add in final carry

; Code to check for how close we are to incurring a round-off error

echk:	fld	_MAXERR			; Load old maximum error
	mov	ebx, norm_col_mults	; Address of the multipliers
	mov	esi, _DESTARG		; Address of multiplied number
	mov	ecx, _FFTLEN		; Loop counter
elp:	simple_error_check		; Error check eight values
	lea	esi, [esi+dist1]	; Next source
	lea	ebx, [ebx+NMD]		; Next multiplier
	dec	ecx
	jnz	short elp		; Loop til done
	fstp	_MAXERR			; Store maximum error

; Normalize the results

noechk:	fldz				; Init SUMOUT
	fld	BIGVAL			; Start process with no carry
	mov	ebx, norm_col_mults	; Address of the multipliers
	mov	esi, _DESTARG		; Address of multiplied number
	mov	edi, scaled_numlit	; Used to compute big vs little words
	mov	ebp, scaled_numbig
	mov	eax, -1			; First word is a big word (-1)
	mov	ecx, _FFTLEN		; Load loop counter
nlp:	simple_normalize		; Normalize 8 values
	lea	esi, [esi+8*dist1]	; Next source
	lea	ebx, [ebx+8*NMD]	; Next set of 8 multipliers
	sub	ecx, 8			; Test loop counter
	JNZ_X	nlp			; Loop til done
					; Add any carry to first values
ndn:	mov	esi, _DESTARG		; Address of multiplied number
	mov	ebx, norm_col_mults	; Restart the multipliers
	cmp	_PLUS1, 0
	jz	short nlp1
	fsub	BIGVAL
	fchs
	fadd	BIGVAL
nlp1:	simple_normalize_012
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
_simpnorm ENDP

_TEXT32	ENDS
END
