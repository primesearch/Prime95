; Copyright 1995-2000 Just For Fun Software, Inc., all rights reserved
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
INCLUDE	lucas1.mac
INCLUDE pfa.mac
INCLUDE mult.mac
INCLUDE fft3.mac
INCLUDE memory.mac

IFDEF PPRO
INCLUDE	lucas1p.mac
ENDIF

EXTRNP	pass1_6_levels_fft
EXTRNP	pass1_6_levels_unfft
EXTRNP	pass1_6a_levels_fft
EXTRNP	pass1_6a_levels_unfft
EXTRNP	pass1_6_levels_fftp
EXTRNP	pass1_6_levels_unfftp
EXTRNP	pass1_7_levels_fft
EXTRNP	pass1_7_levels_unfft
EXTRNP	pass1_7a_levels_fft
EXTRNP	pass1_7a_levels_unfft
EXTRNP	pass1_7_levels_fftp
EXTRNP	pass1_7_levels_unfftp
EXTRNP	pass2_8_levels_type_123
EXTRNP	pass2_8_levels_type_123p
EXTRNP	pass2_8_levels_type_4
EXTRNP	pass2_8_levels_type_4p

	very_convoluted_distances

;; All the FFT routines for each FFT length

_gw_ffts PROC NEAR
	EXPANDING = 2
	fft	80K
	fft	96K
	fft	112K
	fft	128K
	fft	128K p
	fft	160K
	fft	192K
	fft	224K
	fft	256K
	fft	256K p
	fft	320K
	fft	384K
	fft	448K
	fft	512K
	fft	512K p
	fft	640K
	fft	768K
	fft	896K
	fft	1024K
	fft	1024K p
	fft	1280K
	fft	1536K
	fft	1792K
	fft	2048K
	fft	2048K p
	fft	2560K
	fft	3072K
	fft	3584K
	fft	4096K
	fft	4096K p
_gw_ffts ENDP

INCLUDE normal.mac
PURGE simple_normalize
PURGE simple_normalize_012
PURGE simple_normalize_times_3
PURGE simple_norm_op
PURGE simple_norm_addsub
PURGE simple_error_check

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
	lea	esi, [esi+dist8]	; U - Next source
	lea	ebx, [ebx+8*NMD]	; V - Next set of 8 column multipliers
	add	ecx, 65536/16*65536	; U - Decrement loop counter
	JNC_X	plp2			; V - Loop til done
	lea	esi, [esi-16*dist8+dist128]; U - Next source
	lea	edx, [edx+NMD]		; V - Next group multiplier
	fcompp				; UV - Pop group multipliers
	fcompp				; UV - Pop group multipliers
	test	esi, 31			; U - Test loop counter
	JNZ_X	plp3			; V - Loop til done
	sub	ecx, 32768-1		; U - Decrement loop counter
	JS_X	ndn			; V - Loop til done
	sub	edi, scaling_ff		;UU - Make edi accurate again
	test	ecx, 2-1		;UV - Test next loop counter
	JNZ_X	plp3			;*V - Loop til done
	test	ecx, 16*2-1		;UV - Test next loop counter
	lea	esi, [esi-8*dist128+dist1K]; U - Next source
	JNZ_X	plp3			; V - Loop til done
	test	ecx, 8*16*2-1		;UV - Test next loop counter
	lea	esi, [esi-16*dist1K+dist16K]; U - Next source
	JNZ_X	plp3			; V - Loop til done
	sub	edi, scaling_ff2	;UU - Make edi accurate again
	test	ecx, 8*8*16*2-1		;UV - Test next loop counter
	lea	esi, [esi-8*dist16K+dist128K]; U - Next source
	JNZ_X	plp3			; V - Loop til done
	lea	esi, [esi-8*dist128K+dist1M]; U - Next source
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
	add	ecx, 65536/16*65536	; U - Decrement loop counter
	JNC_X	elp2			; V - Loop til done
	lea	esi, [esi-16*dist8+dist128]; U - Next source
	lea	edx, [edx+NMD]		; V - Next group multiplier
	fcompp				; UV - Pop group multipliers
	test	esi, 31			; U - Test loop counter
	JNZ_X	elp3			; V - Loop til done
	sub	ecx, 32768-1		; U - Decrement loop counter
	js	short edn		; V - Loop til done
	sub	edi, scaling_ff		;UU - Make edi accurate again
	test	ecx, 2-1		;UV - Test next loop counter
	JNZ_X	elp3			;*V - Loop til done
	test	ecx, 16*2-1		;UV - Test next loop counter
	lea	esi, [esi-8*dist128+dist1K]; U - Next source
	JNZ_X	elp3			; V - Loop til done
	test	ecx, 8*16*2-1		;UV - Test next loop counter
	lea	esi, [esi-16*dist1K+dist16K]; U - Next source
	JNZ_X	elp3			; V - Loop til done
	sub	edi, scaling_ff2	;UU - Make edi accurate again
	test	ecx, 8*8*16*2-1		;UV - Test next loop counter
	lea	esi, [esi-8*dist16K+dist128K]; U - Next source
	JNZ_X	elp3			; V - Loop til done
	lea	esi, [esi-8*dist128K+dist1M]; U - Next source
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
	lea	esi, [esi+dist8]	; U - Next source
	lea	ebx, [ebx+8*NMD]	; V - Next set of 8 column multipliers
	add	ecx, 65536/16*65536	; U - Decrement loop counter
	JNC_X	nlp2			; V - Loop til done
	lea	esi, [esi-16*dist8+dist128]; U - Next source
	lea	edx, [edx+NMD]		; V - Next group multiplier
	fcompp				; UV - Pop group multipliers
	fcompp				; UV - Pop group multipliers
	test	esi, 31			; U - Test loop counter
	JNZ_X	nlp3			; V - Loop til done
	sub	ecx, 32768-1		; U - Decrement loop counter
	js	short ndn		; V - Loop til done
	sub	edi, scaling_ff		;UU - Make edi accurate again
	test	ecx, 2-1		;UV - Test next loop counter
	JNZ_X	nlp3			;*V - Loop til done
	test	ecx, 16*2-1		;UV - Test next loop counter
	lea	esi, [esi-8*dist128+dist1K]; U - Next source
	JNZ_X	nlp3			; V - Loop til done
	test	ecx, 8*16*2-1		;UV - Test next loop counter
	lea	esi, [esi-16*dist1K+dist16K]; U - Next source
	JNZ_X	nlp3			; V - Loop til done
	sub	edi, scaling_ff2	;UU - Make edi accurate again
	test	ecx, 8*8*16*2-1		;UV - Test next loop counter
	lea	esi, [esi-8*dist16K+dist128K]; U - Next source
	JNZ_X	nlp3			; V - Loop til done
	lea	esi, [esi-8*dist128K+dist1M]; U - Next source
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
