; Copyright 2001-2003 Just For Fun Software, Inc., all rights reserved
; Author:  George Woltman
; Email: woltman@alum.mit.edu
;
; This routine implements a discrete weighted transform to quickly multiply
; two numbers.
;
; This code uses Pentium 4's SSE2 instructions for very fast FFTs.
; FFT sizes of 40K and above are supported.
; This code does two passes, 11 levels on the second pass.
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
INCLUDE xfft3.mac
INCLUDE	xlucas.mac
INCLUDE xmult.mac
INCLUDE xnormal.mac

EXTRN	xpass2_11_levels:PROC

;; All the FFT routines for each FFT length

_xmm_gw_ffts3 PROC NEAR
	EXPANDING = 2
	xfft	40K
	xfft	48K
	xfft	56K
	xfft	64K
	xfft	80K
	xfft	96K
	xfft	112K
	xfft	128K
	xfft	160K
	xfft	192K
	xfft	224K
	xfft	256K
	xfft	320K
	xfft	384K
	xfft	448K
	xfft	512K
	EXPANDING = 3
	xfftclm	640K, 4, 1
	xfftclm	640K, 2, 0
	xfftclm	768K, 4, 1
	xfftclm	768K, 2, 0
	xfftclm	896K, 4, 1
	xfftclm	896K, 2, 1
	xfftclm	896K, 1, 0
	xfftclm	1024K, 4, 1
	xfftclm	1024K, 2, 1
	xfftclm	1024K, 1, 0
	xfftclm	1280K, 4, 1
	xfftclm	1280K, 2, 1
	xfftclm	1280K, 1, 0
	xfftclm	1536K, 4, 1
	xfftclm	1536K, 2, 1
	xfftclm	1536K, 1, 0
	xfftclm	1792K, 4, 1
	xfftclm	1792K, 2, 1
	xfftclm	1792K, 1, 1
	xfftclm	1792K, 0, 0
	xfftclm	2048K, 4, 1
	xfftclm	2048K, 2, 1
	xfftclm	2048K, 1, 1
	xfftclm	2048K, 0, 0
	xfftclm	2560K, 2, 1
	xfftclm	2560K, 1, 1
	xfftclm	2560K, 0, 0
	xfftclm	3072K, 2, 1
	xfftclm	3072K, 1, 1
	xfftclm	3072K, 0, 0
	xfftclm	3584K, 2, 1
	xfftclm	3584K, 1, 1
	xfftclm	3584K, 0, 0
	xfftclm	4096K, 2, 1
	xfftclm	4096K, 1, 1
	xfftclm	4096K, 0, 0

; Common code to finish the FFT by restoring state and returning.

gw_finish_fft_3:
	mov	DWORD PTR [esi-28], 3	; Set has-been-FFTed flags
	xfft_1_ret

; Common code to finish up multiplies

; Split the accumulated carries into two carries - a high carry and a
; low carry.  Handle both the with and without two-to-phi array cases.

gw_split_carries_3:
	sub	ecx, ecx		; Clear big/little flag
	mov	edx, count3		; Load 3 section counts
	mov	loopcount1, edx		; Save for later
	mov	edx, addcount1		; Load count of carry rows
	mov	esi, carries		; Addr of the carries
	mov	ebp, edx		; Compute addr of the high carries
	shl	ebp, 6
	add	ebp, esi
	mov	edi, norm_biglit_array	; Addr of the big/little flags array
	mov	edx, norm_grp_mults	; Addr of the group multipliers
	xnorm012_2d_part1
ilp0:	mov	eax, loopcount1		; Get list of counts
	mov	ebx, eax		; Form count for this section
	and	ebx, 03FFh
	JZ_X	spldn			; No rows to do.  We're all done!
	mov	loopcount2, ebx		; Save count of carry rows this section
	shr	eax, 10			; Move counts list along
	mov	loopcount1, eax
	shl	ebx, 6			; Compute addr of the last carries row
	add	ebx, esi
	sub	eax, eax		; Clear big/little flag
	xnorm012_2d_part2
ilp1:	mov	ebx, norm_col_mults	; Addr of the column multipliers
	xnorm012_2d			; Split carries for one cache line
	mov	ebx, count2		; Cache lines in each pass1 loop
	lea	esi, [esi+64]		; Next carries pointer
	lea	ebp, [ebp+64]		; Next high carries pointer
	cmp	_NUMLIT, 0		; Don't bump these two pointers
	je	short iskip		; for rational FFTs
	lea	edx, [edx+128]		; Next group multiplier
	lea	edi, [edi+ebx*4]	; Next big/little flags pointer
iskip:	sub	loopcount2, 1		; Test loop counter
	JNZ_X	ilp1			; Next carry row in section
	JMP_X	ilp0			; Next section
spldn:	ret

; Finish the multiply by adding the high and low carries (possibly FFTed)
; into the FFT data.

gw_finish_mult_3:
	mov	ebp, carries		; Addr of the low carries
	mov	esi, _DESTARG		; Addr of FFT data
	mov	edx, addcount1		; Load loop counter
	mov	ebx, edx		; Compute addr of the high carries
	shl	ebx, 6
	add	ebx, ebp
	movapd	xmm4, XMM_BIGVAL
	cmp	zero_fft, 1		; Is this a zero high words fft?
	JE_X	idz			; Yes, do special add in of carries
ilp2:	xnorm012_2d_addin noexec	; Add carries to FFT data
	add	esi, pass1blkdst	; Next carries pointer
	lea	ebp, [ebp+64]		; Next carries pointer
	lea	ebx, [ebx+64]		; Next high carries pointer
	sub	edx, 1			; Test loop counter
	JNZ_X	ilp2
	JMP_X	cmnend
idz:	xnorm012_2d_addin exec		; Add carries to FFT data
	add	esi, pass1blkdst	; Next carries pointer
	lea	ebp, [ebp+64]		; Next carries pointer
	lea	ebx, [ebx+64]		; Next high carries pointer
	sub	edx, 1			; Test loop counter
	JNZ_X	idz
	mov	zero_fft, 0		; Clear zero-high-words-fft flag
	mov	esi, _DESTARG		; Addr of FFT data
	sub	_FFTZERO, esi		; Restore FFTZERO addresses
	sub	_FFTZERO+4, esi
	sub	_FFTZERO+8, esi
	sub	_FFTZERO+12, esi
	sub	_FFTZERO+16, esi
	sub	_FFTZERO+20, esi
	sub	_FFTZERO+24, esi
	sub	_FFTZERO+28, esi

; Clear needs-normalize counter, set FFT-started flag

cmnend:	mov	esi, _DESTARG		; Addr of FFT data
	mov	DWORD PTR [esi-4], 0
	mov	eax, _POSTFFT		; Set FFT started flag
	mov	DWORD PTR [esi-28], eax

; Normalize SUMOUT value by multiplying by 1 / (fftlen/2).

	_movsd	xmm7, XMM_SUMOUT	; Add together the two partial sumouts
	addsd	xmm7, XMM_SUMOUT+8
	mulsd	xmm7, ttmp_ff_inv
	_movsd	[esi-24], xmm7		; Save sum of FFT outputs
	_movsd	xmm6, XMM_MAXERR	; Compute new maximum error
	maxsd	xmm6, XMM_MAXERR+8
	_movsd	_MAXERR, xmm6

; Test if the sum of the output values is an error (such as infinity or NaN)

	fld	QWORD PTR [esi-24]	; Load sumout
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

_xmm_gw_ffts3 ENDP

_TEXT32	ENDS
END
