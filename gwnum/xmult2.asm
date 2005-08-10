; Copyright 2001-2005 Just For Fun Software, Inc., all rights reserved
; Author:  George Woltman
; Email: woltman@alum.mit.edu
;
; This routine implements a discrete weighted transform to quickly multiply
; two numbers.
;
; This code uses Pentium 4's SSE2 instructions for very fast FFTs.
; FFT sizes between than 5K and 128K doubles are supported.
; This code does two passes, 8 levels on the second pass.
;
; You will not stand a chance of understanding any of this code without
; thoroughly familiarizing yourself with fast fourier transforms.  This
; code was adapted from an algorithm described in Richard Crandall's article
; on Discrete Weighted Transforms and Large-Integer Arithmetic.
;

	TITLE   setup

IFNDEF X86_64
	.686
	.XMM
	.MODEL	FLAT
ENDIF

INCLUDE	unravel.mac
INCLUDE extrn.mac
INCLUDE xfft2.mac
INCLUDE	xlucas.mac
INCLUDE xmult.mac
INCLUDE xnormal.mac
INCLUDE xpass1.mac
INCLUDE xpass1sc.mac

IFNDEF AMD
PUBLIC	xgw_finish_fft
PUBLIC	xgw_carries
PUBLIC	xgw_finish_mult
ELSE
EXTRN	xgw_finish_fft:PROC
EXTRN	xgw_carries:PROC
EXTRN	xgw_finish_mult:PROC
ENDIF

EXTRNP	xpass2_8_levels
EXTRNP	xpass2_8_levels_p

_TEXT SEGMENT

;; Distance between two pass 2 data blocks.  Pass 2 does 8 FFT levels
;; 2 sets of data (2 * 2^8 complex values = 2^10 doubles = 8KB).

PREFETCHING = 1
blkdst = (8192+128+GAP2_8_4)
xpass2_levels = 8

;; All the FFT routines for each FFT length

PROCP	_xmm_gw_ffts2
	EXPANDING = 2
;	xfft	5K
;	xfft	6K
;	xfft	7K
;	xfft	8K
	xfft	10K
	xfft	12K
	xfft	12Kp
	xfft	14K
	xfft	16K
	xfft	16Kp
	xfft	20K
	xfft	24K
	xfft	24Kp
	xfft	28K
	xfft	32K
	xfft	32Kp
;	xfft	40K
;	xfft	48K
;	xfft	48Kp
;	xfft	56K
;	xfft	64K
;	xfft	64Kp

;	xfftclm	80K, 4
;	xfftclm	96K, 4
;	xfftclm	96Kp, 4
;	xfftclm	112K, 4
;	xfftclm	128K, 4
;	xfftclm	128Kp, 4
;	xfftclm	160K, 4
;	xfftclm	192K, 4
;	xfftclm	192Kp, 4
;	xfftclm	224K, 4
;	xfftclm	256K, 4
;	xfftclm	256Kp, 4
;	xfftclm	320K, 4
;	xfftclm	384K, 4
;	xfftclm	384Kp, 4
;	xfftclm	448K, 4
;	xfftclm	512K, 4
;	xfftclm	512Kp, 4
allfft	xfftclm	640K, 1
allfft	xfftclm	768K, 1
allfft	xfftclm	768Kp, 1
allfft	xfftclm	896K, 1
allfft	xfftclm	1024K, 1
allfft	xfftclm	1024Kp, 1

IFNDEF AMD

; Common code to finish the FFT by restoring state and returning.

xgw_finish_fft:
	mov	DWORD PTR [rsi-28], 3	; Set has-been-FFTed flags
	xfft_1_ret

; Common code to finish up multiplies

; Split the accumulated carries into two carries - a high carry and a
; low carry.  Handle both the with and without two-to-phi array cases.
; Add these carries back into the FFT data.

xgw_carries:
	cmp	_ZERO_PADDED_FFT, 0	; Special case the zero padded FFT case
	jne	xgw_carries_zpad
	mov	rsi, carries		; Addr of the carries
	mov	ebx, addcount1		; Compute addr after carries
	shl	rbx, 6
	add	rbx, rsi
	xnorm012_2d_part1
	mov	rbp, _DESTARG		; Addr of the FFT data
	mov	rdi, norm_biglit_array	; Addr of the big/little flags array
	mov	rdx, norm_grp_mults	; Addr of the group multipliers
	sub	rcx, rcx		; Clear big/little flag
	mov	eax, count3		; Load 3 section counts
	mov	loopcount1, eax		; Save for later
ilp0:	mov	eax, loopcount1		; Get list of counts
	mov	ebx, eax		; Form count for this section
	and	rbx, 07FFh
	jz	spldn			; No rows to do.  We're all done!
	mov	loopcount2, ebx		; Save count of carry rows this section
	shr	eax, 11			; Move counts list along
	mov	loopcount1, eax
	shl	rbx, 6			; Compute addr of the last carries row
	add	rbx, rsi
	xnorm012_2d_part2
	sub	rax, rax		; Clear big/little flag
ilp1:	mov	rbx, norm_col_mults	; Addr of the column multipliers
	xnorm012_2d			; Split carries for one cache line
	mov	ebx, cache_line_multiplier; Cache lines in each pass1 loop
	lea	rsi, [rsi+64]		; Next carries pointer
	add	rbp, pass1blkdst	; Next FFT data pointer
	cmp	_RATIONAL_FFT, 0	; Don't bump these two pointers
	jne	short iskip		; for rational FFTs
	lea	rdx, [rdx+128]		; Next group multiplier
	lea	rdi, [rdi+rbx*4]	; Next big/little flags pointer
iskip:	sub	loopcount2, 1		; Test loop counter
	jnz	ilp1			; Next carry row in section
	jmp	ilp0			; Next section
spldn:	mov	zero_fft, 0		; Clear zero-high-words-fft flag
	ret

xgw_carries_zpad:
	mov	rsi, carries		; Addr of the carries
	mov	ebx, addcount1		; Compute addr after carries
	shl	rbx, 6
	add	rbx, rsi
	mov	rsi, _DESTARG		; Addr of the FFT data
	mov	rdi, norm_biglit_array	; Addr of the big/little flags array
	mov	rbp, norm_col_mults	; Addr of the group multipliers
	cmp	const_fft, 0		; Call correct part1 macro
	je	c2a			; Jump if not const
	xnorm012_2d_zpad_part1 exec
	jmp	c2b
c2a:	xnorm012_2d_zpad_part1 noexec
c2b:	mov	rsi, carries		; Addr of the carries
	mov	rbp, _DESTARG		; Addr of the FFT data
	mov	rdi, norm_biglit_array	; Addr of the big/little flags array
	mov	rdx, norm_grp_mults	; Addr of the group multipliers
	mov	eax, count3		; Load 3 section counts
	mov	loopcount1, eax		; Save for later
zlp0:	mov	eax, loopcount1		; Get list of counts
	mov	ebx, eax		; Form count for this section
	and	rbx, 07FFh
	jz	zpldn			; No rows to do.  We're all done!
	mov	loopcount2, ebx		; Save count of carry rows this section
	shr	eax, 11			; Move counts list along
	mov	loopcount1, eax
	shl	rbx, 6			; Compute addr of the last carries row
	add	rbx, rsi
	xnorm012_2d_zpad_part2
	sub	rax, rax		; Clear big/little flag
zlp1:	mov	rbx, norm_col_mults	; Addr of the column multipliers
	cmp	const_fft, 0		; Call correct zpad macro
	je	c2c			; Jump if not const
	xnorm012_2d_zpad exec		; Split carries for one cache line
	jmp	c2d
c2c:	xnorm012_2d_zpad noexec		; Split carries for one cache line
c2d:	mov	ebx, cache_line_multiplier; Cache lines in each pass1 loop
	lea	rsi, [rsi+64]		; Next carries pointer
	add	rbp, pass1blkdst	; Next FFT data pointer
	cmp	_RATIONAL_FFT, 0	; Don't bump these two pointers
	jne	short zskip		; for rational FFTs
	lea	rdx, [rdx+128]		; Next group multiplier
	lea	rdi, [rdi+rbx*4]	; Next big/little flags pointer
zskip:	sub	loopcount2, 1		; Test loop counter
	jnz	zlp1			; Next carry row in section
	jmp	zlp0			; Next section
zpldn:	mov	const_fft, 0		; Clear mul-by-const-fft flag
	ret

; Finish the multiply

xgw_finish_mult:

; Set FFT-started flag

	mov	rsi, _DESTARG		; Addr of FFT data
	mov	eax, _POSTFFT		; Set FFT started flag
	mov	DWORD PTR [rsi-28], eax

; Normalize SUMOUT value by multiplying by 1 / (fftlen/2).

	_movsd	xmm7, XMM_SUMOUT	; Add together the two partial sumouts
	addsd	xmm7, XMM_SUMOUT+8
	mulsd	xmm7, ttmp_ff_inv
	_movsd	[rsi-24], xmm7		; Save sum of FFT outputs
	_movsd	xmm6, XMM_MAXERR	; Compute new maximum error
	maxsd	xmm6, XMM_MAXERR+8
	_movsd	_MAXERR, xmm6

; Return

exit:	xrestore_registers
	ret
ENDIF

ENDPP	_xmm_gw_ffts2

_TEXT	ENDS
END
