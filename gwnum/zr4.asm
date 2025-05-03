; Copyright 2011-2024 Mersenne Research, Inc.  All rights reserved
; Author:  George Woltman
; Email: woltman@alum.mit.edu
;
; Assemble the small one pass AVX-512 FFTs.
;

	TITLE   setup

zfft_type TEXTEQU <r4>

INCLUDE	unravel.mac
INCLUDE extrn.mac
INCLUDE zarch.mac
INCLUDE zbasics.mac
INCLUDE zmult.mac
INCLUDE zonepass.mac
INCLUDE zr4.mac
INCLUDE zr4dwpnpass1sc.mac

_TEXT SEGMENT

EXTRN zcomplex_square_opcode:PROC
EXTRN zcomplex_mult_opcode:PROC
EXTRN zcomplex_mulf_opcode:PROC

;; To reduce the executable size, a common routine handles the middle section of one pass FFTs
;; On input and output the 128-byte cache line pairs hold these data values:
;;	0	+incr	...	8	+incr	...
;;	...
;;	7
;;	16	...
;;	32	...
;;	...
;; rcx = number of FFT sections
;; rsi = FFT data pointer (preserved or preserved with DEST2ARG added in)

zreal_onepass_middle PROC
	mov	rdi, rsi		;; Save source pointer
	cmp	ffttype, 2		;; Case off of FFT type
	je	type2r
	ja	type34r

	;; Type 1 FFT (forward FFT only)
	zr8_sixteen_reals_eight_complex_fft_final rsi, 8*128, 128, 2*128, 4*128
	zr8_eight_complex_fft_final_preload
	jmp	loop1r			;; Jump to complex middle code

	;; Type 2 FFT (squaring)
type2r:	mov	r8, DEST2ARG
	zr8_sixteen_reals_eight_complex_with_square rsi, 8*128, 128, 2*128, 4*128
	zr8_eight_complex_with_square_preload
	jmp	loop2r			;; Jump to complex middle code

	;; FFT types 3 and 4
type34r:mov	rbp, DIST_TO_MULSRCARG
	cmp	ffttype, 4
	je	type4r

	;; Type 3 FFT (forward FFT and multiply)
	mov	r8, DEST2ARG
	zr8_sixteen_reals_eight_complex_with_mult rsi, 8*128, 128, 2*128, 4*128
	zr8_eight_complex_with_mult_preload
	jmp	loop3r			;; Jump to complex middle code

	;; Type 4 FFT (forward FFT and multiply)
type4r:	zr8_sixteen_reals_eight_complex_with_mulf rsi, 8*128, 128, 2*128, 4*128
	zr8_eight_complex_with_mulf_preload
	jmp	loop4r			;; Jump to complex middle code
zreal_onepass_middle ENDP

zcomplex_onepass_middle PROC
	mov	rdi, rsi		;; Save source pointer
	cmp	ffttype, 2		;; Case off of FFT type
	je	type2c
	ja	type34c

	;; Type 1 FFT (forward FFT only)
	zr8_eight_complex_fft_final_preload
loop1c:	zr8_eight_complex_fft_final rsi, 8*128, 128, 2*128, 4*128
loop1r:	dec	rcx			;; Test loop counter
	jnz	loop1c
	pop	rax			;; Ignore return address
	zfft_1_ret

	;; Type 2 FFT (squaring)
type2c:	mov	r8, DEST2ARG
	zr8_eight_complex_with_square_preload
loop2c:	zr8_eight_complex_with_square rsi, 8*128, 128, 2*128, 4*128
loop2r:	dec	rcx			;; Test loop counter
	jnz	loop2c
	lea	rsi, [rdi+r8]		;; Restore source pointer
	mov	DESTARG, rsi
	ret

	;; FFT types 3 and 4
type34c:mov	rbp, DIST_TO_MULSRCARG
	cmp	ffttype, 4
	je	type4c

	;; Type 3 FFT (forward FFT and multiply)
	mov	r8, DEST2ARG
	zr8_eight_complex_with_mult_preload
loop3c:	zr8_eight_complex_with_mult rsi, 8*128, 128, 2*128, 4*128
loop3r:	dec	rcx			;; Test loop counter
	jnz	loop3c
	lea	rsi, [rdi+r8]		;; Restore source pointer
	mov	DESTARG, rsi
	ret

	;; Type 4 FFT (forward FFT and multiply)
type4c:	zr8_eight_complex_with_mulf_preload
loop4c:	zr8_eight_complex_with_mulf rsi, 8*128, 128, 2*128, 4*128
loop4r:	dec	rcx			;; Test loop counter
	jnz	loop4c
	mov	rsi, rdi		;; Restore source pointer
	ret
zcomplex_onepass_middle ENDP

;; Implement the small custom-coded one pass FFTs

buildfor SKX,	zonepass 128, 0
buildfor SKX,	zonepass 256, 0
buildfor SKX,	zonepass 384, 0
buildfor SKX,	zonepass 512, 0
buildfor SKX,	zonepass 640, 0
buildfor SKX,	zonepass 768, 0
buildfor SKX,	zonepass 896, 0
buildfor SKX,	zonepass 1K, 0
buildfor SKX,	zonepass 1152, 0
buildfor SKX,	zonepass 1280, 0
buildfor SKX,	zonepass 1536, 0
buildfor SKX,	zonepass 1792, 0
buildfor SKX,	zonepass 2K, 0
buildfor SKX,	zonepass 2304, 0
buildfor SKX,	zonepass 2560, 0
buildfor SKX,	zonepass 3K, 0
buildfor SKX,	zonepass 3584, 0
buildfor SKX,	zonepass 3840, 0
buildfor SKX,	zonepass 4K, 0
buildfor SKX,	zonepass 4608, 0
buildfor SKX,	zonepass 5K, 0
buildfor SKX,	zonepass 5376, 0
buildfor SKX,	zonepass 6K, 0
buildfor SKX,	zonepass 7K, 0
buildfor SKX,	zonepass 8K, 0
buildfor SKX,	zonepass 9K, 0
buildfor SKX,	zonepass 10K, 0
buildfor SKX,	zonepass 12K, 0
buildfor SKX,	zonepass 14K, 0
buildfor SKX,	zonepass 15K, 0
buildfor SKX,	zonepass 16K, 0
buildfor SKX,	zonepass 18K, 0
buildfor SKX,	zonepass 20K, 0
buildfor SKX,	zonepass 21K, 0
buildfor SKX,	zonepass 24K, 0
buildfor SKX,	zonepass 25K, 0
buildfor SKX,	zonepass 28K, 0
buildfor SKX,	zonepass 30K, 0
buildfor SKX,	zonepass 32K, 0
buildfor SKX,	zonepass 35K, 0
buildfor SKX,	zonepass 36K, 0

buildfor SKX,	zonepass 128, 1
buildfor SKX,	zonepass 256, 1
buildfor SKX,	zonepass 384, 1
buildfor SKX,	zonepass 512, 1
buildfor SKX,	zonepass 640, 1
buildfor SKX,	zonepass 768, 1
buildfor SKX,	zonepass 896, 1
buildfor SKX,	zonepass 1K, 1
buildfor SKX,	zonepass 1152, 1
buildfor SKX,	zonepass 1280, 1
buildfor SKX,	zonepass 1536, 1
buildfor SKX,	zonepass 1792, 1
buildfor SKX,	zonepass 2K, 1
buildfor SKX,	zonepass 2304, 1
buildfor SKX,	zonepass 2560, 1
buildfor SKX,	zonepass 3K, 1
buildfor SKX,	zonepass 3584, 1
buildfor SKX,	zonepass 3840, 1
buildfor SKX,	zonepass 4K, 1
buildfor SKX,	zonepass 4608, 1
buildfor SKX,	zonepass 5K, 1
buildfor SKX,	zonepass 5376, 1
buildfor SKX,	zonepass 6K, 1
buildfor SKX,	zonepass 7K, 1
buildfor SKX,	zonepass 8K, 1
buildfor SKX,	zonepass 9K, 1
buildfor SKX,	zonepass 10K, 1
buildfor SKX,	zonepass 12K, 1
buildfor SKX,	zonepass 14K, 1
buildfor SKX,	zonepass 15K, 1
buildfor SKX,	zonepass 16K, 1
buildfor SKX,	zonepass 18K, 1
buildfor SKX,	zonepass 20K, 1
buildfor SKX,	zonepass 21K, 1
buildfor SKX,	zonepass 24K, 1
buildfor SKX,	zonepass 25K, 1
buildfor SKX,	zonepass 28K, 1
buildfor SKX,	zonepass 30K, 1
buildfor SKX,	zonepass 32K, 1
buildfor SKX,	zonepass 35K, 1
buildfor SKX,	zonepass 36K, 1

_TEXT	ENDS
END
