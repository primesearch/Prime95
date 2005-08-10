; Copyright 2001-2005 Just For Fun Software, Inc., all rights reserved
; Author:  George Woltman
; Email: woltman@alum.mit.edu
;
; This routine implements part of a discrete weighted transform to
; quickly multiply two numbers.
;
; This code handles the last 8 levels of two pass FFTs that use the
; SSE2 instructions.
;

	TITLE   setup

IFNDEF X86_64
	.686
	.XMM
	.MODEL	FLAT
ENDIF

INCLUDE	unravel.mac
INCLUDE extrn.mac
INCLUDE xmult.mac
INCLUDE	xlucas.mac
INCLUDE xpass2.mac
INCLUDE xnormal.mac

_TEXT SEGMENT

PREFETCHING = 1

;; Routines to do the last 8 levels in a two-pass FFT

PUBLICP	xpass2_8_levels
PUBLICP	xpass2_8_levels_p

PROCP	xpass2_8_levels
	start_timer 2
	xpass2_8_levels_real 0
	add	rsi, pass1blkdst
	end_timer 2
LABELP	xpass2_8_levels_p
	mov	ecx, count1		; Number of complex iterations
	mov	rdx, pass2_premults	; Address of the group multipliers
p8lp:	xpass2_8_levels_complex 0
	add	rsi, pass1blkdst
	sub	ecx, 1
	jnz	p8lp
	mov	rsi, _DESTARG		; Restore source pointer
	ret
ENDPP	xpass2_8_levels

;;
;; Routines to do the normalization after a multiply
;;


PROCP _xnorm2

;; When doing zero-padded FFTs, the multiplied 7 words around the halfway point
;; must be subtracted from the bottom of the FFT.  This must be done before
;; normalization multiplies the FFT data by k.  This macro does that.

xsub_7_words MACRO
	LOCAL	nozpad, zlp
	IFDEF X86_64
	mov	rax, OFFSET XMM_ZPAD6
	cmp	zpad_addr, rax
	ELSE
	cmp	zpad_addr, OFFSET XMM_ZPAD6;; Have we subtracted all 7 words?
	ENDIF
	jg	short nozpad		;; Yes, skip this code
	mov	eax, cache_line_multiplier ;; Load loop counter
	push	rsi
	mov	rdi, zpad_addr		;; Addr of next zpad element to process
zlp:	_movsd	xmm0, [rsi]		;; Load FFT word
	_movsd	xmm1, [rdi]		;; Load ZPAD data
	mulsd	xmm1, XMM_NORM012_FF	;; Scale by FFTLEN/2
	subsd	xmm0, xmm1
	addsd	xmm7, xmm1		;; Adjust sumout
	_movsd	[rsi], xmm0
	lea	rsi, [rsi+64]		;; Bump pointers
	lea	rdi, [rdi+8]
	dec	eax			;; Iterate 2*clm (up to 8) times
	jnz	short zlp		;; Loop if necessary
	pop	rsi			;; Restore source ptr
	mov	zpad_addr, rdi
nozpad:
	ENDM

; Macro to loop through all the FFT values and apply the proper normalization
; routine.

inorm	MACRO	lab, ttp, zero, echk, const
	LOCAL	noadd, setlp, ilp0, ilp1, ilexit
	PUBLICP	lab
	LABELP	lab
zero	mov	zero_fft, 1		;; Set flag saying zero upper half
	movapd	xmm7, XMM_SUMOUT	;; Load SUMOUT
	movapd	xmm6, XMM_MAXERR	;; Load maximum error
no zero	cmp	edx, _ADDIN_ROW		;; Is this the time to do our addin?
no zero	jne	short noadd		;; Jump if addin does not occur now
no zero	mov	edi, _ADDIN_OFFSET	;; Get address to add value into
no zero	_movsd	xmm0, [rsi][rdi]	;; Get the value
no zero	addsd	xmm0, _ADDIN_VALUE	;; Add in the requested value
no zero	_movsd	[rsi][rdi], xmm0	;; Save the new value
no zero	subsd	xmm7, _ADDIN_VALUE	;; Do not include addin in sumout
noadd:	push	rdx
	push	rsi

	mov	rbx, norm_ptr2		;; Load column multipliers ptr
ttp	mov	eax, cache_line_multiplier ;; Load inner loop counter
	mov	rdi, OFFSET XMM_COL_MULTS ;; Load col mult scratch area
setlp:	xnorm_2d_setup ttp
ttp	lea	rdi, [rdi+512]		;; Next scratch area section
ttp	lea	rbx, [rbx+32]		;; Next column multiplier
ttp	sub	al, 1			;; Each cache line has its own col mult
ttp	jnz	setlp
ttp	mov	norm_ptr2, rbx		;; Save column multipliers ptr

	mov	rdx, norm_grp_mults	;; Addr of the group multipliers
	mov	rbp, carries		;; Addr of the carries
	mov	rdi, norm_ptr1		;; Load big/little flags array ptr
	mov	eax, addcount1		;; Load loop counter
	mov	loopcount2, eax		;; Save loop counter
	mov	loopcount3, 0		;; Clear outermost loop counter
	sub	rax, rax		;; Clear big/lit flags
	sub	rcx, rcx
ttp	mov	al, [rdi+0]		;; Load big vs. little flags
ttp	mov	cl, [rdi+1]		;; Load big vs. little flags
ilp0:	mov	ebx, cache_line_multiplier ;; Load inner loop counter
	mov	loopcount1, ebx		;; Save loop counter
	mov	rbx, OFFSET XMM_COL_MULTS ;; Load col mult scratch area
	xprefetcht1 [rdx+128]		;; Prefetch group multiplier
ilp1:	xprefetchw [rsi+64]
	xnorm_2d ttp, zero, echk, const	;; Normalize 8 values
	lea	rsi, [rsi+64]		;; Next cache line
ttp	lea	rbx, [rbx+512]		;; Next column multipliers
ttp	lea	rdi, [rdi+4]		;; Next big/little flags
	sub	loopcount1, 1		;; Test loop counter
	jnz	ilp1			;; Loop til done
	add	rsi, normblkdst		;; Skip gap in blkdst or clmblkdst
	lea	rbp, [rbp+64]		;; Next set of carries
ttp	lea	rdx, [rdx+128]		;; Next set of 8 group multipliers
	sub	loopcount2, 1		;; Test loop counter
	jz	ilexit			;; Jump when loop complete
	add	loopcount3, 80000000h/4 ;; 8 iterations
	jnc	ilp0
	add	rsi, normblkdst8	;; Add 128 every 8 clmblkdsts
	jmp	ilp0			;; Iterate
ilexit:	movapd	XMM_SUMOUT, xmm7	;; Save SUMOUT
	movapd	XMM_MAXERR, xmm6	;; Save maximum error
ttp	mov	norm_ptr1, rdi		;; Save big/little flags array ptr
	pop	rsi
	pop	rdx
	sub	rbx, rbx
	cmp	edx, 65536+256		;; Check for last iteration
	je	xtop_carry_adjust	;; Top carry may require adjusting
	ret
	ENDM

zpnorm	MACRO	lab, ttp, echk, const
	LOCAL	setlp, ilp0, ilp1, ilexit
	PUBLICP	lab
	LABELP	lab
const	mov	const_fft, 1		;; Set flag saying mul-by-const
	movapd	xmm7, XMM_SUMOUT	;; Load SUMOUT
	movapd	xmm6, XMM_MAXERR	;; Load maximum error
	xsub_7_words
	push	rdx
	push	rsi

	mov	rbx, norm_ptr2		;; Load column multipliers ptr
ttp	mov	eax, cache_line_multiplier ;; Load inner loop counter
	mov	rdi, OFFSET XMM_COL_MULTS ;; Load col mult scratch area
setlp:	xnorm_2d_setup ttp
ttp	lea	rdi, [rdi+512]		;; Next scratch area section
ttp	lea	rbx, [rbx+32]		;; Next column multiplier
ttp	sub	al, 1			;; Each cache line has its own col mult
ttp	jnz	setlp
ttp	mov	norm_ptr2, rbx		;; Save column multipliers ptr

	mov	rdx, norm_grp_mults	;; Addr of the group multipliers
	mov	rbp, carries		;; Addr of the carries
	mov	rdi, norm_ptr1		;; Load big/little flags array ptr
	mov	eax, addcount1		;; Load loop counter
	mov	loopcount2, eax		;; Save loop counter
	mov	loopcount3, 0		;; Clear outermost loop counter
	sub	rax, rax		;; Clear big/lit flags
ttp	mov	al, [rdi+0]		;; Load big vs. little flags
ilp0:	mov	ebx, cache_line_multiplier ;; Load inner loop counter
	mov	loopcount1, ebx		;; Save loop counter
	mov	rbx, OFFSET XMM_COL_MULTS ;; Load col mult scratch area
	xprefetcht1 [rdx+128]		;; Prefetch group multiplier
ilp1:	xprefetchw [rsi+64]
	xnorm_2d_zpad ttp, echk, const	;; Normalize 8 values
	lea	rsi, [rsi+64]		;; Next cache line
ttp	lea	rbx, [rbx+512]		;; Next column multipliers
ttp	lea	rdi, [rdi+4]		;; Next big/little flags
	sub	loopcount1, 1		;; Test loop counter
	jnz	ilp1			;; Loop til done
	add	rsi, normblkdst		;; Skip gap in blkdst or clmblkdst
	lea	rbp, [rbp+64]		;; Next set of carries
ttp	lea	rdx, [rdx+128]		;; Next set of 8 group multipliers
	sub	loopcount2, 1		;; Test loop counter
	jz	ilexit			;; Jump when loop complete
	add	loopcount3, 80000000h/4 ;; 8 iterations
	jnc	ilp0
	add	rsi, normblkdst8	;; Add 128 every 8 clmblkdsts
	jmp	ilp0			;; Iterate
ilexit:	movapd	XMM_SUMOUT, xmm7	;; Save SUMOUT
	movapd	XMM_MAXERR, xmm6	;; Save maximum error
ttp	mov	norm_ptr1, rdi		;; Save big/little flags array ptr
	pop	rsi
	pop	rdx
	sub	rbx, rbx
	ret
	ENDM

; The 16 different normalization routines.  One for each combination of
; rational/irrational, zeroing/no zeroing, error check/no error check, and
; mul by const/no mul by const.

	inorm	xr2, noexec, noexec, noexec, noexec
	inorm	xr2e, noexec, noexec, exec, noexec
	inorm	xr2c, noexec, noexec, noexec, exec
	inorm	xr2ec, noexec, noexec, exec, exec
	inorm	xr2z, noexec, exec, noexec, noexec
	inorm	xr2ze, noexec, exec, exec, noexec
	inorm	xi2, exec, noexec, noexec, noexec
	inorm	xi2e, exec, noexec, exec, noexec
	inorm	xi2c, exec, noexec, noexec, exec
	inorm	xi2ec, exec, noexec, exec, exec
	inorm	xi2z, exec, exec, noexec, noexec
	inorm	xi2ze, exec, exec, exec, noexec

	zpnorm	xr2zp, noexec, noexec, noexec
	zpnorm	xr2zpe, noexec, exec, noexec
	zpnorm	xr2zpc, noexec, noexec, exec
	zpnorm	xr2zpec, noexec, exec, exec
	zpnorm	xi2zp, exec, noexec, noexec
	zpnorm	xi2zpe, exec, exec, noexec
	zpnorm	xi2zpc, exec, noexec, exec
	zpnorm	xi2zpec, exec, exec, exec

; Special code to handle adjusting the carry out of the topmost FFT word

xtop_carry_adjust:
	xnorm_top_carry
	ret

ENDPP _xnorm2

_TEXT	ENDS
END
