; Copyright 2001 Just For Fun Software, Inc., all rights reserved
; Author:  George Woltman
; Email: woltman@alum.mit.edu
;

	TITLE   setup

	.686
	.XMM

_TEXT32 SEGMENT PARA USE32 PUBLIC 'DATA'

	ASSUME  CS: _TEXT32, DS: _TEXT32, SS: _TEXT32, ES: _TEXT32

INCLUDE extrn.mac
INCLUDE	unravel.mac
INCLUDE xmult.mac
INCLUDE memory.mac
INCLUDE xnormal.mac

;;
;; Add two FFTed numbers
;;

	PUBLIC	gwxaddf1
gwxaddf1 PROC NEAR
	mov	ecx, _SRCARG		; Address of first number
	mov	edx, _SRC2ARG		; Address of second number
	push	esi			; Save registers
	push	edi			; Save registers
	mov	eax, [ecx-4]		; Load needs-normalize counter
	mov	edi, [edx-4]		; Load needs-normalize counter
	mov	esi, _DESTARG		; Address of destination
	lea	eax, [eax+edi+1]	; Set new needs-normalize counter
	mov	[esi-4], eax		; Store needs-normalize counter
	push	ebx
	mov	bl, 2			; Load loop counter
	mov	eax, addcount1		; Load loop counter
faddlp:	movapd	xmm0, [edx]		; Load second number
	addpd	xmm0, [ecx]		; Add in first number
	movapd	xmm1, [edx+16]		; Load second number
	addpd	xmm1, [ecx+16]		; Add in first number
	movapd	xmm2, [edx+32]		; Load second number
	addpd	xmm2, [ecx+32]		; Add in first number
	movapd	xmm3, [edx+48]		; Load second number
	addpd	xmm3, [ecx+48]		; Add in first number
	movapd	[esi], xmm0		; Save result
	movapd	[esi+16], xmm1		; Save result
	movapd	[esi+32], xmm2		; Save result
	movapd	[esi+48], xmm3		; Save result
	lea	ecx, [ecx+64]		; Next source
	lea	edx, [edx+64]		; Next source
	lea	esi, [esi+64]		; Next dest
	sub	ax, 1			; Check loop counter
	jnz	short faddlp		; Loop if necessary
	shr	eax, 16			; Optionally bump pointers
	add	ecx, eax
	add	edx, eax
	add	esi, eax
	mov	eax, addcount2
	sub	bl, 1			; Test loop counter
	jnz	short faddlp
	pop	ebx
	pop	edi			; Restore registers
	pop	esi			; Restore registers
	ret
gwxaddf1 ENDP


;;
;; Subtract two FFTed numbers
;;

	PUBLIC	gwxsubf1
gwxsubf1 PROC NEAR
	mov	ecx, _SRCARG		; Address of first number
	mov	edx, _SRC2ARG		; Address of second number
	push	esi			; Save registers
	push	edi			; Save registers
	mov	eax, [ecx-4]		; Load needs-normalize counter
	mov	edi, [edx-4]		; Load needs-normalize counter
	mov	esi, _DESTARG		; Address of destination
	lea	eax, [eax+edi+1]	; Set new needs-normalize counter
	mov	[esi-4], eax		; Store needs-normalize counter
	push	ebx
	mov	bl, 2			; Load loop counter
	mov	eax, addcount1		; Load loop counter
fsublp:	movapd	xmm0, [edx]		; Load second number
	subpd	xmm0, [ecx]		; Subtract first number
	movapd	xmm1, [edx+16]		; Load second number
	subpd	xmm1, [ecx+16]		; Subtract first number
	movapd	xmm2, [edx+32]		; Load second number
	subpd	xmm2, [ecx+32]		; Subtract first number
	movapd	xmm3, [edx+48]		; Load second number
	subpd	xmm3, [ecx+48]		; Subtract first number
	movapd	[esi], xmm0		; Save result
	movapd	[esi+16], xmm1		; Save result
	movapd	[esi+32], xmm2		; Save result
	movapd	[esi+48], xmm3		; Save result
	lea	ecx, [ecx+64]		; Next source
	lea	edx, [edx+64]		; Next source
	lea	esi, [esi+64]		; Next dest
	sub	ax, 1			; Check loop counter
	jnz	short fsublp		; Loop if necessary
	shr	eax, 16			; Optionally bump pointers
	add	ecx, eax
	add	edx, eax
	add	esi, eax
	mov	eax, addcount2
	sub	bl, 1			; Test loop counter
	jnz	short fsublp
	pop	ebx
	pop	edi			; Restore registers
	pop	esi			; Restore registers
	ret
gwxsubf1 ENDP


;;
;; Add and subtract two FFTed numbers
;;

	PUBLIC	gwxaddsubf1
gwxaddsubf1 PROC NEAR
	mov	ecx, _SRCARG		; Address of first number
	mov	edx, _SRC2ARG		; Address of second number
	push	esi			; Save registers
	push	edi			; Save registers
	mov	eax, [ecx-4]		; Load needs-normalize counter
	mov	edi, [edx-4]		; Load needs-normalize counter
	push	ebp			; Save registers
	mov	esi, _DESTARG		; Address of destination #1
	mov	ebp, _DEST2ARG  	; Address of destination #2
	lea	eax, [eax+edi+1]	; Set new needs-normalize counter
	mov	[esi-4], eax		; Store needs-normalize counter
	mov	[ebp-4], eax		; Store needs-normalize counter
	push	ebx
	mov	bl, 2			; Load loop counter
	mov	eax, addcount1		; Load loop counter
faddsublp:
	movapd	xmm0, [ecx]		; Load first number
	movapd	xmm1, xmm0		; Dup first number
	addpd	xmm0, [edx]		; Add in second number
	subpd	xmm1, [edx]		; Subtract out second number
	movapd	xmm2, [ecx+16]		; Load first number
	movapd	xmm3, xmm2		; Dup first number
	addpd	xmm2, [edx+16]		; Add in second number
	subpd	xmm3, [edx+16]		; Subtract out second number
	movapd	xmm4, [ecx+32]		; Load first number
	movapd	xmm5, xmm4		; Dup first number
	addpd	xmm4, [edx+32]		; Add in second number
	subpd	xmm5, [edx+32]		; Subtract out second number
	movapd	xmm6, [ecx+48]		; Load first number
	movapd	xmm7, xmm6		; Dup first number
	addpd	xmm6, [edx+48]		; Add in second number
	subpd	xmm7, [edx+48]		; Subtract out second number
	movapd	[esi], xmm0		; Save result
	movapd	[ebp], xmm1		; Save result
	movapd	[esi+16], xmm2		; Save result
	movapd	[ebp+16], xmm3		; Save result
	movapd	[esi+32], xmm4		; Save result
	movapd	[ebp+32], xmm5		; Save result
	movapd	[esi+48], xmm6		; Save result
	movapd	[ebp+48], xmm7		; Save result
	lea	ecx, [ecx+64]		; Next source
	lea	edx, [edx+64]		; Next source
	lea	esi, [esi+64]		; Next dest
	lea	ebp, [ebp+64]		; Next dest
	sub	ax, 1			; Check loop counter
	JNZ_X	faddsublp		; Loop if necessary
	shr	eax, 16			; Optionally bump pointers
	add	ecx, eax
	add	edx, eax
	add	esi, eax
	add	ebp, eax
	mov	eax, addcount2
	sub	bl, 1			; Test loop counter
	JNZ_X	faddsublp
	pop	ebx
	pop	ebp			; Restore registers
	pop	edi			; Restore registers
	pop	esi			; Restore registers
	ret
gwxaddsubf1 ENDP



;;
;; Add two numbers without carry propogation.  Caller can use this for
;; consecutive add or subtract operations.  However, the last operation
;; before a multiply must use the routine that will normalize data.
;;

	PUBLIC	gwxaddq1
gwxaddq1 PROC NEAR
	mov	ecx, _SRCARG		; Address of first number
	mov	edx, _SRC2ARG		; Address of second number
	push	esi			; Save registers
	push	edi			; Save registers
	mov	eax, [ecx-4]		; Load needs-normalize counter
	mov	edi, [edx-4]		; Load needs-normalize counter
	mov	esi, _DESTARG		; Address of destination
	lea	eax, [eax+edi+1]	; Set new needs-normalize counter
	mov	[esi-4], eax		; Store needs-normalize counter
	jmp	short uadd		; Do the unnormalized add
gwxaddq1 ENDP


;;
;; Add two numbers.  This will be done without carry propogation
;; (for extra speed) if it can be done safely.
;;

	PUBLIC	gwxadd1
gwxadd1	PROC NEAR

; See if add can be done without normalization

	mov	ecx, _SRCARG		; Address of first number
	mov	edx, _SRC2ARG		; Address of second number
	push	esi			; Save registers
	push	edi			; Save registers
	mov	eax, [ecx-4]		; Load needs-normalize counter
	mov	edi, [edx-4]		; Load needs-normalize counter
	mov	esi, _DESTARG		; Address of destination
	lea	eax, [eax+edi+1]	; Set new needs-normalize counter
	mov	[esi-4], eax		; Store needs-normalize counter
	cmp	eax, extra_bits		; Is normalization needed?
	jg	short nadd		; Yes, do a normalized add

; Do an unnormalized add
; ecx = src arg #1
; edx = src arg #2
; esi = dest arg

uadd:	mov	eax, normcount1		; Load loop counter
uaddlp:	movapd	xmm0, [edx]		; Load second number
	addpd	xmm0, [ecx]		; Add in first number
	movapd	xmm1, [edx+16]		; Load second number
	addpd	xmm1, [ecx+16]		; Add in first number
	movapd	xmm2, [edx+32]		; Load second number
	addpd	xmm2, [ecx+32]		; Add in first number
	movapd	xmm3, [edx+48]		; Load second number
	addpd	xmm3, [ecx+48]		; Add in first number
	movapd	[esi], xmm0		; Save result
	movapd	[esi+16], xmm1		; Save result
	movapd	[esi+32], xmm2		; Save result
	movapd	[esi+48], xmm3		; Save result
	lea	ecx, [ecx+64]		; Next source
	lea	edx, [edx+64]		; Next source
	lea	esi, [esi+64]		; Next dest
	sub	ax, 1			; Check loop counter
	jnz	short uaddlp		; Loop if necessary
uaddlp2:sub	eax, 65536		; Check loop counter
	js	short uadddn		; Loop if necessary
	movapd	xmm0, [edx]		; Load second number
	addpd	xmm0, [ecx]		; Add in first number
	movapd	xmm1, [edx+16]		; Load second number
	addpd	xmm1, [ecx+16]		; Add in first number
	movapd	[esi], xmm0		; Save result
	movapd	[esi+16], xmm1		; Save result
	lea	ecx, [ecx+64]		; Next source
	lea	edx, [edx+64]		; Next source
	lea	esi, [esi+64]		; Next dest
	jmp	short uaddlp2		; Loop
uadddn:	pop	edi			; Restore registers
	pop	esi			; Restore registers
	ret

; Do a normalized add
; ecx = src arg #1
; edx = src arg #2
; esi = dest arg

nadd:	push	ebp			; Save registers
	push	ebx			; Save registers
	movapd	xmm2, XMM_BIGVAL	; Start process with no carry
	movapd	xmm3, xmm2
	movapd	xmm6, xmm2
	movapd	xmm7, xmm2
	mov	eax, normcount1		; Load loop counter
	mov	loopcount1, eax		; Save loop counter
	mov	ebp, norm_col_mults	; Address of the multipliers
	mov	edi, norm_biglit_array	; Addr of the big/little flags array
	mov	eax, 48			; Rational FFTs are all big words
	mov	ebx, 48
	cmp	_NUMLIT, 0		; Test for irrational FFTs
	JNE_X	naddlp			; Yes, use two-to-phi multipliers
raddlp:	xnorm_op_1d addpd, 8, noexec	; Add and normalize 8 values
	sub	WORD PTR loopcount1, 1	; Decrement loop counter
	JNZ_X	raddlp			; Loop til done
raddlp2:sub	loopcount1, 65536	; Test loop counter
	JS_X	nadddn			; Jump if done
	xnorm_op_1d addpd, 4, noexec	; Add and normalize 4 values
	JMP_X	raddlp2			; Loop til done
naddlp:	xnorm_op_1d addpd, 8, exec	; Add and normalize 8 values
	sub	WORD PTR loopcount1, 1	; Decrement loop counter
	JNZ_X	naddlp			; Loop til done
naddlp2:sub	loopcount1, 65536	; Test loop counter
	JS_X	nadddn			; Jump if done
	xnorm_op_1d addpd, 4, exec	; Add and normalize 4 values
	JMP_X	naddlp2			; Loop til done
nadddn:	mov	ecx, _DESTARG		; Address of result
	mov	ebp, norm_col_mults	; Address of the multipliers
	xnorm_op_1d_cleanup
	mov	DWORD PTR [ecx-4], 0	; Clear needs-normalize flag
	pop	ebx			; Restore registers
	pop	ebp			; Restore registers
	pop	edi			; Restore registers
	pop	esi			; Restore registers
	ret
gwxadd1	ENDP


;;
;; Subtract two numbers without carry propogation.  Caller can use this for
;; consecutive add or subtract operations.  However, the last operation
;; before a multiply must use the routine that will normalize data.
;;

	PUBLIC	gwxsubq1
gwxsubq1	PROC NEAR
	mov	ecx, _SRCARG		; Address of first number
	mov	edx, _SRC2ARG		; Address of second number
	push	esi			; Save registers
	push	edi			; Save registers
	mov	eax, [ecx-4]		; Load needs-normalize counter
	mov	edi, [edx-4]		; Load needs-normalize counter
	mov	esi, _DESTARG		; Address of destination
	lea	eax, [eax+edi+1]	; Set new needs-normalize counter
	mov	[esi-4], eax		; Store needs-normalize counter
	jmp	short usub		; Do the unnormalized sub
gwxsubq1	ENDP

;;
;; Subtract two numbers.  This will be done without carry propogation
;; (for extra speed) if it can be done safely.
;;

	PUBLIC	gwxsub1
gwxsub1	PROC NEAR

; See if subtract can be done without normalization

	mov	ecx, _SRCARG		; Address of first number
	mov	edx, _SRC2ARG		; Address of second number
	push	esi			; Save registers
	push	edi			; Save registers
	mov	eax, [ecx-4]		; Load needs-normalize counter
	mov	edi, [edx-4]		; Load needs-normalize counter
	mov	esi, _DESTARG		; Address of destination
	lea	eax, [eax+edi+1]	; Set new needs-normalize counter
	mov	[esi-4], eax		; Store needs-normalize counter
	cmp	eax, extra_bits		; Is normalization needed?
	jg	short nsub		; Yes, do a normalized subtract

; Do an unnormalized subtract
; ecx = src arg #1
; edx = src arg #2
; esi = dest arg

usub:	mov	eax, normcount1		; Load loop counter
usublp:	movapd	xmm0, [edx]		; Load second number
	subpd	xmm0, [ecx]		; Subtract first number
	movapd	xmm1, [edx+16]		; Load second number
	subpd	xmm1, [ecx+16]		; Subtract first number
	movapd	xmm2, [edx+32]		; Load second number
	subpd	xmm2, [ecx+32]		; Subtract first number
	movapd	xmm3, [edx+48]		; Load second number
	subpd	xmm3, [ecx+48]		; Subtract first number
	movapd	[esi], xmm0		; Save result
	movapd	[esi+16], xmm1		; Save result
	movapd	[esi+32], xmm2		; Save result
	movapd	[esi+48], xmm3		; Save result
	lea	ecx, [ecx+64]		; Next source
	lea	edx, [edx+64]		; Next source
	lea	esi, [esi+64]		; Next dest
	sub	ax, 1			; Check loop counter
	jnz	short usublp		; Loop if necessary
usublp2:sub	eax, 65536		; Check loop counter
	js	short usubdn		; Loop if necessary
	movapd	xmm0, [edx]		; Load second number
	subpd	xmm0, [ecx]		; Subtract first number
	movapd	xmm1, [edx+16]		; Load second number
	subpd	xmm1, [ecx+16]		; Subtract first number
	movapd	[esi], xmm0		; Save result
	movapd	[esi+16], xmm1		; Save result
	lea	ecx, [ecx+64]		; Next source
	lea	edx, [edx+64]		; Next source
	lea	esi, [esi+64]		; Next dest
	jmp	short usublp2		; Loop
usubdn:	pop	edi			; Restore registers
	pop	esi			; Restore registers
	ret

; Do a normalized subtract
; ecx = src arg #1
; edx = src arg #2
; esi = dest arg

nsub:	push	ebp			; Save registers
	push	ebx			; Save registers
	movapd	xmm2, XMM_BIGVAL	; Start process with no carry
	movapd	xmm3, xmm2
	mov	eax, normcount1		; Load loop counter
	mov	loopcount1, eax		; Save loop counter
	mov	ebp, norm_col_mults	; Address of the multipliers
	mov	edi, norm_biglit_array	; Addr of the big/little flags array
	mov	eax, 48			; Rational FFTs are all big words
	mov	ebx, 48
	cmp	_NUMLIT, 0		; Test for irrational FFTs
	JNE_X	nsublp			; Yes, use two-to-phi multipliers
rsublp:	xnorm_op_1d subpd, 8, noexec	; Subtract and normalize 8 values
	sub	WORD PTR loopcount1, 1	; Decrement loop counter
	JNZ_X	rsublp			; Loop til done
rsublp2:sub	loopcount1, 65536	; Test loop counter
	JS_X	nsubdn			; Jump if done
	xnorm_op_1d subpd, 4, noexec	; Subtract and normalize 4 values
	JMP_X	rsublp2			; Loop til done
nsublp:	xnorm_op_1d subpd, 8, exec	; Subtract and normalize 8 values
	sub	WORD PTR loopcount1, 1	; Decrement loop counter
	JNZ_X	nsublp			; Loop til done
nsublp2:sub	loopcount1, 65536	; Test loop counter
	JS_X	nsubdn			; Jump if done
	xnorm_op_1d subpd, 4, exec	; Subtract and normalize 4 values
	JMP_X	nsublp2			; Loop til done
nsubdn:	mov	ecx, _DESTARG		; Address of result
	mov	ebp, norm_col_mults	; Address of the multipliers
	xnorm_op_1d_cleanup
	mov	DWORD PTR [ecx-4], 0	; Clear needs-normalize flag
	pop	ebx
	pop	ebp
	pop	edi
	pop	esi
	ret
gwxsub1	ENDP


;;
;; Add and subtract two numbers without carry propogation.
;;

	PUBLIC	gwxaddsubq1
gwxaddsubq1 PROC NEAR
	mov	ecx, _SRCARG		; Address of first number
	mov	edx, _SRC2ARG		; Address of second number
	push	esi			; Save registers
	push	edi			; Save registers
	mov	eax, [ecx-4]		; Load needs-normalize counter
	mov	edi, [edx-4]		; Load needs-normalize counter
	push	ebp			; Save registers
	mov	esi, _DESTARG		; Address of destination #1
	mov	ebp, _DEST2ARG  	; Address of destination #2
	lea	eax, [eax+edi+1]	; Set new needs-normalize counter
	mov	[esi-4], eax		; Store needs-normalize counter
	mov	[ebp-4], eax		; Store needs-normalize counter
	jmp	short uaddsub		; Do the unnormalized add
gwxaddsubq1 ENDP

;;
;; Add and subtract two numbers.  This will be done without carry propogation
;; (for extra speed) if it can be done safely.
;;

	PUBLIC	gwxaddsub1
gwxaddsub1 PROC NEAR

; See if add can be done without normalization

	mov	ecx, _SRCARG		; Address of first number
	mov	edx, _SRC2ARG		; Address of second number
	push	esi			; Save registers
	push	edi			; Save registers
	mov	eax, [ecx-4]		; Load needs-normalize counter
	mov	edi, [edx-4]		; Load needs-normalize counter
	push	ebp			; Save registers
	mov	esi, _DESTARG		; Address of destination
	mov	ebp, _DEST2ARG  	; Address of destination #2
	lea	eax, [eax+edi+1]	; Set new needs-normalize counter
	mov	[esi-4], eax		; Store needs-normalize counter
	mov	[ebp-4], eax		; Store needs-normalize counter
	cmp	eax, extra_bits		; Is normalization needed?
	JG_X	naddsub			; Yes, do a normalized add

; Do an unnormalized add
; ecx = src arg #1
; edx = src arg #2
; esi = dest arg #1
; ebp = dest arg #2

uaddsub:mov	eax, normcount1		; Load loop counter
uaddsublp:
	movapd	xmm0, [ecx]		; Load first number
	movapd	xmm1, xmm0		; Dup first number
	addpd	xmm0, [edx]		; Add in second number
	subpd	xmm1, [edx]		; Subtract out second number
	movapd	xmm2, [ecx+16]		; Load first number
	movapd	xmm3, xmm2		; Dup first number
	addpd	xmm2, [edx+16]		; Add in second number
	subpd	xmm3, [edx+16]		; Subtract out second number
	movapd	xmm4, [ecx+32]		; Load first number
	movapd	xmm5, xmm4		; Dup first number
	addpd	xmm4, [edx+32]		; Add in second number
	subpd	xmm5, [edx+32]		; Subtract out second number
	movapd	xmm6, [ecx+48]		; Load first number
	movapd	xmm7, xmm6		; Dup first number
	addpd	xmm6, [edx+48]		; Add in second number
	subpd	xmm7, [edx+48]		; Subtract out second number
	movapd	[esi], xmm0		; Save result
	movapd	[ebp], xmm1		; Save result
	movapd	[esi+16], xmm2		; Save result
	movapd	[ebp+16], xmm3		; Save result
	movapd	[esi+32], xmm4		; Save result
	movapd	[ebp+32], xmm5		; Save result
	movapd	[esi+48], xmm6		; Save result
	movapd	[ebp+48], xmm7		; Save result
	lea	ecx, [ecx+64]		; Next source
	lea	edx, [edx+64]		; Next source
	lea	esi, [esi+64]		; Next dest
	lea	ebp, [ebp+64]		; Next dest
	sub	ax, 1			; Check loop counter
	JNZ_X	uaddsublp		; Loop if necessary
uaddsublp2:
	sub	eax, 65536		; Check loop counter
	js	short uaddsubdn		; Loop if necessary
	movapd	xmm0, [ecx]		; Load first number
	movapd	xmm1, xmm0		; Dup first number
	addpd	xmm0, [edx]		; Add in second number
	subpd	xmm1, [edx]		; Subtract out second number
	movapd	xmm2, [ecx+16]		; Load first number
	movapd	xmm3, xmm2		; Dup first number
	addpd	xmm2, [edx+16]		; Add in second number
	subpd	xmm3, [edx+16]		; Subtract out second number
	movapd	[esi], xmm0		; Save result
	movapd	[ebp], xmm1		; Save result
	movapd	[esi+16], xmm2		; Save result
	movapd	[ebp+16], xmm3		; Save result
	lea	ecx, [ecx+64]		; Next source
	lea	edx, [edx+64]		; Next source
	lea	esi, [esi+64]		; Next dest
	lea	ebp, [ebp+64]		; Next dest
	jmp	short uaddsublp2	; Loop
uaddsubdn:
	pop	ebp			; Restore registers
	pop	edi			; Restore registers
	pop	esi			; Restore registers
	ret

; Do a normalized add and subtract
; ecx = src arg #1
; edx = src arg #2
; esi = dest arg #1
; ebp = dest arg #2

naddsub:push	ebx			; Save registers
	movapd	xmm2, XMM_BIGVAL	; Start process with no carry
	movapd	xmm3, xmm2
	mov	eax, normcount1		; Load loop counter
	mov	loopcount1, eax		; Save loop counter
	mov	ebx, norm_col_mults	; Address of the multipliers
	mov	edi, norm_biglit_array	; Addr of the big/little flags array
	mov	eax, 48			; Rational FFTs are all big words
	cmp	_NUMLIT, 0		; Test for irrational FFTs
	JNE_X	naddsublp		; Yes, use two-to-phi multipliers
raddsublp:
	xnorm_addsub_1d 8, noexec	; Add/sub and normalize 8 values
	sub	WORD PTR loopcount1, 1	; Decrement loop counter
	JNZ_X	raddsublp		; Loop til done
raddsublp2:
	sub	loopcount1, 65536	; Test loop counter
	JS_X	naddsubdn		; Jump if done
	xnorm_addsub_1d 4, noexec	; Add/sub and normalize 4 values
	JMP_X	raddsublp2		; Loop til done
naddsublp:
	xnorm_addsub_1d 8, exec		; Add/sub and normalize 8 values
	sub	WORD PTR loopcount1, 1	; Decrement loop counter
	JNZ_X	naddsublp		; Loop til done
naddsublp2:
	sub	loopcount1, 65536	; Test loop counter
	JS_X	naddsubdn		; Jump if done
	xnorm_addsub_1d 4, exec		; Add/sub and normalize 4 values
	JMP_X	naddsublp2		; Loop til done
naddsubdn:
	mov	esi, _DESTARG		; Address of result #1
	mov	ebp, _DEST2ARG		; Address of result #2
	mov	ebx, norm_col_mults	; Address of the multipliers
	xnorm_addsub_1d_cleanup
	mov	DWORD PTR [esi-4], 0	; Clear needs-normalize flag
	mov	DWORD PTR [ebp-4], 0	; Clear needs-normalize flag
	pop	ebx			; Restore registers
	pop	ebp			; Restore registers
	pop	edi			; Restore registers
	pop	esi			; Restore registers
	ret
gwxaddsub1 ENDP

;;
;; Copy one number and zero some low order words.
;;

	PUBLIC	gwxcopyzero1
gwxcopyzero1 PROC NEAR
	push	esi			; Save registers
	push	edi			; Save registers
	mov	esi, _SRCARG		; Address of first number
	mov	edi, _DESTARG		; Address of destination
	mov	eax, [esi-4]		; Load needs-normalize counter
	mov	[edi-4], eax		; Store needs-normalize counter
	sub	ecx, ecx		; Offset to compare to COPYZERO
	mov	eax, normcount1		; Load loop counter
cz1:	xcopyzero 8			; Copy/zero 8 values
	lea	esi, [esi+64]		; Next source
	lea	edi, [edi+64]		; Next dest
	lea	ecx, [ecx+64]		; Next compare offset
	sub	ax, 1			; Test loop counter
	JNZ_X	cz1			; Loop if necessary
cz2:	sub	eax, 65536		; Test loop counter
	js	short cz3		; Jump if done
	xcopyzero 4			; Copy/zero 8 values
	lea	esi, [esi+64]		; Next source
	lea	edi, [edi+64]		; Next dest
	lea	ecx, [ecx+64]		; Next compare offset
	jmp	short cz2		; Loop if necessary
cz3:	pop	edi			; Restore registers
	pop	esi			; Restore registers
	ret
gwxcopyzero1 ENDP

;;
;; Do a mod k*2^n+/-1
;;

	PUBLIC	gwxprothmod1
gwxprothmod1 PROC NEAR
	push	esi			; Save registers
	push	edi
	push	ebp
	push	ebx
	mov	esi, _SRCARG		; Address of first number
	prothmod_upper_prep_0d		; Prepare for prothmod stage 1
	sub	eax, eax		; Compute total num of cache lines
	mov	ax, WORD PTR normcount1
	mov	ebx, eax		; Compute # bytes in these cache lines
	shl	ebx, 6
	mov	ebp, eax		; Load loop counter
pmulp:	prothmod_upper_0d 16, 64, 2*64, pmudn ; Divide 8 upper values by k
	lea	esi, [esi-4*64]		; Next source word
	lea	edi, [edi+8*8]		; Next scratch words
	sub	ebp, 4			; Test loop counter
	JNZ_X	pmulp			; Iterate
	mov	ebp, eax		; Load loop counter
	add	esi, ebx		; Restore source pointer
	sub	esi, 8			; Next source pointer
	add	eax, 80000000h		; Test loop counter
	JNC_X	pmulp
	sub	esi, ebx
	sub	eax, eax		; Compute total num of cache lines
	mov	ax, WORD PTR normcount1
	add	ax, WORD PTR normcount1+2
	mov	ebx, eax		; Compute # bytes in these cache lines
	shl	ebx, 6
	mov	ebp, eax		; Load loop counter
	add	esi, ebx		; Restore source pointer
	sub	esi, 16			; Next source pointer
	JMP_X	pmulp			; Loop til done
pmudn:	prothmod_lower_prep_0d		; Prepare for prothmod stage 2
	mov	ebp, eax		; Load loop counter
pmllp:	prothmod_lower_0d 16, 64, 2*64	; Add 8 shifted quotients to lower data
	lea	esi, [esi+4*64]		; Next source word
	lea	edi, [edi-8*8]		; Next scratch words
	dec	edx
	jz	short pmldn
	sub	ebp, 4			; Test loop counter
	JNZ_X	pmllp			; Iterate
	mov	ebp, eax		; Load loop counter
	sub	esi, ebx		; Restore source pointer
	add	esi, 8			; Next source pointer
	add	eax, 80000000h		; Test loop counter
	JNC_X	pmllp
	mov	ebp, eax		; Load loop counter
	add	esi, 16			; Next source pointer
	JMP_X	pmllp			; Loop til done
pmldn:	prothmod_final_0d 16, 64, 2*64	; Finish off the mod
	pop	ebx			; Restore registers
	pop	ebp
	pop	edi
	pop	esi
	ret
gwxprothmod1 ENDP

;;
;; Routines to do the normalization after a multiply
;;

_xsimpnorm PROC NEAR

; Macro to loop through all the FFT values and apply the proper normalization
; routine.

inorm	MACRO	lab, ttp, zero, echk, const
	LOCAL	ilp1, ilp2
	PUBLIC	lab
lab:	mov	esi, _DESTARG		;; Addr of multiplied number
zero	mov	zero_fft, 1		; Set flag saying we've initialized
zero	add	_FFTZERO, esi
zero	add	_FFTZERO+4, esi
zero	add	_FFTZERO+8, esi
zero	add	_FFTZERO+12, esi
zero	add	_FFTZERO+16, esi
zero	add	_FFTZERO+20, esi
zero	add	_FFTZERO+24, esi
zero	add	_FFTZERO+28, esi
	subpd	xmm7, xmm7		;; Init SUMOUT
	mov	edi, _ADDIN_OFFSET	;; Get address to add value into
	_movsd	xmm0, [esi][edi]	;; Get the value
	addsd	xmm0, _ADDIN_VALUE	;; Add in the requested value
	_movsd	[esi][edi], xmm0	;; Save the new value
	subsd	xmm7, _ADDIN_VALUE	;; Do not include addin in sumout
	movapd	xmm2, XMM_BIGVAL	;; Start process with no carry
	movapd	xmm3, xmm2
	movlpd	xmm6, _MAXERR		;; Current maximum error
	movhpd	xmm6, _MAXERR
	mov	ebp, norm_col_mults	;; Addr of the multipliers
	mov	edi, norm_biglit_array	;; Addr of the big/little flags array
	mov	eax, 48			;; Rational FFTs are all big words
	mov	ecx, 48
	mov	ebx, normcount1		;; Load loop counter
ilp1:	xnorm_1d 8, ttp, zero, echk, const ;; Normalize 8 values
	lea	esi, [esi+64]		;; Next cache line
ttp	lea	ebp, [ebp+128]		;; Next set of 8 multipliers
ttp	lea	edi, [edi+4]		;; Next big/little flags
	sub	bx, 1			;; Test loop counter
	JNZ_X	ilp1			;; Loop til done
ilp2:	sub	ebx, 65536		;; Test loop counter
	JS_X	idn
	xnorm_1d 4, ttp, zero, echk, const ;; Normalize 4 values
	lea	esi, [esi+64]		;; Next cache line
ttp	lea	ebp, [ebp+64]		;; Next set of 4 multipliers
ttp	lea	edi, [edi+2]		;; Next big/little flags
	JMP_X	ilp2			;; Loop til done
	ENDM

; The 16 different normalization routines.  One for each combination of
; rational/irrational, zeroing/no zeroing, error check/no error check, and
; mul by const/no mul by const.

	inorm	xr1, noexec, noexec, noexec, noexec
	inorm	xr1e, noexec, noexec, exec, noexec
	inorm	xr1c, noexec, noexec, noexec, exec
	inorm	xr1ec, noexec, noexec, exec, exec
	inorm	xr1z, noexec, exec, noexec, noexec
	inorm	xr1ze, noexec, exec, exec, noexec
	inorm	xi1, exec, noexec, noexec, noexec
	inorm	xi1e, exec, noexec, exec, noexec
	inorm	xi1c, exec, noexec, noexec, exec
	inorm	xi1ec, exec, noexec, exec, exec
	inorm	xi1z, exec, exec, noexec, noexec
	inorm	xi1ze, exec, exec, exec, noexec

; Finish off the normalization process by add any carry to first values.
; Handle both the with and without two-to-phi array cases.

idn:	mov	esi, _DESTARG		; Address of squared number
	mov	edi, norm_biglit_array	; Address of the big/little flags array
	mov	ebp, norm_col_mults	; Restart the column multipliers
	cmp	zero_fft, 1		; Is this a zero high words fft?
	JE_X	idz			; Yes, do special add in of carries
	xnorm012_1d noexec		; Add in carries
	JMP_X	cmnend			; All done, go cleanup
idz:	xnorm012_1d exec		; Add in carries
	mov	zero_fft, 0		; Clear zero-high-words-fft flag
	sub	_FFTZERO, esi		; Restore FFTZERO addresses
	sub	_FFTZERO+4, esi
	sub	_FFTZERO+8, esi
	sub	_FFTZERO+12, esi
	sub	_FFTZERO+16, esi
	sub	_FFTZERO+20, esi
	sub	_FFTZERO+24, esi
	sub	_FFTZERO+28, esi

; Clear needs-normalize counter

cmnend:	mov	DWORD PTR [esi-4], 0

; Normalize SUMOUT value by multiplying by 1 / (fftlen/2).

	movapd	XMM_TMP1, xmm7		; Add together the two partial sumouts
	addsd	xmm7, XMM_TMP1+8
	mulsd	xmm7, ttmp_ff_inv
	_movsd	[esi-24], xmm7		; Save sum of FFT outputs
	movapd	XMM_TMP1, xmm6		; Compute new maximum error
	maxsd	xmm6, XMM_TMP1+8
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
_xsimpnorm ENDP

_TEXT32	ENDS
END
