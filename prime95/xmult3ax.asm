; Copyright 2001-2003 Just For Fun Software, Inc., all rights reserved
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

blkrows	EQU	1024

gwx3procs PROC NEAR

;;
;; Add two numbers without carry propogation.  Caller can use this for
;; consecutive add or subtract operations.  However, the last operation
;; before a multiply must use the routine that will normalize data.
;;

	PUBLIC	gwxaddf3
	PUBLIC	gwxaddq3
gwxaddf3:
gwxaddq3:
	mov	ecx, _SRCARG		; Address of first number
	mov	edx, _SRC2ARG		; Address of second number
	push	esi			; Save registers
	push	edi			; Save registers
	mov	eax, [ecx-4]		; Load needs-normalize counter
	mov	edi, [edx-4]		; Load needs-normalize counter
	mov	esi, _DESTARG		; Address of destination
	lea	eax, [eax+edi+1]	; Set new needs-normalize counter
	mov	edi, [ecx-28]		; Load has-been-FFTed flag
	mov	[esi-4], eax		; Store needs-normalize counter
	mov	[esi-28], edi		; Save has-been-FFTed flag
	jmp	short uadd		; Do the unnormalized add

;;
;; Add two numbers.  This will be done without carry propogation
;; (for extra speed) if it can be done safely.
;;

	PUBLIC	gwxadd3
gwxadd3:

; See if add can be done without normalization

	mov	ecx, _SRCARG		; Address of first number
	mov	edx, _SRC2ARG		; Address of second number
	push	esi			; Save registers
	push	edi			; Save registers
	mov	eax, [ecx-4]		; Load needs-normalize counter
	mov	edi, [edx-4]		; Load needs-normalize counter
	mov	esi, _DESTARG		; Address of destination
	lea	eax, [eax+edi+1]	; Set new needs-normalize counter
	mov	DWORD PTR [esi-28], 0	; Save has-been-FFTed flag
	cmp	eax, extra_bits		; Is normalization needed?
	jg	nadd			; Yes, do a normalized add
	mov	[esi-4], eax		; Store needs-normalize counter

; Do an unnormalized add
; ecx = src arg #1
; edx = src arg #2
; esi = dest arg

uadd:	mov	eax, addcount1		; Load loop counter
uadd1:	mov	edi, blkrows		; Cache lines in a block
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
	sub	edi, 1			; Test inner counter
	jnz	short uaddlp		; Loop if necessary
	lea	ecx, [ecx-blkrows*64]	; Next source
	lea	edx, [edx-blkrows*64]	; Next source
	lea	esi, [esi-blkrows*64]	; Next dest
	add	ecx, pass1blkdst	; Next source
	add	edx, pass1blkdst	; Next source
	add	esi, pass1blkdst	; Next dest
	sub	eax, 1			; Check loop counter
	jnz	short uadd1		; Loop if necessary
	pop	edi			; Restore registers
	pop	esi			; Restore registers
	ret

; Do a normalized add
; ecx = src arg #1
; edx = src arg #2
; esi = dest arg

nadd:	mov	DWORD PTR [esi-4], 0	; Clear needs-normalize flag
	push	ebp			; Save registers
	push	ebx			; Save registers
	mov	eax, addcount1		; Get number of blocks
	mov	loopcount1, eax		; Save outer loop counter
	mov	eax, norm_grp_mults	; Addr of the group multipliers
	mov	edi, norm_biglit_array	; Addr of the big/little flags array
	mov	ebp, carries		; Get pointer to carries
	cmp	_NUMLIT, 0		; Test for irrational FFTs
	JNE_X	nadd0   		; Yes, use two-to-phi multipliers

radd0:	mov	loopcount3, blkrows	; Save inner loop count
radd2:	xnorm_op_2d addpd, noexec	; Add and normalize 8 values
	sub	loopcount3, 1		; Decrement inner loop counter
	JNZ_X	radd2 			; Loop til done
	lea	ecx, [ecx-blkrows*64]	; Next source
	lea	edx, [edx-blkrows*64]	; Next source
	lea	esi, [esi-blkrows*64]	; Next dest
	add	ecx, pass1blkdst	; Next source
	add	edx, pass1blkdst	; Next source
	add	esi, pass1blkdst	; Next dest
	lea	ebp, [ebp+64]		; Next set of carries
	sub	loopcount1, 1		; Decrement outer loop counter
	JNZ_X	radd0 			; Loop til done
	JMP_X	nadddn			; Skip to finishing code

nadd0:	mov	ebx, normval1		; Load middle loop count
	mov	loopcount2, ebx		; Save middle loop count
	mov	ebx, norm_col_mults	; Addr of the column multipliers
nadd1:	push	esi
	mov	esi, count2		; Load inner loop count
	mov	loopcount3, esi		; Save inner loop count
	pop	esi
nadd2:	xnorm_op_2d addpd, exec		; Add and normalize 8 values
	sub	loopcount3, 1		; Decrement inner loop counter
	JNZ_X	nadd2 			; Loop til done
	add	edi, normval2		; Adjust ptr to little/big flags
	sub	loopcount2, 1		; Decrement middle loop counter
	JNZ_X	nadd1			; Loop til done
	lea	ecx, [ecx-blkrows*64]	; Next source
	lea	edx, [edx-blkrows*64]	; Next source
	lea	esi, [esi-blkrows*64]	; Next dest
	add	ecx, pass1blkdst	; Next source
	add	edx, pass1blkdst	; Next source
	add	esi, pass1blkdst	; Next dest
	lea	eax, [eax+128]		; Next set of group multipliers
	lea	ebp, [ebp+64]		; Next set of carries
	add	edi, normval3		; Adjust little/big flags ptr
	sub	loopcount1, 1		; Decrement outer loop counter
	JNZ_X	nadd0 			; Loop til done

nadddn:	CALL_X	gw_norm_cleanup		; Split the carries
	pop	ebx			; Restore registers
	pop	ebp			; Restore registers
	pop	edi			; Restore registers
	pop	esi			; Restore registers
	ret


; Rotate the accumulated carries and add them into the FFT data.
; This is very similar to gw_split_carries_3 except that after an
; add or subtract the carries will be very small and only need to
; be rotated, not split into high and low halves.

gw_norm_cleanup:
	mov	edx, count3		; Load 3 section counts
	mov	loopcount1, edx		; Save for later
	mov	edx, addcount1		; Load count of carry rows
	mov	esi, carries		; Addr of the carries
	mov	ebp, edx		; Compute addr of the high carries
	shl	ebp, 6
	add	ebp, esi
	mov	edx, norm_grp_mults	; Addr of the group multipliers
	xnorm012_2d_part1
	mov	ebp, _DESTARG		; Addr of FFT data
	movapd	xmm6, XMM_BIGVAL
clp0:	mov	eax, loopcount1		; Get list of counts
	mov	ebx, eax		; Form count for this section
	and	ebx, 03FFh
	JZ_X	cdn			; No rows to do.  We're all done!
	mov	loopcount2, ebx		; Save count of carry rows this section
	shr	eax, 10			; Move counts list along
	mov	loopcount1, eax
	shl	ebx, 6			; Compute addr of the last carries row
	add	ebx, esi
	xnorm012_2d_part2
clp1:	xnorm_op012_2d			; Split carries for one cache line
	lea	esi, [esi+64]		; Next carries pointer
	add	ebp, pass1blkdst	; Next FFT data pointer
	cmp	_NUMLIT, 0		; Don't bump group mult pointer
	je	short cskip		; for rational FFTs
	lea	edx, [edx+128]		; Next group multiplier
cskip:	sub	loopcount2, 1		; Test loop counter
	JNZ_X	clp1			; Next carry row in section
	JMP_X	clp0			; Next section
cdn:	ret




;;
;; Subtract two numbers without carry propogation.  Caller can use this for
;; consecutive add or subtract operations.  However, the last operation
;; before a multiply must use the routine that will normalize data.
;;

	PUBLIC	gwxsubf3
	PUBLIC	gwxsubq3
gwxsubf3:
gwxsubq3:
	mov	ecx, _SRCARG		; Address of first number
	mov	edx, _SRC2ARG		; Address of second number
	push	esi			; Save registers
	push	edi			; Save registers
	mov	eax, [ecx-4]		; Load needs-normalize counter
	mov	edi, [edx-4]		; Load needs-normalize counter
	mov	esi, _DESTARG		; Address of destination
	lea	eax, [eax+edi+1]	; Set new needs-normalize counter
	mov	edi, [ecx-28]		; Load has-been-FFTed flag
	mov	[esi-28], edi		; Save has-been-FFTed flag
	jmp	short usub		; Do the unnormalized sub
	mov	[esi-4], eax		; Store needs-normalize counter

;;
;; Subtract two numbers.  This will be done without carry propogation
;; (for extra speed) if it can be done safely.
;;

	PUBLIC	gwxsub3
gwxsub3:

; See if subtract can be done without normalization

	mov	ecx, _SRCARG		; Address of first number
	mov	edx, _SRC2ARG		; Address of second number
	push	esi			; Save registers
	push	edi			; Save registers
	mov	eax, [ecx-4]		; Load needs-normalize counter
	mov	edi, [edx-4]		; Load needs-normalize counter
	mov	esi, _DESTARG		; Address of destination
	lea	eax, [eax+edi+1]	; Set new needs-normalize counter
	mov	DWORD PTR [esi-28], 0	; Save has-been-FFTed flag
	cmp	eax, extra_bits		; Is normalization needed?
	jg	nsub			; Yes, do a normalized subtract
	mov	[esi-4], eax		; Store needs-normalize counter

; Do an unnormalized subtract
; ecx = src arg #1
; edx = src arg #2
; esi = dest arg

usub:	mov	eax, addcount1		; Load loop counter
usub1:	mov	edi, blkrows		; Cache lines in a block
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
	sub	edi, 1			; Test inner counter
	jnz	short usublp		; Loop if necessary
	lea	ecx, [ecx-blkrows*64]	; Next source
	lea	edx, [edx-blkrows*64]	; Next source
	lea	esi, [esi-blkrows*64]	; Next dest
	add	ecx, pass1blkdst	; Next source
	add	edx, pass1blkdst	; Next source
	add	esi, pass1blkdst	; Next dest
	sub	eax, 1			; Check loop counter
	jnz	short usub1		; Loop if necessary
	pop	edi			; Restore registers
	pop	esi			; Restore registers
	ret

; Do a normalized subtract
; ecx = src arg #1
; edx = src arg #2
; esi = dest arg

nsub:	mov	DWORD PTR [esi-4], 0	; Clear needs-normalize flag
	push	ebp			; Save registers
	push	ebx			; Save registers
	mov	eax, addcount1		; Get number of blocks
	mov	loopcount1, eax		; Save outer loop counter
	mov	eax, norm_grp_mults	; Addr of the group multipliers
	mov	edi, norm_biglit_array	; Addr of the big/little flags array
	mov	ebp, carries		; Get pointer to carries
	cmp	_NUMLIT, 0		; Test for irrational FFTs
	JNE_X	nsub0   		; Yes, use two-to-phi multipliers

rsub0:	mov	loopcount3, blkrows	; Save inner loop count
rsub2:	xnorm_op_2d subpd, noexec	; Add and normalize 8 values
	sub	loopcount3, 1		; Decrement inner loop counter
	JNZ_X	rsub2 			; Loop til done
	lea	ecx, [ecx-blkrows*64]	; Next source
	lea	edx, [edx-blkrows*64]	; Next source
	lea	esi, [esi-blkrows*64]	; Next dest
	add	ecx, pass1blkdst	; Next source
	add	edx, pass1blkdst	; Next source
	add	esi, pass1blkdst	; Next dest
	lea	ebp, [ebp+64]		; Next set of carries
	sub	loopcount1, 1		; Decrement outer loop counter
	JNZ_X	rsub0 			; Loop til done
	JMP_X	nsubdn			; Skip to finishing code

nsub0:	mov	ebx, normval1		; Load middle loop count
	mov	loopcount2, ebx		; Save middle loop count
	mov	ebx, norm_col_mults	; Addr of the column multipliers
nsub1:	push	esi
	mov	esi, count2		; Load inner loop count
	mov	loopcount3, esi		; Save inner loop count
	pop	esi
nsub2:	xnorm_op_2d subpd, exec		; Add and normalize 8 values
	sub	loopcount3, 1		; Decrement inner loop counter
	JNZ_X	nsub2 			; Loop til done
	add	edi, normval2		; Adjust ptr to little/big flags
	sub	loopcount2, 1		; Decrement middle loop counter
	JNZ_X	nsub1			; Loop til done
	lea	ecx, [ecx-blkrows*64]	; Next source
	lea	edx, [edx-blkrows*64]	; Next source
	lea	esi, [esi-blkrows*64]	; Next dest
	add	ecx, pass1blkdst	; Next source
	add	edx, pass1blkdst	; Next source
	add	esi, pass1blkdst	; Next dest
	lea	eax, [eax+128]		; Next set of group multipliers
	lea	ebp, [ebp+64]		; Next set of carries
	add	edi, normval3		; Adjust little/big flags ptr
	sub	loopcount1, 1		; Decrement outer loop counter
	JNZ_X	nsub0 			; Loop til done

nsubdn:	CALL_X	gw_norm_cleanup		; Split the carries
	pop	ebx			; Restore registers
	pop	ebp			; Restore registers
	pop	edi			; Restore registers
	pop	esi			; Restore registers
	ret


;;
;; Add and subtract two numbers without carry propogation.
;;

	PUBLIC	gwxaddsubf3
	PUBLIC	gwxaddsubq3
gwxaddsubf3:
gwxaddsubq3:
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
	mov	edi, [ecx-28]		; Load has-been-FFTed flag
	mov	[esi-4], eax		; Store needs-normalize counter
	mov	[ebp-4], eax		; Store needs-normalize counter
	mov	[esi-28], edi		; Save has-been-FFTed flag
	mov	[ebp-28], edi		; Save has-been-FFTed flag
	jmp	short uaddsub		; Do the unnormalized add

;;
;; Add and subtract two numbers.  This will be done without carry propogation
;; (for extra speed) if it can be done safely.
;;

	PUBLIC	gwxaddsub3
gwxaddsub3:

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
	mov	DWORD PTR [esi-28], 0	; Save has-been-FFTed flag
	mov	DWORD PTR [ebp-28], 0	; Save has-been-FFTed flag
	cmp	eax, extra_bits		; Is normalization needed?
	JG_X	naddsub			; Yes, do a normalized add
	mov	[esi-4], eax		; Store needs-normalize counter
	mov	[ebp-4], eax		; Store needs-normalize counter

; Do an unnormalized add
; ecx = src arg #1
; edx = src arg #2
; esi = dest arg #1
; ebp = dest arg #2

uaddsub:mov	eax, addcount1		; Load loop counter
uaddsub1:mov	edi, blkrows		; Cache lines in a block
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
	sub	edi, 1			; Test inner counter
	jnz	uaddsublp		; Loop if necessary
	lea	ecx, [ecx-blkrows*64]	; Next source
	lea	edx, [edx-blkrows*64]	; Next source
	lea	esi, [esi-blkrows*64]	; Next dest
	lea	ebp, [ebp-blkrows*64]	; Next dest
	add	ecx, pass1blkdst	; Next source
	add	edx, pass1blkdst	; Next source
	add	esi, pass1blkdst	; Next dest
	add	ebp, pass1blkdst	; Next dest
	sub	eax, 1			; Check loop counter
	JNZ_X	uaddsub1		; Loop if necessary
	pop	ebp			; Restore registers
	pop	edi			; Restore registers
	pop	esi			; Restore registers
	ret

; Do a normalized add and subtract
; ecx = src arg #1
; edx = src arg #2
; esi = dest arg #1
; ebp = dest arg #2

naddsub:mov	DWORD PTR [esi-4], 0	; Clear needs-normalize counter
	mov	DWORD PTR [ebp-4], 0	; Clear needs-normalize counter
	push	ebx			; Save registers
	mov	eax, addcount1		; Get number of blocks
	mov	loopcount1, eax		; Save outer loop counter
	mov	eax, norm_grp_mults	; Addr of the group multipliers
	mov	normgrpptr, eax
	mov	edi, norm_biglit_array	; Addr of the big/little flags array
	mov	eax, carries		; Get pointer to carries
	cmp	_NUMLIT, 0		; Test for irrational FFTs
	JNE_X	nas0	   		; Yes, use two-to-phi multipliers

ras0:	mov	loopcount3, blkrows	; Save inner loop count
	movapd	xmm6, XMM_BIGVAL	; clear carries
	movapd	[eax], xmm6
	movapd	[eax+16], xmm6
	movapd	[eax+32], xmm6
	movapd	[eax+48], xmm6
	movapd	[eax+64], xmm6
	movapd	[eax+80], xmm6
	movapd	[eax+96], xmm6
	movapd	[eax+112], xmm6
ras2:	xnorm_addsub_2d noexec		; Add and normalize 8 values
	sub	loopcount3, 1		; Decrement inner loop counter
	JNZ_X	ras2 			; Loop til done
	lea	ecx, [ecx-blkrows*64]	; Next source
	lea	edx, [edx-blkrows*64]	; Next source
	lea	esi, [esi-blkrows*64]	; Next dest
	lea	ebp, [ebp-blkrows*64]	; Next dest
	add	ecx, pass1blkdst	; Next source
	add	edx, pass1blkdst	; Next source
	add	esi, pass1blkdst	; Next dest
	add	ebp, pass1blkdst	; Next dest
	lea	eax, [eax+128]		; Next set of carries
	sub	loopcount1, 1		; Decrement outer loop counter
	JNZ_X	ras0 			; Loop til done
	JMP_X	nasdn			; Skip to finishing code

nas0:	mov	ebx, normval1		; Load middle loop count
	mov	loopcount2, ebx		; Save middle loop count
	mov	ebx, norm_col_mults	; Addr of the column multipliers
	movapd	xmm6, XMM_BIGVAL	; clear carries
	movapd	[eax], xmm6
	movapd	[eax+16], xmm6
	movapd	[eax+32], xmm6
	movapd	[eax+48], xmm6
	movapd	[eax+64], xmm6
	movapd	[eax+80], xmm6
	movapd	[eax+96], xmm6
	movapd	[eax+112], xmm6
nas1:	push	esi
	mov	esi, count2		; Load inner loop count
	mov	loopcount3, esi		; Save inner loop count
	pop	esi
nas2:	xnorm_addsub_2d exec		; Add and normalize 8 values
	sub	loopcount3, 1		; Decrement inner loop counter
	JNZ_X	nas2 			; Loop til done
	add	edi, normval2		; Adjust ptr to little/big flags
	sub	loopcount2, 1		; Decrement middle loop counter
	JNZ_X	nas1			; Loop til done
	lea	ecx, [ecx-blkrows*64]	; Next source
	lea	edx, [edx-blkrows*64]	; Next source
	lea	esi, [esi-blkrows*64]	; Next dest
	lea	ebp, [ebp-blkrows*64]	; Next dest
	add	ecx, pass1blkdst	; Next source
	add	edx, pass1blkdst	; Next source
	add	esi, pass1blkdst	; Next dest
	add	ebp, pass1blkdst	; Next dest
	add	normgrpptr, 128		; Next set of group multipliers
	lea	eax, [eax+128]		; Next set of carries
	add	edi, normval3		; Adjust little/big flags ptr
	sub	loopcount1, 1		; Decrement outer loop counter
	JNZ_X	nas0 			; Loop til done

; Rotate the accumulated carries and add them into the FFT data.
; This is very similar to gw_split_carries_3 except that after an
; add or subtract the carries will be very small and only need to
; be rotated, not split into high and low halves.

nasdn:	mov	edx, count3		; Load 3 section counts
	mov	loopcount1, edx		; Save for later
	mov	edx, addcount1		; Load count of carry rows
	mov	esi, carries		; Addr of the carries
	mov	ebp, edx		; Compute addr of the high carries
	shl	ebp, 7
	add	ebp, esi
	mov	edx, norm_grp_mults	; Addr of the group multipliers
	xnorm_as012_2d_part1
	mov	ebp, _DESTARG		; Addr of FFT data
	mov	ecx, _DEST2ARG
	movapd	xmm6, XMM_BIGVAL
asclp0:	mov	eax, loopcount1		; Get list of counts
	mov	ebx, eax		; Form count for this section
	and	ebx, 03FFh
	JZ_X	ascdn			; No rows to do.  We're all done!
	mov	loopcount2, ebx		; Save count of carry rows this section
	shr	eax, 10			; Move counts list along
	mov	loopcount1, eax
	shl	ebx, 7			; Compute addr of the last carries row
	add	ebx, esi
	xnorm_as012_2d_part2
asclp1:	xnorm_as012_2d			; Split carries for one cache line
	lea	esi, [esi+128]		; Next carries pointer
	add	ebp, pass1blkdst	; Next FFT data pointer
	add	ecx, pass1blkdst	; Next FFT data pointer
	cmp	_NUMLIT, 0		; Don't bump group mult pointer
	je	short acskip		; for rational FFTs
	lea	edx, [edx+128]		; Next group multiplier
acskip:	sub	loopcount2, 1		; Test loop counter
	JNZ_X	asclp1			; Next carry row in section
	JMP_X	asclp0			; Next section
ascdn:	pop	ebx			; Restore registers
	pop	ebp			; Restore registers
	pop	edi			; Restore registers
	pop	esi			; Restore registers
	ret

gwx3procs ENDP

;;
;; Copy one number and zero some low order words.
;;

	PUBLIC	gwxcopyzero3
gwxcopyzero3 PROC NEAR
	push	esi			; Save registers
	push	edi			; Save registers
	mov	esi, _SRCARG		; Address of first number
	mov	edi, _DESTARG		; Address of destination
	mov	eax, [esi-4]		; Load needs-normalize counter
	mov	[edi-4], eax		; Store needs-normalize counter
	sub	ecx, ecx		; Offset to compare to COPYZERO
	mov	DWORD PTR [edi-28], ecx	; Save has-been-FFTed flag
	mov	eax, addcount1		; Get number of blocks
cz1:	mov	edx, blkrows		; Number of cache lines in one block
cz2:	xcopyzero 8
	lea	esi, [esi+64]		; Next source
	lea	edi, [edi+64]		; Next dest
	lea	ecx, [ecx+64]		; Next compare offset
	sub	edx, 1			; Test loop counter
	JNZ_X	cz2			; Loop if necessary
	lea	esi, [esi-blkrows*64]
	lea	edi, [edi-blkrows*64]
	lea	ecx, [ecx-blkrows*64]
	add	esi, pass1blkdst
	add	edi, pass1blkdst
	add	ecx, pass1blkdst
	sub	eax, 1			; Test loop counter
	JNZ_X	cz1			; Loop if necessary
	pop	edi			; Restore registers
	pop	esi			; Restore registers
	ret
gwxcopyzero3 ENDP

;;
;; Do a mod k*2^n+/-1.  This is grossly inefficient because only
;; one-fourth of a cache line is processed before it leaves the L2 cache.
;;

	PUBLIC	gwxprothmod3
gwxprothmod3 PROC NEAR
	push	esi			; Save registers
	push	edi
	push	ebp
	push	ebx
	mov	esi, _SRCARG		; Address of first number
	prothmod_upper_prep_0d		; Prepare for prothmod stage 1
	mov	ecx, count3		; Load 3 section counts
	shl	ecx, 2
pmu0:	sub	eax, eax		; Form count for this section
	shld	eax, ecx, 10
	shl	ecx, 10			; Move counts list along
	and	eax, eax		; Make sure we got a non-zero count
	jz	short pmu0
	mov	ebx, eax		; Compute size of FFT in bytes
	imul	ebx, pass1blkdst
	mov	loopcount1, eax
pmu1:	mov	ebp, blkrows		; Cache lines in a block
pmu2:	prothmod_upper_0d 64, 2*64, 4*64, pmudn ; Divide 8 upper values by k
	lea	esi, [esi-8*64]		; Next source word
	lea	edi, [edi+8*8]		; Next scratch words
	sub	ebp, 8			; Test loop counter
	JNZ_X	pmu2			; Iterate
	lea	esi, [esi+blkrows*64-8]	; Next source pointer
	add	eax, 80000000h		; Test loop counter
	JNC_X	pmu1			; Iterate if necessary
	lea	esi, [esi+2*8]		; Next source pointer
	sub	esi, pass1blkdst	; Next source pointer
	sub	eax, 1			; Test loop counter
	JNZ_X	pmu1			; Iterate if necessary
	add	esi, ebx		; Restore source address
	sub	esi, 16			; Next source address
	mov	eax, loopcount1		; Restore count of blocks
	add	cl, 256/4		; Test loop counter
	JNC_X	pmu1			; Iterate if necessary
	lea	esi, [esi+4*16]		; Restore source pointer
	sub	esi, ebx		; Next source pointer
	JMP_X	pmu0			; Do next set of blocks
pmudn:	prothmod_lower_prep_0d		; Prepare for prothmod stage 2
	mov	ecx, count3		; Load 3 section counts
pml0:	mov	eax, ecx		; Form count for this section
	and	eax, 03FFh
	shr	ecx, 10			; Move counts list along
	mov	ebx, eax		; Compute size of FFT in bytes
	imul	ebx, pass1blkdst
	mov	loopcount1, eax
pml1:	mov	ebp, blkrows		; Cache lines in a block
pml2:	prothmod_lower_0d 64, 2*64, 4*64; Add 8 shifted quotients to lower data
	lea	esi, [esi+8*64]		; Next source word
	lea	edi, [edi-8*8]		; Next scratch words
	dec	edx
	jz	short pmldn
	sub	ebp, 8			; Test loop counter
	JNZ_X	pml2			; Iterate
	lea	esi, [esi-blkrows*64+8]	; Next source pointer
	add	eax, 80000000h		; Test loop counter
	JNC_X	pml1			; Iterate if necessary
	lea	esi, [esi-2*8]		; Next source pointer
	add	esi, pass1blkdst	; Next source pointer
	sub	eax, 1			; Test loop counter
	JNZ_X	pml1			; Iterate if necessary
	sub	esi, ebx		; Restore source address
	add	esi, 16			; Next source address
	mov	eax, loopcount1		; Restore count of blocks
	add	ecx, 40000000h		; Test loop counter
	JNC_X	pml1			; Iterate if necessary
	lea	esi, [esi-4*16]		; Restore source pointer
	add	esi, ebx		; Next source pointer
	JMP_X	pml0			; Next set of blocks
pmldn:	prothmod_final_0d 64, 2*64, 4*64; Finish off the mod
	pop	ebx			; Restore registers
	pop	ebp
	pop	edi
	pop	esi
	ret
gwxprothmod3 ENDP

;;
;; Routines to do the normalization after a multiply
;;

_xnorm3 PROC NEAR

; Macro to loop through all the FFT values and apply the proper normalization
; routine.

inorm	MACRO	lab, ttp, zero, echk, const
	LOCAL	noadd, ilp0, ilp1, ilexit
	PUBLIC	lab
lab:
zero	cmp	zero_fft, 0		;; Is the first call to inorm?
zero	JE_X	zinit			;; Yes, do special initialization
	movapd	xmm7, XMM_SUMOUT	;; Load SUMOUT
	movapd	xmm6, XMM_MAXERR	;; Load maximum error
	cmp	edx, _ADDIN_ROW		;; Is this the time to do our addin?
	jne	short noadd		;; Jump if addin does not occur now
	mov	edi, _ADDIN_OFFSET	;; Get address to add value into
	_movsd	xmm0, [esi][edi]	;; Get the value
	addsd	xmm0, _ADDIN_VALUE	;; Add in the requested value
	_movsd	[esi][edi], xmm0	;; Save the new value
	subsd	xmm7, _ADDIN_VALUE	;; Do not include addin in sumout
noadd:	push	edx
	mov	edx, norm_grp_mults	;; Addr of the group multipliers
	mov	ebp, carries		;; Addr of the carries
	mov	edi, norm_ptr1		;; Load big/little flags array ptr
	mov	eax, addcount1		;; Load loop counter
	mov	loopcount2, eax		;; Save loop counter
	mov	loopcount3, 0		;; Clear outermost loop counter
	mov	eax, 48			;; Rational FFTs are all big words
	mov	ecx, 48
ilp0:	mov	ebx, count2		;; Load inner loop counter
	mov	loopcount1, ebx		;; Save loop counter
	mov	ebx, norm_ptr2		;; Load column multipliers ptr
ilp1:	xnorm_2d ttp, zero, echk, const	;; Normalize 8 values
	lea	esi, [esi+64]		;; Next cache line
ttp	lea	ebx, [ebx+32]		;; Next column multipliers
ttp	lea	edi, [edi+4]		;; Next big/little flags
	sub	loopcount1, 1		;; Test loop counter
	JNZ_X	ilp1			;; Loop til done
	add	esi, normblkdst		;; Skip gap in blkdst or clmblkdst
	lea	ebp, [ebp+64]		;; Next set of carries
ttp	lea	edx, [edx+128]		;; Next set of 8 group multipliers
	sub	loopcount2, 1		;; Test loop counter
	JZ_X	ilexit			;; Jump when loop complete
	add	loopcount3, 80000000h/4 ;; 8 iterations
	jnc	ilp0
	add	esi, normblkdst8	;; Add 128 every 8 clmblkdsts
	jmp	ilp0			;; Iterate
ilexit:	movapd	XMM_SUMOUT, xmm7	;; Save SUMOUT
	movapd	XMM_MAXERR, xmm6	;; Save maximum error
ttp	mov	norm_ptr1, edi		;; Save big/little flags array ptr
ttp	mov	norm_ptr2, ebx		;; Save column multipliers ptr
	pop	edx
	sub	ebx, ebx
	ret
	ENDM

; The 16 different normalization routines.  One for each combination of
; rational/irrational, zeroing/no zeroing, error check/no error check, and
; mul by const/no mul by const.

	inorm	xr3, noexec, noexec, noexec, noexec
	inorm	xr3e, noexec, noexec, exec, noexec
	inorm	xr3c, noexec, noexec, noexec, exec
	inorm	xr3ec, noexec, noexec, exec, exec
	inorm	xr3z, noexec, exec, noexec, noexec
	inorm	xr3ze, noexec, exec, exec, noexec
	inorm	xi3, exec, noexec, noexec, noexec
	inorm	xi3e, exec, noexec, exec, noexec
	inorm	xi3c, exec, noexec, noexec, exec
	inorm	xi3ec, exec, noexec, exec, exec
	inorm	xi3z, exec, exec, noexec, noexec
	inorm	xi3ze, exec, exec, exec, noexec

; Special init code for FFTs that zero the high words

zinit:	mov	zero_fft, 1		; Set flag saying we've initialized
	mov	eax, _DESTARG		; Add dest addr to 8 fftzero ptrs
	add	_FFTZERO, eax
	add	_FFTZERO+4, eax
	add	_FFTZERO+8, eax
	add	_FFTZERO+12, eax
	add	_FFTZERO+16, eax
	add	_FFTZERO+20, eax
	add	_FFTZERO+24, eax
	add	_FFTZERO+28, eax
	mov	eax, _NORMRTN		; Go back and do the normalize
	jmp	eax

_xnorm3 ENDP

_TEXT32	ENDS
END
