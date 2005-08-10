; Copyright 2001-2005 Just For Fun Software, Inc., all rights reserved
; Author:  George Woltman
; Email: woltman@alum.mit.edu
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
INCLUDE memory.mac
INCLUDE xnormal.mac

_TEXT SEGMENT

gwx2procs PROC

;;
;; Add two numbers without carry propogation.  Caller can use this for
;; consecutive add or subtract operations.  However, the last operation
;; before a multiply must use the routine that will normalize data.
;;

	PUBLIC	gwxaddf2
gwxaddf2:
	PUBLIC	gwxaddq2
gwxaddq2:
	push	rsi			; Save registers
	push	rdi			; Save registers
	push	rbx
	mov	rcx, _SRCARG		; Address of first number
	mov	rdx, _SRC2ARG		; Address of second number
	mov	rsi, _DESTARG		; Address of destination
	mov	ebx, addcount1		; Load blk count
uadd0:	mov	eax, normval4		; Load count of 8KB chunks in a block
uaddlp:	movapd	xmm0, [rdx]		; Load second number
	addpd	xmm0, [rcx]		; Add in first number
	movapd	xmm1, [rdx+16]		; Load second number
	addpd	xmm1, [rcx+16]		; Add in first number
	movapd	xmm2, [rdx+32]		; Load second number
	addpd	xmm2, [rcx+32]		; Add in first number
	movapd	xmm3, [rdx+48]		; Load second number
	addpd	xmm3, [rcx+48]		; Add in first number
	movapd	[rsi], xmm0		; Save result
	movapd	[rsi+16], xmm1		; Save result
	movapd	[rsi+32], xmm2		; Save result
	movapd	[rsi+48], xmm3		; Save result
	lea	rcx, [rcx+64]		; Next source
	lea	rdx, [rdx+64]		; Next source
	lea	rsi, [rsi+64]		; Next dest
	add	eax, 80000000h/64	; 128 cache lines in a 8KB chunk
	jnc	short uaddlp		; Loop if necessary
	lea	rcx, [rcx+128]		; Skip 128 bytes every 8KB
	lea	rdx, [rdx+128]		; Skip 128 bytes every 8KB
	lea	rsi, [rsi+128]		; Skip 128 bytes every 8KB
	dec	rax			; Check middle loop counter
	jnz	short uaddlp		; Loop if necessary
	add	rcx, pass2gapsize	; Next source
	add	rdx, pass2gapsize	; Next source
	add	rsi, pass2gapsize	; Next dest
	dec	rbx			; Check loop counter
	jnz	uadd0			; Loop if necessary
	pop	rbx
	pop	rdi			; Restore registers
	pop	rsi			; Restore registers
	ret

;;
;; Add two numbers with carry propogation
;;

	PUBLIC	gwxadd2
gwxadd2:
	push	rsi			; Save registers
	push	rdi			; Save registers
	push	rbp			; Save registers
	push	rbx			; Save registers
	IFDEF X86_64
	sub	rsp, 16
	_movsd	[rsp+8], xmm6
	_movsd	[rsp], xmm7
	ENDIF
	mov	rcx, _SRCARG		; Address of first number
	mov	rdx, _SRC2ARG		; Address of second number
	mov	rsi, _DESTARG		; Address of destination
	mov	rbp, norm_grp_mults	; Addr of the group multipliers
	mov	rdi, norm_biglit_array	; Addr of the big/little flags array
	movapd	xmm7, XMM_BIGVAL	; Init 4 carries
	movapd	XMM_TMP1, xmm7
	movapd	XMM_TMP2, xmm7
	movapd	XMM_TMP3, xmm7
	movapd	XMM_TMP4, xmm7
	mov	eax, count3		; Load 3 section counts

	;; Do a section

asec:	mov	loopcount1, eax		; Save section counts
	and	eax, 07FFh		; Form block count for this section
	mov	loopcount2, eax		; Save count of blocks in this section
	mov	norm_ptr1, rsi		; Save section start address
	mov	norm_ptr2, rbp		; Save section group ttp address

	;; Do a block

ablk:	mov	eax, normval4		; Load count of 8KB chunks in a blk
	mov	loopcount3, eax
	mov	rbx, norm_col_mults	; Addr of the column multipliers
add0:	cmp	_RATIONAL_FFT, 0	; Test for irrational FFTs
	je	iadd0   		; Yes, use two-to-phi multipliers

	mov	loopcount4, 128		; Cache lines in 8KB chunk
	sub	rax, rax		; Clear big/lit flag
radd1:	xnorm_op_2d addpd, noexec	; Add and normalize 8 values
	dec	loopcount4		; Decrement inner loop counter
	jnz	radd1 			; Loop til done
	jmp	achunkdn		; Jump to chunk done code

iadd0:	mov	eax, normval1		; Load count of clms in 8KB chunk
	mov	loopcount4, eax
iadd1:	mov	eax, cache_line_multiplier ; Load inner loop count
	mov	loopcount5, eax		; Save inner loop count
	sub	rax, rax		; Clear big/lit flag
iadd2:	xnorm_op_2d addpd, exec		; Add and normalize 8 values
	dec	loopcount5		; Decrement inner loop counter
	jnz	iadd2 			; Loop til done
	add	rdi, normval2		; Adjust ptr to little/big flags
	dec	loopcount4		; Decrement middle loop counter
	jnz	iadd1			; Loop til done

	;; Chunk done

achunkdn:lea	rcx, [rcx+128]		; Skip 128 bytes every 8KB
	lea	rdx, [rdx+128]		; Skip 128 bytes every 8KB
	lea	rsi, [rsi+128]		; Skip 128 bytes every 8KB
	dec	loopcount3		; Test loop counter
	jnz	add0
	add	rcx, pass2gapsize	; Next source
	add	rdx, pass2gapsize	; Next source
	add	rsi, pass2gapsize	; Next dest

	;; Block done

ablkdn:	sub	rsi, pass1blkdst	; Restore start of block ptr
					; Add 4 carries to start of block
	xnorm_op_2d_blk rsi, rbp, XMM_TMP1, XMM_TMP2, XMM_TMP3, XMM_TMP4
	add	rsi, pass1blkdst
	lea	rbp, [rbp+128]		; Next set of group multipliers
	add	rdi, normval3		; Adjust little/big flags ptr
	dec	loopcount2		; Decrement outer loop counter
	jnz	ablk 			; Loop til done

	;; Section done

	mov	rax, norm_ptr1		; Reload section start ptr
	mov	rbx, norm_ptr2		; Reload section group ttp ptr
					; Add 2 carries to start of section
	xnorm_op_2d_sec XMM_TMP1, XMM_TMP2, XMM_TMP3, XMM_TMP4
	mov	eax, loopcount1		; Get list of counts
	shr	eax, 11			; Move counts list along
	jnz	asec			; Do next section if necessary

	;; All sections done

	mov	rsi, _DESTARG		; Addr of FFT data
	movapd	xmm7, XMM_TMP3		; Load wraparound carry
	xnorm_top_carry_cmn rsi, xmm7, 2
	mov	rbp, norm_grp_mults	; Addr of the group multipliers
	movapd	xmm6, XMM_TMP1		; Load non-wraparound carry
	xnorm_op_2d_fft			; Add 2 carries to start of fft

	IFDEF X86_64
	_movsd	xmm7, [rsp]
	_movsd	xmm6, [rsp+8]
	add	rsp, 16
	ENDIF
	pop	rbx			; Restore registers
	pop	rbp			; Restore registers
	pop	rdi			; Restore registers
	pop	rsi			; Restore registers
	ret



;;
;; Subtract two numbers without carry propogation.  Caller can use this for
;; consecutive add or subtract operations.  However, the last operation
;; before a multiply must use the routine that will normalize data.
;;

	PUBLIC	gwxsubf2
gwxsubf2:
	PUBLIC	gwxsubq2
gwxsubq2:
	push	rsi			; Save registers
	push	rdi			; Save registers
	push	rbx
	mov	rcx, _SRCARG		; Address of first number
	mov	rdx, _SRC2ARG		; Address of second number
	mov	rsi, _DESTARG		; Address of destination
	mov	ebx, addcount1		; Load blk count
usub0:	mov	eax, normval4		; Load count of 8KB chunks in a block
usublp:	movapd	xmm0, [rdx]		; Load second number
	subpd	xmm0, [rcx]		; Subtract first number
	movapd	xmm1, [rdx+16]		; Load second number
	subpd	xmm1, [rcx+16]		; Subtract first number
	movapd	xmm2, [rdx+32]		; Load second number
	subpd	xmm2, [rcx+32]		; Subtract first number
	movapd	xmm3, [rdx+48]		; Load second number
	subpd	xmm3, [rcx+48]		; Subtract first number
	movapd	[rsi], xmm0		; Save result
	movapd	[rsi+16], xmm1		; Save result
	movapd	[rsi+32], xmm2		; Save result
	movapd	[rsi+48], xmm3		; Save result
	lea	rcx, [rcx+64]		; Next source
	lea	rdx, [rdx+64]		; Next source
	lea	rsi, [rsi+64]		; Next dest
	add	eax, 80000000h/64	; 128 cache lines in a 8KB chunk
	jnc	short usublp		; Loop if necessary
	lea	rcx, [rcx+128]		; Skip 128 bytes every 8KB
	lea	rdx, [rdx+128]		; Skip 128 bytes every 8KB
	lea	rsi, [rsi+128]		; Skip 128 bytes every 8KB
	dec	rax			; Check middle loop counter
	jnz	short usublp		; Loop if necessary
	add	rcx, pass2gapsize	; Next source
	add	rdx, pass2gapsize	; Next source
	add	rsi, pass2gapsize	; Next dest
	dec	rbx			; Check loop counter
	jnz	usub0			; Loop if necessary
	pop	rbx
	pop	rdi			; Restore registers
	pop	rsi			; Restore registers
	ret

;;
;; Subtract two numbers with carry propogation
;;

	PUBLIC	gwxsub2
gwxsub2:
	push	rsi			; Save registers
	push	rdi			; Save registers
	push	rbp			; Save registers
	push	rbx			; Save registers
	IFDEF X86_64
	sub	rsp, 16
	_movsd	[rsp+8], xmm6
	_movsd	[rsp], xmm7
	ENDIF
	mov	rcx, _SRCARG		; Address of first number
	mov	rdx, _SRC2ARG		; Address of second number
	mov	rsi, _DESTARG		; Address of destination
	mov	rbp, norm_grp_mults	; Addr of the group multipliers
	mov	rdi, norm_biglit_array	; Addr of the big/little flags array
	movapd	xmm7, XMM_BIGVAL	; Init 4 carries
	movapd	XMM_TMP1, xmm7
	movapd	XMM_TMP2, xmm7
	movapd	XMM_TMP3, xmm7
	movapd	XMM_TMP4, xmm7
	mov	eax, count3		; Load 3 section counts

	;; Do a section

ssec:	mov	loopcount1, eax		; Save section counts
	and	eax, 07FFh		; Form block count for this section
	mov	loopcount2, eax		; Save count of blocks in this section
	mov	norm_ptr1, rsi		; Save section start address
	mov	norm_ptr2, rbp		; Save section group ttp address

	;; Do a block

sblk:	mov	eax, normval4		; Load count of 8KB chunks in a blk
	mov	loopcount3, eax
	mov	rbx, norm_col_mults	; Addr of the column multipliers
sub0:	cmp	_RATIONAL_FFT, 0	; Test for irrational FFTs
	je	isub0   		; Yes, use two-to-phi multipliers

	mov	loopcount4, 128		; Cache lines in 8KB chunk
	sub	rax, rax		; Clear big/lit flag
rsub1:	xnorm_op_2d subpd, noexec	; Subtract and normalize 8 values
	dec	loopcount4		; Decrement inner loop counter
	jnz	rsub1 			; Loop til done
	jmp	schunkdn		; Jump to chunk done code

isub0:	mov	eax, normval1		; Load count of clms in 8KB chunk
	mov	loopcount4, eax
isub1:	mov	eax, cache_line_multiplier ; Load inner loop count
	mov	loopcount5, eax		; Save inner loop count
	sub	rax, rax		; Clear big/lit flag
isub2:	xnorm_op_2d subpd, exec		; Subtract and normalize 8 values
	dec	loopcount5		; Decrement inner loop counter
	jnz	isub2 			; Loop til done
	add	rdi, normval2		; Adjust ptr to little/big flags
	dec	loopcount4		; Decrement middle loop counter
	jnz	isub1			; Loop til done

	;; Chunk done

schunkdn:lea	rcx, [rcx+128]		; Skip 128 bytes every 8KB
	lea	rdx, [rdx+128]		; Skip 128 bytes every 8KB
	lea	rsi, [rsi+128]		; Skip 128 bytes every 8KB
	dec	loopcount3		; Test loop counter
	jnz	sub0
	add	rcx, pass2gapsize	; Next source
	add	rdx, pass2gapsize	; Next source
	add	rsi, pass2gapsize	; Next dest

	;; Block done

sblkdn:	sub	rsi, pass1blkdst	; Restore start of block ptr
					; Add 4 carries to start of block
	xnorm_op_2d_blk rsi, rbp, XMM_TMP1, XMM_TMP2, XMM_TMP3, XMM_TMP4
	add	rsi, pass1blkdst
	lea	rbp, [rbp+128]		; Next set of group multipliers
	add	rdi, normval3		; Adjust little/big flags ptr
	dec	loopcount2		; Decrement outer loop counter
	jnz	sblk 			; Loop til done

	;; Section done

	mov	rax, norm_ptr1		; Reload section start ptr
	mov	rbx, norm_ptr2		; Reload section group ttp ptr
					; Add 2 carries to start of section
	xnorm_op_2d_sec XMM_TMP1, XMM_TMP2, XMM_TMP3, XMM_TMP4
	mov	eax, loopcount1		; Get list of counts
	shr	eax, 11			; Move counts list along
	jnz	ssec			; Do next section if necessary

	;; All sections done

	mov	rsi, _DESTARG		; Addr of FFT data
	movapd	xmm7, XMM_TMP3		; Load wraparound carry
	xnorm_top_carry_cmn rsi, xmm7, 2
	mov	rbp, norm_grp_mults	; Addr of the group multipliers
	movapd	xmm6, XMM_TMP1		; Load non-wraparound carry
	xnorm_op_2d_fft			; Add 2 carries to start of fft

	IFDEF X86_64
	_movsd	xmm7, [rsp]
	_movsd	xmm6, [rsp+8]
	add	rsp, 16
	ENDIF
	pop	rbx			; Restore registers
	pop	rbp			; Restore registers
	pop	rdi			; Restore registers
	pop	rsi			; Restore registers
	ret


;;
;; Add and subtract two numbers without carry propogation.
;;

	PUBLIC	gwxaddsubf2
gwxaddsubf2:
	PUBLIC	gwxaddsubq2
gwxaddsubq2:
	push	rsi			; Save registers
	push	rdi			; Save registers
	push	rbp			; Save registers
	push	rbx
	IFDEF X86_64
	sub	rsp, 16
	_movsd	[rsp+8], xmm6
	_movsd	[rsp], xmm7
	ENDIF
	mov	rcx, _SRCARG		; Address of first number
	mov	rdx, _SRC2ARG		; Address of second number
	mov	rsi, _DESTARG		; Address of destination #1
	mov	rbp, _DEST2ARG  	; Address of destination #2
	mov	ebx, addcount1		; Load blk count
uaddsub0:mov	eax, normval4		; Load count of 8KB chunks in a block
uaddsublp:
	movapd	xmm0, [rcx]		; Load first number
	movapd	xmm1, xmm0		; Dup first number
	addpd	xmm0, [rdx]		; Add in second number
	subpd	xmm1, [rdx]		; Subtract out second number
	movapd	xmm2, [rcx+16]		; Load first number
	movapd	xmm3, xmm2		; Dup first number
	addpd	xmm2, [rdx+16]		; Add in second number
	subpd	xmm3, [rdx+16]		; Subtract out second number
	movapd	xmm4, [rcx+32]		; Load first number
	movapd	xmm5, xmm4		; Dup first number
	addpd	xmm4, [rdx+32]		; Add in second number
	subpd	xmm5, [rdx+32]		; Subtract out second number
	movapd	xmm6, [rcx+48]		; Load first number
	movapd	xmm7, xmm6		; Dup first number
	addpd	xmm6, [rdx+48]		; Add in second number
	subpd	xmm7, [rdx+48]		; Subtract out second number
	movapd	[rsi], xmm0		; Save result
	movapd	[rbp], xmm1		; Save result
	movapd	[rsi+16], xmm2		; Save result
	movapd	[rbp+16], xmm3		; Save result
	movapd	[rsi+32], xmm4		; Save result
	movapd	[rbp+32], xmm5		; Save result
	movapd	[rsi+48], xmm6		; Save result
	movapd	[rbp+48], xmm7		; Save result
	lea	rcx, [rcx+64]		; Next source
	lea	rdx, [rdx+64]		; Next source
	lea	rsi, [rsi+64]		; Next dest
	lea	rbp, [rbp+64]		; Next dest
	add	eax, 80000000h/64	; 128 cache lines in a 8KB chunk
	jnc	uaddsublp		; Loop if necessary
	lea	rcx, [rcx+128]		; Skip 128 bytes every 8KB
	lea	rdx, [rdx+128]		; Skip 128 bytes every 8KB
	lea	rsi, [rsi+128]		; Skip 128 bytes every 8KB
	lea	rbp, [rbp+128]		; Skip 128 bytes every 8KB
	dec	rax			; Check middle loop counter
	jnz	uaddsublp		; Loop if necessary
	add	rcx, pass2gapsize	; Next source
	add	rdx, pass2gapsize	; Next source
	add	rsi, pass2gapsize	; Next dest
	add	rbp, pass2gapsize	; Next dest
	dec	rbx			; Check loop counter
	jnz	uaddsub0		; Loop if necessary
	IFDEF X86_64
	_movsd	xmm7, [rsp]
	_movsd	xmm6, [rsp+8]
	add	rsp, 16
	ENDIF
	pop	rbx
	pop	rbp			; Restore registers
	pop	rdi			; Restore registers
	pop	rsi			; Restore registers
	ret

;;
;; Add and subtract two numbers with carry propogation
;;

	PUBLIC	gwxaddsub2
gwxaddsub2:
	push	rsi			; Save registers
	push	rdi			; Save registers
	push	rbp			; Save registers
	push	rbx			; Save registers
	IFDEF X86_64
	sub	rsp, 16
	_movsd	[rsp+8], xmm6
	_movsd	[rsp], xmm7
	ENDIF
	mov	rcx, _SRCARG		; Address of first number
	mov	rdx, _SRC2ARG		; Address of second number
	mov	rsi, _DESTARG		; Address of destination
	mov	rbp, _DEST2ARG  	; Address of destination #2
	mov	rbx, norm_grp_mults	; Addr of the group multipliers
	mov	rdi, norm_biglit_array	; Addr of the big/little flags array
	movapd	xmm7, XMM_BIGVAL	; Init 8 carries
	movapd	XMM_TMP1, xmm7
	movapd	XMM_TMP2, xmm7
	movapd	XMM_TMP3, xmm7
	movapd	XMM_TMP4, xmm7
	movapd	XMM_TMP5, xmm7
	movapd	XMM_TMP6, xmm7
	movapd	XMM_TMP7, xmm7
	movapd	XMM_TMP8, xmm7
	mov	eax, count3		; Load 3 section counts

	;; Do a section

assec:	mov	loopcount1, eax		; Save section counts
	and	eax, 07FFh		; Form block count for this section
	mov	loopcount2, eax		; Save count of blocks in this section
	mov	norm_ptr1, rsi		; Save section start address
	mov	norm_ptr2, rbp		; Save section start address
	mov	norm_ptr3, rbx		; Save section group ttp address

	;; Do a block

asblk:	mov	eax, normval4		; Load count of 8KB chunks in a blk
	mov	loopcount3, eax
	mov	rax, norm_col_mults	; Addr of the column multipliers
as0:	cmp	_RATIONAL_FFT, 0	; Test for irrational FFTs
	je	ias0	   		; Yes, use two-to-phi multipliers

	mov	loopcount4, 128		; Cache lines in an 8KB chunk
	sub	rax, rax		; Clear big/lit flag
ras1:	xnorm_addsub_2d noexec		; Add & subtract and normalize 8 values
	dec	loopcount4		; Decrement inner loop counter
	jnz	ras1 			; Loop til done
	jmp	aschunkdn		; Jump to chunk done code

ias0:	push	rax
	mov	eax, normval1		; Load count of clms in 8KB chunk
	mov	loopcount4, eax
	pop	rax
ias1:	push	rax
	mov	eax, cache_line_multiplier ; Load inner loop count
	mov	loopcount5, eax		; Save inner loop count
	pop	rax
ias2:	xnorm_addsub_2d exec		; Add & subtract and normalize 8 values
	dec	loopcount5		; Decrement inner loop counter
	jnz	ias2 			; Loop til done
	add	rdi, normval2		; Adjust ptr to little/big flags
	dec	loopcount4		; Decrement middle loop counter
	jnz	ias1			; Loop til done

	;; Chunk done

aschunkdn:lea	rcx, [rcx+128]		; Skip 128 bytes every 8KB
	lea	rdx, [rdx+128]		; Skip 128 bytes every 8KB
	lea	rsi, [rsi+128]		; Skip 128 bytes every 8KB
	lea	rbp, [rbp+128]		; Skip 128 bytes every 8KB
	dec	loopcount3		; Test loop counter
	jnz	as0
	add	rcx, pass2gapsize	; Next source
	add	rdx, pass2gapsize	; Next source
	add	rsi, pass2gapsize	; Next dest
	add	rbp, pass2gapsize	; Next dest

	;; Block done

asblkdn:sub	rsi, pass1blkdst	; Restore start of block ptr
					; Add 4 carries to start of block
	xnorm_op_2d_blk rsi, rbx, XMM_TMP1, XMM_TMP2, XMM_TMP3, XMM_TMP4
	sub	rbp, pass1blkdst	; Restore start of block ptr
					; Add 4 carries to start of block
	xnorm_op_2d_blk rbp, rbx, XMM_TMP5, XMM_TMP6, XMM_TMP7, XMM_TMP8
	add	rsi, pass1blkdst
	add	rbp, pass1blkdst
	lea	rbx, [rbx+128]		; Next set of group multipliers
	add	rdi, normval3		; Adjust little/big flags ptr
	dec	loopcount2		; Decrement outer loop counter
	jnz	asblk 			; Loop til done

	;; Section done

	push	rbx
	mov	rax, norm_ptr1		; Reload section start ptr for dest #1
	mov	rbx, norm_ptr3		; Reload section group ttp ptr
					; Add 2 carries to start of section
	xnorm_op_2d_sec XMM_TMP1, XMM_TMP2, XMM_TMP3, XMM_TMP4
	mov	rax, norm_ptr2		; Reload section start ptr for dest #2
					; Add 2 carries to start of section
	xnorm_op_2d_sec XMM_TMP5, XMM_TMP6, XMM_TMP7, XMM_TMP8
	pop	rbx
	mov	eax, loopcount1		; Get list of counts
	shr	eax, 11			; Move counts list along
	jnz	assec			; Do next section if necessary

	;; All sections done

	mov	rsi, _DESTARG		; Addr of FFT data
	movapd	xmm7, XMM_TMP3		; Load wraparound carry
	xnorm_top_carry_cmn rsi, xmm7, 2
	mov	rbp, norm_grp_mults	; Addr of the group multipliers
	movapd	xmm6, XMM_TMP1		; Load non-wraparound carry
	xnorm_op_2d_fft			; Add 2 carries to start of fft

	mov	rsi, _DEST2ARG		; Addr of FFT data
	movapd	xmm7, XMM_TMP7		; Load wraparound carry
	xnorm_top_carry_cmn rsi, xmm7, 2
	movapd	xmm6, XMM_TMP5		; Load non-wraparound carry
	xnorm_op_2d_fft			; Add 2 carries to start of fft

	IFDEF X86_64
	_movsd	xmm7, [rsp]
	_movsd	xmm6, [rsp+8]
	add	rsp, 16
	ENDIF
	pop	rbx			; Restore registers
	pop	rbp			; Restore registers
	pop	rdi			; Restore registers
	pop	rsi			; Restore registers
	ret

gwx2procs ENDP

;;
;; Copy one number and zero some low order words.
;;

	PUBLIC	gwxcopyzero2
gwxcopyzero2 PROC
	push	rsi			; Save registers
	push	rdi			; Save registers
	push	rbx
	mov	rsi, _SRCARG		; Address of first number
	mov	rdi, _DESTARG		; Address of destination
	sub	ecx, ecx		; Offset to compare to COPYZERO
	mov	ebx, addcount1		; Get number of blocks
cz1:	mov	eax, normval4		; Load count of 8KB chunks in a block
cz2:	xcopyzero
	lea	rsi, [rsi+64]		; Next source
	lea	rdi, [rdi+64]		; Next dest
	lea	rcx, [rcx+64]		; Next compare offset
	add	eax, 80000000h/64	; 128 cache lines in a 8KB chunk
	jnc	short cz2		; Loop if necessary
	lea	rsi, [rsi+128]		; Skip 128 bytes every 8KB
	lea	rdi, [rdi+128]		; Skip 128 bytes every 8KB
	lea	rcx, [rcx+128]		; Skip 128 bytes every 8KB
	dec	rax			; Test loop counter
	jnz	cz2			; Loop if necessary
	add	rsi, pass2gapsize
	add	rdi, pass2gapsize
	add	rcx, pass2gapsize
	dec	rbx			; Test loop counter
	jnz	cz1			; Loop if necessary
	pop	rbx
	pop	rdi			; Restore registers
	pop	rsi			; Restore registers
	ret
gwxcopyzero2 ENDP

_TEXT	ENDS
END
