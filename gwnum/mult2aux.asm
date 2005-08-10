; Copyright 1995-2005 Just For Fun Software, Inc., all rights reserved
; Author:  George Woltman
; Email: woltman@alum.mit.edu
;

	TITLE   setup

	.686
	.XMM
	.MODEL	FLAT

INCLUDE	unravel.mac
INCLUDE extrn.mac
INCLUDE mult.mac
INCLUDE memory.mac
INCLUDE normal.mac

_TEXT SEGMENT

	flat_distances

;;
;; Add two numbers without carry propogation.  Caller can use this for
;; consecutive add or subtract operations.  However, the last operation
;; before a multiply must use the routine that will normalize data.
;;

IFNDEF PPRO
	PUBLIC	gwaddq2
gwaddq2	PROC NEAR
	push	esi			; Save registers
	push	edi			; Save registers
	push	ebx			; Save registers
	mov	ecx, _SRCARG		; Address of first number
	mov	edx, _SRC2ARG		; Address of second number
	mov	esi, _DESTARG		; Address of destination
	mov	ebx, addcount1		; Load blk count
uadd0:	mov	eax, normval4		; Load count of 4KB pages in a block
uaddlp:	fld	QWORD PTR [ecx]		; Load first number
	fadd	QWORD PTR [edx]		; Add in second number
	fld	QWORD PTR [ecx+8]	; Load first number
	fadd	QWORD PTR [edx+8]	; Add in second number
	fld	QWORD PTR [ecx+16]	; Load first number
	fadd	QWORD PTR [edx+16]	; Add in second number
	fxch	st(1)
	fld	QWORD PTR [ecx+24]	; Load first number
	fadd	QWORD PTR [edx+24]	; Add in second number
	fxch	st(3)
	fstp	QWORD PTR [esi]		; Save result
	fstp	QWORD PTR [esi+8]	; Save result
	fstp	QWORD PTR [esi+16]	; Save result
	fstp	QWORD PTR [esi+24]	; Save result
	lea	ecx, [ecx+32]		; Bump source pointer
	lea	edx, [edx+32]		; Bump source pointer
	lea	esi, [esi+32]		; Bump dest pointer
	add	eax, 80000000h/64	; 128 cache lines in a 4KB page
	jnc	short uaddlp		; Loop if necessary
	lea	ecx, [ecx+64]		; Skip 64 bytes every 4KB
	lea	edx, [edx+64]		; Skip 64 bytes every 4KB
	lea	esi, [esi+64]		; Skip 64 bytes every 4KB
	dec	eax			; Check middle loop counter
	jnz	short uaddlp		; Loop if necessary
	lea	ecx, [ecx+64]		; Skip 64 bytes every blk
	lea	edx, [edx+64]		; Skip 64 bytes every blk
	lea	esi, [esi+64]		; Skip 64 bytes every blk
	dec	ebx			; Check outer/restore inner count
	jnz	short uadd0		; Loop if necessary
	pop	ebx
	pop	edi
	pop	esi
	ret
gwaddq2	ENDP


;;
;; Add two numbers with carry propogation
;;

	PUBLIC	gwadd2
gwadd2	PROC NEAR
	push	esi			; Save registers
	push	edi			; Save registers
	push	ebp			; Save registers
	push	ebx			; Save registers
	mov	ecx, _SRCARG		; Address of first number
	mov	edx, _SRC2ARG		; Address of second number
	mov	esi, _DESTARG		; Address of destination
	mov	ebp, norm_grp_mults	; Address of group multipliers
	mov	edi, norm_biglit_array	; Addr of the big/little flags array
	fninit
	fld	BIGVAL			; Start process with no carry
	fld	BIGVAL
	mov	eax, addcount1		; Load block count
	mov	loopcount1, eax
ablk:	mov	eax, normval4		; Load outer count (4KB pages in a blk)
	mov	loopcount2, eax
	mov	ebx, norm_col_mults	; Addr of the column multipliers
iadd0:	mov	eax, normval1		; Load middle count (clms in 4KB page)
	mov	loopcount3, eax		; Save middle count
iadd1:	mov	eax, cache_line_multiplier ; Load inner loop count (clm)
	mov	loopcount4, eax		; Save inner loop count
	sub	eax, eax		; Clear big/lit flag
iadd2:	norm_op_2d fadd			; Add and normalize 4 values
	sub	loopcount4, 1		; Decrement inner loop counter
	jnz	iadd2 			; Loop til done
	add	edi, normval2		; Adjust ptr to little/big flags
	sub	loopcount3, 1		; Decrement middle loop counter
	jnz	iadd1			; Loop til done
	lea	ecx, [ecx+64]		; Skip 64 bytes every 4KB
	lea	edx, [edx+64]		; Skip 64 bytes every 4KB
	lea	esi, [esi+64]		; Skip 64 bytes every 4KB
	sub	loopcount2, 1		; Decrement outer loop counter
	jnz	iadd0			; Loop til done
	lea	ecx, [ecx+64]		; Skip 64 bytes every blk
	lea	edx, [edx+64]		; Skip 64 bytes every blk
	lea	esi, [esi+64]		; Skip 64 bytes every blk
	lea	ebp, [ebp+2*16]		; Next set of 2 group multipliers
	add	edi, normval3		; Adjust little/big flags ptr
	sub	loopcount1, 1		; Decrement outer loop counter
	jnz	ablk 			; Loop til done

	;; All blocks done

	mov	esi, _DESTARG		; Addr of FFT data
	mov	ebp, norm_grp_mults	; Addr of the group multipliers
	norm_op_2d_cleanup		; Add 2 carries to start of fft

	pop	ebx			; Restore registers
	pop	ebp
	pop	edi
	pop	esi
	ret
gwadd2	ENDP


;;
;; Subtract two numbers without carry propogation.  Caller can use this for
;; consecutive add or subtract operations.  However, the last operation
;; before a multiply must use the routine that will normalize data.
;;

	PUBLIC	gwsubq2
gwsubq2	PROC NEAR
	push	esi			; Save registers
	push	edi			; Save registers
	push	ebx			; Save registers
	mov	ecx, _SRCARG		; Address of first number
	mov	edx, _SRC2ARG		; Address of second number
	mov	esi, _DESTARG		; Address of destination
	mov	ebx, addcount1		; Load blk count
usub0:	mov	eax, normval4		; Load count of 4KB pages in a block
usublp:	fld	QWORD PTR [edx]		; Load second number
	fsub	QWORD PTR [ecx]		; Subtract first number
	fld	QWORD PTR [edx+8]	; Load second number
	fsub	QWORD PTR [ecx+8]	; Subtract first number
	fld	QWORD PTR [edx+16]	; Load second number
	fsub	QWORD PTR [ecx+16]	; Subtract first number
	fxch	st(1)
	fld	QWORD PTR [edx+24]	; Load second number
	fsub	QWORD PTR [ecx+24]	; Subtract first number
	fxch	st(3)
	fstp	QWORD PTR [esi]		; Save result
	fstp	QWORD PTR [esi+8]	; Save result
	fstp	QWORD PTR [esi+16]	; Save result
	fstp	QWORD PTR [esi+24]	; Save result
	lea	ecx, [ecx+32]		; Bump source pointer
	lea	edx, [edx+32]		; Bump source pointer
	lea	esi, [esi+32]		; Bump dest pointer
	add	eax, 80000000h/64	; 128 cache lines in a 4KB page
	jnc	short usublp		; Loop if necessary
	lea	ecx, [ecx+64]		; Skip 64 bytes every 4KB
	lea	edx, [edx+64]		; Skip 64 bytes every 4KB
	lea	esi, [esi+64]		; Skip 64 bytes every 4KB
	dec	eax			; Check middle loop counter
	jnz	short usublp		; Loop if necessary
	lea	ecx, [ecx+64]		; Skip 64 bytes every blk
	lea	edx, [edx+64]		; Skip 64 bytes every blk
	lea	esi, [esi+64]		; Skip 64 bytes every blk
	dec	ebx			; Check outer/restore inner count
	jnz	short usub0		; Loop if necessary
	pop	ebx
	pop	edi
	pop	esi
	ret
gwsubq2	ENDP

;;
;; Subtract two numbers with carry propogation
;;

	PUBLIC	gwsub2
gwsub2	PROC NEAR
	push	esi			; Save registers
	push	edi			; Save registers
	push	ebp			; Save registers
	push	ebx			; Save registers
	mov	ecx, _SRCARG		; Address of first number
	mov	edx, _SRC2ARG		; Address of second number
	mov	esi, _DESTARG		; Address of destination
	mov	ebp, norm_grp_mults	; Address of group multipliers
	mov	edi, norm_biglit_array	; Addr of the big/little flags array
	fninit
	fld	BIGVAL			; Start process with no carry
	fld	BIGVAL
	mov	eax, addcount1		; Load block count
	mov	loopcount1, eax

sblk:	mov	eax, normval4		; Load outer count (4KB pages in a blk)
	mov	loopcount2, eax
	mov	ebx, norm_col_mults	; Addr of the column multipliers
isub0:	mov	eax, normval1		; Load middle count (clms in 4KB page)
	mov	loopcount3, eax		; Save middle count
isub1:	mov	eax, cache_line_multiplier ; Load inner loop count (clm)
	mov	loopcount4, eax		; Save inner loop count
	sub	eax, eax		; Clear big/lit flag
isub2:	norm_op_2d fsub			; Add and normalize 4 values
	sub	loopcount4, 1		; Decrement inner loop counter
	jnz	isub2 			; Loop til done
	add	edi, normval2		; Adjust ptr to little/big flags
	sub	loopcount3, 1		; Decrement middle loop counter
	jnz	isub1			; Loop til done
	lea	ecx, [ecx+64]		; Skip 64 bytes every 4KB
	lea	edx, [edx+64]		; Skip 64 bytes every 4KB
	lea	esi, [esi+64]		; Skip 64 bytes every 4KB
	sub	loopcount2, 1		; Decrement outer loop counter
	jnz	isub0			; Loop til done
	lea	ecx, [ecx+64]		; Skip 64 bytes every blk
	lea	edx, [edx+64]		; Skip 64 bytes every blk
	lea	esi, [esi+64]		; Skip 64 bytes every blk
	lea	ebp, [ebp+2*16]		; Next set of 2 group multipliers
	add	edi, normval3		; Adjust little/big flags ptr
	sub	loopcount1, 1		; Decrement outer loop counter
	jnz	sblk 			; Loop til done

	;; All blocks done

	mov	esi, _DESTARG		; Addr of FFT data
	mov	ebp, norm_grp_mults	; Addr of the group multipliers
	norm_op_2d_cleanup		; Add 2 carries to start of fft

	pop	ebx			; Restore registers
	pop	ebp
	pop	edi
	pop	esi
	ret
gwsub2	ENDP

;;
;; Add and subtract two numbers without carry propogation.
;;

	PUBLIC	gwaddsubq2
gwaddsubq2 PROC NEAR
	push	esi			; Save registers
	push	edi			; Save registers
	push	ebp			; Save registers
	push	ebx			; Save registers
	mov	ecx, _SRCARG		; Address of first number
	mov	edx, _SRC2ARG		; Address of second number
	mov	esi, _DESTARG		; Address of destination
	mov	ebp, _DEST2ARG  	; Address of destination #2
	mov	ebx, addcount1		; Load blk count
uaddsub0:
	mov	eax, normval4		; Load count of 4KB pages in a block
uaddsublp:
	fld	QWORD PTR [ecx]		; Load first number
	fld	st(0)			; Dup first number
	fadd	QWORD PTR [edx]		; Add in second number
	fxch	st(1)			; S0,A0
	fsub	QWORD PTR [edx]		; Subtract out second number
	fld	QWORD PTR [ecx+8]	; Load first number
	fld	st(0)			; Dup first number
	fadd	QWORD PTR [edx+8]	; Add in second number
	fxch	st(1)			; S1,A1,S0,A0
	fsub	QWORD PTR [edx+8]	; Subtract out second number
	fld	QWORD PTR [ecx+16]	; Load first number
	fld	st(0)			; Dup first number
	fadd	QWORD PTR [edx+16]	; Add in second number
	fxch	st(1)			; S2,A2,S1,A1,S0,A0
	fsub	QWORD PTR [edx+16]	; Subtract out second number
	fld	QWORD PTR [ecx+24]	; Load first number
	fld	st(0)			; Dup first number
	fadd	QWORD PTR [edx+24]	; Add in second number
	fxch	st(7)			; A0,S3,S2,A2,S1,A1,S0,A3
	fstp	QWORD PTR [esi]		; Save result
	fsub	QWORD PTR [edx+24]	; Subtract out second number
	fxch	st(5)			; S0,S2,A2,S1,A1,S3,A3
	fstp	QWORD PTR [ebp]		; Save result
	fstp	QWORD PTR [ebp+16]	; Save result
	fstp	QWORD PTR [esi+16]	; Save result
	fstp	QWORD PTR [ebp+8]	; Save result
	fstp	QWORD PTR [esi+8]	; Save result
	fstp	QWORD PTR [ebp+24]	; Save result
	fstp	QWORD PTR [esi+24]	; Save result
	lea	ecx, [ecx+32]		; Bump source pointer
	lea	edx, [edx+32]		; Bump source pointer
	lea	esi, [esi+32]		; Bump dest pointer
	lea	ebp, [ebp+32]		; Bump dest pointer
	add	eax, 80000000h/64	; 128 cache lines in a 4KB page
	jnc	short uaddsublp		; Loop if necessary
	lea	ecx, [ecx+64]		; Skip 64 bytes every 4KB
	lea	edx, [edx+64]		; Skip 64 bytes every 4KB
	lea	esi, [esi+64]		; Skip 64 bytes every 4KB
	lea	ebp, [ebp+64]		; Skip 64 bytes every 4KB
	dec	eax			; Check middle loop counter
	jnz	short uaddsublp		; Loop if necessary
	lea	ecx, [ecx+64]		; Skip 64 bytes every blk
	lea	edx, [edx+64]		; Skip 64 bytes every blk
	lea	esi, [esi+64]		; Skip 64 bytes every blk
	lea	ebp, [ebp+64]		; Skip 64 bytes every blk
	dec	ebx			; Check outer/restore inner count
	jnz	short uaddsub0		; Loop if necessary
	pop	ebx
	pop	ebp
	pop	edi
	pop	esi
	ret
gwaddsubq2 ENDP

;;
;; Add and subtract two numbers with carry propogation
;;

	PUBLIC	gwaddsub2
gwaddsub2 PROC NEAR
	push	esi			; Save registers
	push	edi			; Save registers
	push	ebp			; Save registers
	push	ebx			; Save registers
	mov	ecx, _SRCARG		; Address of first number
	mov	edx, _SRC2ARG		; Address of second number
	mov	esi, _DESTARG		; Address of destination
	mov	ebp, _DEST2ARG  	; Address of destination #2
	mov	eax, norm_grp_mults	; Address of group multipliers
	mov	edi, norm_biglit_array	; Addr of the big/little flags array
	fninit
	fld	BIGVAL			; Start process with no carry
	fld	BIGVAL
	fld	BIGVAL
	fld	BIGVAL
	mov	ebx, addcount1		; Load block count
	mov	loopcount1, ebx
asblk:	mov	ebx, normval4		; Load outer count (4KB pages in a blk)
	mov	loopcount2, ebx
	mov	ebx, norm_col_mults	; Addr of the column multipliers
ias0:	push	eax
	mov	eax, normval1		; Load middle count (clms in 4KB page)
	mov	loopcount3, eax		; Save middle count
	pop	eax
ias1:	push	eax
	mov	eax, cache_line_multiplier ; Load inner loop count (clm)
	mov	loopcount4, eax		; Save inner loop count
	pop	eax
ias2:	norm_addsub_2d			; Add and normalize 4 values
	sub	loopcount4, 1		; Decrement inner loop counter
	jnz	ias2 			; Loop til done
	add	edi, normval2		; Adjust ptr to little/big flags
	sub	loopcount3, 1		; Decrement middle loop counter
	jnz	ias1			; Loop til done
	lea	ecx, [ecx+64]		; Skip 64 bytes every 4KB
	lea	edx, [edx+64]		; Skip 64 bytes every 4KB
	lea	esi, [esi+64]		; Skip 64 bytes every 4KB
	lea	ebp, [ebp+64]		; Skip 64 bytes every 4KB
	sub	loopcount2, 1		; Decrement outer loop counter
	jnz	ias0			; Loop til done
	lea	ecx, [ecx+64]		; Skip 64 bytes every blk
	lea	edx, [edx+64]		; Skip 64 bytes every blk
	lea	esi, [esi+64]		; Skip 64 bytes every blk
	lea	ebp, [ebp+64]		; Skip 64 bytes every blk
	lea	eax, [eax+2*16]		; Next set of 2 group multipliers
	add	edi, normval3		; Adjust little/big flags ptr
	sub	loopcount1, 1		; Decrement outer loop counter
	jnz	asblk 			; Loop til done

	;; All blocks done

	mov	esi, _DESTARG		; Addr of add FFT data
	mov	ebp, _DEST2ARG  	; Addr of sub FFT data
	mov	ebx, norm_grp_mults	; Addr of the group multipliers
	norm_addsub_2d_cleanup		; Add 2 carries to start of fft

	pop	ebx			; Restore registers
	pop	ebp
	pop	edi
	pop	esi
	ret
gwaddsub2 ENDP

;;
;; Copy one number and zero some low order words.
;;

	PUBLIC	gwcopyzero2
gwcopyzero2 PROC NEAR
	push	esi			; Save registers
	push	edi			; Save registers
	push	ebx
	mov	esi, _SRCARG		; Address of first number
	mov	edi, _DESTARG		; Address of destination
	mov	ecx, _NUMARG		; Count of words to zero
	mov	ebx, addcount1		; Load blk count
zlp0:	mov	eax, normval4		; Load count of 4KB pages in a block
zlp:	dec	ecx			; Decrement count of words to zero
	js	short c2		; If negative we are now copying
	fldz
	dec	ecx			; Decrement count of words to zero
	js	short c1		; If negative we are now copying
	fldz
	jmp	short c0
c2:	fld	QWORD PTR [esi]		; Load number
c1:	fld	QWORD PTR [esi+16]	; Load number
c0:	fstp	QWORD PTR [edi+16]	; Save result
	fstp	QWORD PTR [edi]		; Save result
	fld	QWORD PTR [esi+8]	; Copy high word
	fstp	QWORD PTR [edi+8]
	fld	QWORD PTR [esi+24]	; Copy high word
	fstp	QWORD PTR [edi+24]
	lea	esi, [esi+32]		; Bump source pointer
	lea	edi, [edi+32]		; Bump dest pointer
	add	eax, 80000000h/64	; 128 cache lines in a 4KB page
	jnc	short zlp		; Loop if necessary
	lea	esi, [esi+64]		; Skip 64 bytes every 4KB
	lea	edi, [edi+64]		; Skip 64 bytes every 4KB
	dec	eax			; Check middle loop counter
	jnz	short zlp		; Loop if necessary
	lea	esi, [esi+64]		; Skip 64 bytes every blk
	lea	edi, [edi+64]		; Skip 64 bytes every blk
	dec	ebx			; Check outer/restore inner count
	jnz	short zlp0		; Loop if necessary
	pop	ebx
	pop	edi
	pop	esi
	ret
gwcopyzero2 ENDP

ENDIF

;; Routines to do the normalization after a multiply

;; When doing zero-padded FFTs, the multiplied 7 words around the halfway point
;; must be subtracted from the bottom of the FFT.  This must be done before
;; normalization multiplies the FFT data by k.  This macro does that.

sub_7_words MACRO
	LOCAL	nozpad, zlp
	cmp	zpad_addr, OFFSET XMM_ZPAD6;; Have we subtracted all 7 words?
	jg	short nozpad		;; Yes, skip this code
	mov	edi, zpad_addr		;; Addr of next zpad element to process
	mov	eax, cache_line_multiplier ;; Load loop counter
	add	eax, eax		;; Two values per cache line
	push	esi
zlp:	fld	QWORD PTR [esi]		;; Load FFT word
	fld	QWORD PTR [edi]		;; Load ZPAD data
	fmul	XMM_NORM012_FF		;; Scale by FFTLEN/2
	fsub	st(1), st
	faddp	st(2), st		;; Adjust sumout
	fstp	QWORD PTR [esi]
	lea	esi, [esi+dist1]	;; Bump pointers
	lea	edi, [edi+8]
	dec	eax			;; Iterate 2*clm (up to 8) times
	jnz	short zlp		;; Loop if necessary
	pop	esi			;; Restore source ptr
	mov	zpad_addr, edi
nozpad:
	ENDM

IFNDEF PPRO
_comp6norm PROC NEAR

; Macro to loop through all the FFT values and apply the proper normalization
; routine.  Used whenever we are using an irrational-base FFT.

inorm	MACRO	lab, ttp, zero, echk, const
	LOCAL	noadd, ilp0, ilp1, ilpdn
	PUBLICP	lab
	LABELP	lab
zero	mov	zero_fft, 1		;; Set flag saying zero upper half
echk	fld	_MAXERR			;; Load maximum error
	fld	XMM_SUMOUT		;; Load SUMOUT
no zero	cmp	edx, _ADDIN_ROW		;; Is this the time to do our addin?
no zero	jne	short noadd		;; Jump if addin does not occur now
no zero	mov	edi, _ADDIN_OFFSET	;; Get address to add value into
no zero	fld	QWORD PTR [esi][edi]	;; Get the value
no zero	fadd	_ADDIN_VALUE		;; Add in the requested value
no zero	fstp	QWORD PTR [esi][edi]	;; Save the new value
no zero	fsub	_ADDIN_VALUE		;; Do not include addin in sumout
noadd:	push	edx
	push	esi
	mov	edx, norm_grp_mults	;; Addr of the group multipliers
	mov	ebp, carries		;; Addr of the carries
	mov	edi, norm_ptr1		;; Load big/little flags array ptr
	mov	eax, addcount1		;; Load loop counter
	mov	loopcount2, eax		;; Save loop counter
	mov	loopcount3, 0		;; Clear outermost loop counter
	sub	eax, eax		;; Clear big/lit flags
ilp0:	mov	ebx, cache_line_multiplier ;; Load inner loop counter
	mov	loopcount1, ebx		;; Save loop counter
	mov	ebx, norm_ptr2		;; Load column multipliers ptr
	fld	QWORD PTR [ebp+0*8]	;; Load carry1
	fadd	BIGVAL			;; c1 = c1 + rounding constant
	fld	QWORD PTR [ebp+1*8]	;; Load carry2
	fadd	BIGVAL			;; c2 = c2 + rounding constant
	push	esi
ilp1:	norm_2d ttp, zero, echk, const	;; Normalize 4 vals (2 grps, 2 cols)
	lea	esi, [esi+4*8]		;; Next FFT source
ttp	lea	ebx, [ebx+2*16]		;; Next column multipliers
ttp	lea	edi, [edi+2]		;; Next big/little flags
	sub	loopcount1, 1		;; Test clm loop counter
	jnz	ilp1			;; Loop til done
	pop	esi
	fsub	BIGVAL			;; c2 = c2 - rounding constant
	fstp	QWORD PTR [ebp+1*8]	;; Save carry2
	fsub	BIGVAL			;; c1 = c1 - rounding constant
	fstp	QWORD PTR [ebp+0*8]	;; Save carry1
	add	esi, normblkdst		;; Next source pointer
	lea	ebp, [ebp+2*8]		;; Next set of carries
ttp	lea	edx, [edx+2*16]		;; Next set of 2 group multipliers
	sub	loopcount2, 1		;; Test outer loop counter
	jz	short ilpdn		;; Iterate
	add	loopcount3, 80000000h/16;; 32 iterations
	jnc	ilp0
	add	esi, normblkdst8	;; Add 64 pad bytes every 32 clmblkdsts
	jmp	ilp0			;; Iterate
ilpdn:	fstp	XMM_SUMOUT		;; Save SUMOUT
echk	fstp	_MAXERR			;; Save maximum error
ttp	mov	norm_ptr1, edi		;; Save big/little flags array ptr
ttp	mov	norm_ptr2, ebx		;; Save column multipliers ptr
	pop	esi
	pop	edx
	sub	ebx, ebx
	cmp	edx, 65536+256		;; Check for last iteration
	je	top_carry_adjust_2	;; Top carry may require adjusting
	ret
	ENDM

zpnorm	MACRO	lab, ttp, echk, const
	LOCAL	ilp0, ilp1, ilpdn
	PUBLICP	lab
	LABELP	lab
const	mov	const_fft, 1		;; Set flag saying mul-by-const
echk	fld	_MAXERR			;; Load maximum error
	fld	XMM_SUMOUT		;; Load SUMOUT
	sub_7_words
	push	edx
	push	esi
	mov	edx, norm_grp_mults	;; Addr of the group multipliers
	mov	ebp, carries		;; Addr of the carries
	mov	edi, norm_ptr1		;; Load big/little flags array ptr
	mov	eax, addcount1		;; Load loop counter
	mov	loopcount2, eax		;; Save loop counter
	mov	loopcount3, 0		;; Clear outermost loop counter
	sub	eax, eax		;; Clear big/lit flags
ilp0:	mov	ebx, cache_line_multiplier ;; Load inner loop counter
	mov	loopcount1, ebx		;; Save loop counter
	mov	ebx, norm_ptr2		;; Load column multipliers ptr
	fld	QWORD PTR [ebp+0*8]	;; Load carry1
	fadd	BIGVAL			;; c1 = c1 + rounding constant
	fld	QWORD PTR [ebp+1*8]	;; Load carry2
	push	esi
ilp1:	norm_2d_zpad ttp, echk, const	;; Normalize 4 vals (2 grps, 2 cols)
	lea	esi, [esi+4*8]		;; Next FFT source
ttp	lea	ebx, [ebx+2*16]		;; Next column multipliers
ttp	lea	edi, [edi+2]		;; Next big/little flags
	sub	loopcount1, 1		;; Test clm loop counter
	jnz	ilp1			;; Loop til done
	pop	esi
	fstp	QWORD PTR [ebp+1*8]	;; Save carry2
	fsub	BIGVAL			;; c1 = c1 - rounding constant
	fstp	QWORD PTR [ebp+0*8]	;; Save carry1
	add	esi, normblkdst		;; Next source pointer
	lea	ebp, [ebp+2*8]		;; Next set of carries
ttp	lea	edx, [edx+2*16]		;; Next set of 2 group multipliers
	sub	loopcount2, 1		;; Test outer loop counter
	jz	short ilpdn		;; Iterate
	add	loopcount3, 80000000h/16;; 32 iterations
	jnc	ilp0
	add	esi, normblkdst8	;; Add 64 pad bytes every 32 clmblkdsts
	jmp	ilp0			;; Iterate
ilpdn:	fstp	XMM_SUMOUT		;; Save SUMOUT
echk	fstp	_MAXERR			;; Save maximum error
ttp	mov	norm_ptr1, edi		;; Save big/little flags array ptr
ttp	mov	norm_ptr2, ebx		;; Save column multipliers ptr
	pop	esi
	pop	edx
	sub	ebx, ebx
	ret
	ENDM

; The 16 different normalization routines

	inorm	r2, noexec, noexec, noexec, noexec
	inorm	r2e, noexec, noexec, exec, noexec
	inorm	r2c, noexec, noexec, noexec, exec
	inorm	r2ec, noexec, noexec, exec, exec
	inorm	r2z, noexec, exec, noexec, noexec
	inorm	r2ze, noexec, exec, exec, noexec
	inorm	i2, exec, noexec, noexec, noexec
	inorm	i2e, exec, noexec, exec, noexec
	inorm	i2c, exec, noexec, noexec, exec
	inorm	i2ec, exec, noexec, exec, exec
	inorm	i2z, exec, exec, noexec, noexec
	inorm	i2ze, exec, exec, exec, noexec

	zpnorm	r2zp, noexec, noexec, noexec
	zpnorm	r2zpe, noexec, exec, noexec
	zpnorm	r2zpc, noexec, noexec, exec
	zpnorm	r2zpec, noexec, exec, exec
	zpnorm	i2zp, exec, noexec, noexec
	zpnorm	i2zpe, exec, exec, noexec
	zpnorm	i2zpc, exec, noexec, exec
	zpnorm	i2zpec, exec, exec, exec

; Special code to handle adjusting the carry out of the topmost FFT word

top_carry_adjust_2:
	norm_top_carry_2d
	ret

_comp6norm ENDP
ENDIF

_TEXT	ENDS
END
