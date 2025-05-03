; Copyright 2024 Mersenne Research, Inc.  All rights reserved
; Author:  George Woltman
; Email: woltman@alum.mit.edu
;

	TITLE   setup

INCLUDE	unravel.mac
INCLUDE extrn.mac
INCLUDE zarch.mac
INCLUDE zbasics.mac
INCLUDE znormal.mac
INCLUDE znormal_zpad.mac

_TEXT SEGMENT

;; Only assemble the add/sub/etc routines when compiled with the CORE architecture

IF (@INSTR(,%zarch,<SKX>) NE 0)

;; General register layout for routines below (mostly compatible with zmult3ax.asm layout):

;; rsi		source1 pointer
;; rcx		source2 pointer
;; rbx		dest pointer
;; rbp		dest2 pointer (addsub only)
;; rdi		big/lit flags ptr
;; r9		saved dest ptr at section start
;; r8		saved big/lit ptr at section start
;; rax		outer loop counter
;; r10		inner loop counter
;; r12		compressed biglit table pointer
;; r13		distance to source #2
;; r14		distance to source #4
;; r15		saved dest ptr #2 at section start (for addsub)
;; zmm0-7	carries
;; zmm8-15	carries (addsub only)

;;
;; Add two numbers without carry propagation.  Caller can use this for
;; consecutive add or subtract operations.  However, the last operation
;; before a multiply must use the routine that will normalize data.
;;

PROCFL	gwzaddq1
	ad_prolog 0,0,rbx,rsi
	mov	rsi, SRCARG			; Address of first number
	mov	rcx, SRC2ARG			; Address of second number
	mov	rbx, DESTARG			; Address of destination
	mov	al, 8				; Eight sections
uaddlp2:mov	edx, addcount1			; Count of double cache lines in a section
;	mov	r9, rsi				; Save section pointers for later carry rotate and add
uaddlp:	vmovapd	zmm0, [rcx]			; Load second number
	vaddpd	zmm0, zmm0, [rsi]		; Add in first number
	vmovapd	zmm1, [rcx+64]			; Load second number
	vaddpd	zmm1, zmm1, [rsi+64]		; Add in first number
	zstore	[rbx], zmm0			; Save result
	zstore	[rbx+64], zmm1			; Save result
	bump	rcx, 128			; Next source
	bump	rsi, 128			; Next source
	bump	rbx, 128			; Next dest
	dec	edx				; Test loop counter
	jnz	short uaddlp			; Loop if necessary
;	bump	rcx, 64				; Pad 64 bytes
;	bump	rsi, 64				; Pad 64 bytes
;	bump	rbx, 64				; Pad 64 bytes
	dec	al				; Check loop counter
	jnz	uaddlp2				; Loop if necessary
	ad_epilog 0,0,rbx,rsi
gwzaddq1 ENDP

;;
;; Add two numbers with carry propagation (eight different versions)
;;

	; Irrational, not zero-padded
PROCFL	gwzadd1
	ad_prolog 0,0,rbx,rbp,rsi,rdi,r8,r9,r10,r12,r13,r14,zmm6,zmm7,zmm8,zmm9,zmm10,zmm11,zmm12,zmm13,zmm14,zmm15
	mov	rsi, SRCARG			; Address of first number
	mov	rcx, SRC2ARG			; Address of second number
	mov	rbx, DESTARG			; Address of destination
	mov	r13, pass1blkdst		; Distance between two sections (there are 8 sections)
	shl	r13, 1				; In znorm_op_1d macro, work on every other section
	lea	r14, [2*r13+r13]		; Distance to section #4
	mov	rdi, norm_biglit_array		; rdi = pointer to big/lit flags
	mov	r12, compressed_biglits		; Get pointer to compressed biglit table
	znorm_op_1d_preload exec		; Preload constants
	mov	al, 2				; Loop counter (even sections, then odd sections)
addsec:	mov	r10d, addcount1			; Count of double cache lines in a section
	mov	r9, rbx				; Save section start pointers for later zr1_add_carry and positioning to next set of sections
	mov	r8, rdi
	push	rsi
	push	rcx
addlp:	znorm_op_1d vaddpd, exec		; Add and normalize 4 double cache lines
	bump	rsi, 128			; Next source/dest pointers
	bump	rcx, 128
	bump	rbx, 128
	bump	rdi, 1				; Next big/little flags
	dec	r10d				; Are we done with the inner loop?
	jnz	addlp				; Loop til done
	call	zi1_add_carry			; Rotate and apply carries, create carries for next four sections
	; Position to next set of 4 sections and loop
	shr	r13, 1				; Compute section distance
	lea	rbx, [r9+r13]			; Calculate next source/dest ptrs
	pop	rcx
	pop	rsi
	add	rcx, r13
	add	rsi, r13
	shl	r13, 1				; Restore every other section distance
	dec	al				; Test outer loop counter
	jnz	addsec				; Iterate
	call	zi1_add_carry_cleanup		; Handle last eight carries out of the last four sections
	ad_epilog 0,0,rbx,rbp,rsi,rdi,r8,r9,r10,r12,r13,r14,zmm6,zmm7,zmm8,zmm9,zmm10,zmm11,zmm12,zmm13,zmm14,zmm15
gwzadd1 ENDP

	; Rational, not zero-padded
PROCFL	gwzaddr1
	ad_prolog 0,0,rbx,rbp,rsi,rdi,r8,r9,r10,r12,r13,r14,zmm6,zmm7,zmm8,zmm9,zmm10,zmm11,zmm12,zmm13,zmm14,zmm15
	mov	rsi, SRCARG			; Address of first number
	mov	rcx, SRC2ARG			; Address of second number
	mov	rbx, DESTARG			; Address of destination
	mov	r13, pass1blkdst		; Distance between two sections (there are 8 sections)
	shl	r13, 1				; In znorm_op_1d macro, work on every other section
	lea	r14, [2*r13+r13]		; Distance to section #4
	znorm_op_1d_preload noexec		; Preload constants
	mov	al, 2				; Loop counter (even sections, then odd sections)
raddsec:mov	r10d, addcount1			; Count of double cache lines in a section
	mov	r9, rbx				; Save section start pointers for later zr1_add_carry and positioning to next set of sections
	push	rsi
	push	rcx
raddlp:	znorm_op_1d vaddpd, noexec		; Add and normalize 4 double cache lines
	bump	rsi, 128			; Next source/dest pointers
	bump	rcx, 128
	bump	rbx, 128
	dec	r10d				; Are we done with the inner loop?
	jnz	raddlp				; Loop til done
	call	zr1_add_carry			; Rotate and apply carries, create carries for next four sections
	;; Position to next set of 4 sections and loop
	shr	r13, 1				; Compute section distance
	lea	rbx, [r9+r13]			; Calculate next source/dest ptrs
	pop	rcx
	pop	rsi
	add	rcx, r13
	add	rsi, r13
	shl	r13, 1				; Restore every other section distance
	dec	al				; Test outer loop counter
	jnz	raddsec				; Iterate
	call	zr1_add_carry_cleanup		; Handle last eight carries out of the last four sections
	ad_epilog 0,0,rbx,rbp,rsi,rdi,r8,r9,r10,r12,r13,r14,zmm6,zmm7,zmm8,zmm9,zmm10,zmm11,zmm12,zmm13,zmm14,zmm15
gwzaddr1 ENDP

	; Irrational, zero-padded
PROCFL	gwzaddzp1
	ad_prolog 0,0,rbx,rbp,rsi,rdi,r8,r9,r10,r12,r13,r14,zmm6,zmm7,zmm8,zmm9,zmm10,zmm11,zmm12,zmm13,zmm14,zmm15
	mov	rsi, SRCARG			; Address of first number
	mov	rcx, SRC2ARG			; Address of second number
	mov	rbx, DESTARG			; Address of destination
	mov	r13, pass1blkdst		; Distance between two sections (there are 8 sections)
	shl	r13, 1				; In znorm_op_1d macro, work on every other section
	lea	r14, [2*r13+r13]		; Distance to section #4
	mov	rdi, norm_biglit_array		; rdi = pointer to big/lit flags
	mov	r12, compressed_biglits		; Get pointer to compressed biglit table
	znorm_op_1d_zpad_preload exec		; Preload constants
	mov	al, 2				; Loop counter (even sections, then odd sections)
zaddsec:mov	r10d, addcount1			; Count of double cache lines in a section
	mov	r9, rbx				; Save section start pointers for later zr1_add_carry and positioning to next set of sections
	mov	r8, rdi
	push	rsi
	push	rcx
	push	r12				; Save section start compressed biglit table pointer
zaddlp:	znorm_op_1d_zpad vaddpd, exec		; Add and normalize 4 double cache lines
	bump	rsi, 128			; Next source/dest pointers
	bump	rcx, 128
	bump	rbx, 128
	mov	rdx, r12			; If r12 has been been bumped by 4, we'll need to unbump r12 and bump rdi
	shr	rdx, 3				; Set carry if r12 is in +4 state
	adc	rdi, 0				; Increment (or not) the biglit table pointer
	xor	r12, 4				; Bump or unbump pointer into the compressed biglit table
	dec	r10d				; Are we done with the inner loop?
	jnz	zaddlp				; Loop til done
	pop	rbp				; Load section start compressed biglit table pointer
	call	zi1_add_carry_op_zpad		; Rotate and apply carries, create carries for next four sections
	; Position to next set of 4 sections and loop
	shr	r13, 1				; Compute section distance
	lea	rbx, [r9+r13]			; Calculate next source/dest ptrs
	pop	rcx
	pop	rsi
	add	rcx, r13
	add	rsi, r13
	shl	r13, 1				; Restore every other section distance
	dec	al				; Test outer loop counter
	jnz	zaddsec				; Iterate
	call	zi1_add_carry_op_zpad_cleanup	; Handle last eight carries out of the last four sections
	ad_epilog 0,0,rbx,rbp,rsi,rdi,r8,r9,r10,r12,r13,r14,zmm6,zmm7,zmm8,zmm9,zmm10,zmm11,zmm12,zmm13,zmm14,zmm15
gwzaddzp1 ENDP

	; Rational, zero-padded
PROCFL	gwzaddrzp1
	ad_prolog 0,0,rbx,rbp,rsi,rdi,r8,r9,r10,r12,r13,r14,zmm6,zmm7,zmm8,zmm9,zmm10,zmm11,zmm12,zmm13,zmm14,zmm15
	mov	rsi, SRCARG			; Address of first number
	mov	rcx, SRC2ARG			; Address of second number
	mov	rbx, DESTARG			; Address of destination
	mov	r13, pass1blkdst		; Distance between two sections (there are 8 sections)
	shl	r13, 1				; In znorm_op_1d macro, work on every other section
	lea	r14, [2*r13+r13]		; Distance to section #4
	znorm_op_1d_zpad_preload noexec		; Preload constants
	mov	al, 2				; Loop counter (even sections, then odd sections)
zraddsec:mov	r10d, addcount1			; Count of double cache lines in a section
	mov	r9, rbx				; Save section start pointers for later zr1_add_carry and positioning to next set of sections
	push	rsi
	push	rcx
zraddlp:znorm_op_1d_zpad vaddpd, noexec		; Add and normalize 4 double cache lines
	bump	rsi, 128			; Next source/dest pointers
	bump	rcx, 128
	bump	rbx, 128
	dec	r10d				; Are we done with the inner loop?
	jnz	zraddlp				; Loop til done
	call	zr1_add_carry_op_zpad		; Rotate and apply carries, create carries for next four sections
	;; Position to next set of 4 sections and loop
	shr	r13, 1				; Compute section distance
	lea	rbx, [r9+r13]			; Calculate next source/dest ptrs
	pop	rcx
	pop	rsi
	add	rcx, r13
	add	rsi, r13
	shl	r13, 1				; Restore every other section distance
	dec	al				; Test outer loop counter
	jnz	zraddsec			; Iterate
	call	zr1_add_carry_op_zpad_cleanup	; Handle last eight carries out of the last four sections
	ad_epilog 0,0,rbx,rbp,rsi,rdi,r8,r9,r10,r12,r13,r14,zmm6,zmm7,zmm8,zmm9,zmm10,zmm11,zmm12,zmm13,zmm14,zmm15
gwzaddrzp1 ENDP


;;
;; Subtract two numbers without carry propagation.  Caller can use this for
;; consecutive add or subtract operations.  However, the last operation
;; before a multiply must use the routine that will normalize data.
;;

PROCFL	gwzsubq1
	ad_prolog 0,0,rbx,rsi
	mov	rsi, SRCARG			; Address of first number
	mov	rcx, SRC2ARG			; Address of second number
	mov	rbx, DESTARG			; Address of destination
	mov	al, 8				; Eight sections
usublp2:mov	edx, addcount1			; Count of double cache lines in a section
;	mov	r9, rsi				; Save section pointers for later carry rotate and add
usublp:	vmovapd	zmm0, [rcx]			; Load second number
	vsubpd	zmm0, zmm0, [rsi]		; Add in first number
	vmovapd	zmm1, [rcx+64]			; Load second number
	vsubpd	zmm1, zmm1, [rsi+64]		; Add in first number
	zstore	[rbx], zmm0			; Save result
	zstore	[rbx+64], zmm1			; Save result
	bump	rcx, 128			; Next source
	bump	rsi, 128			; Next source
	bump	rbx, 128			; Next dest
	dec	edx				; Test loop counter
	jnz	short usublp			; Loop if necessary
;	bump	rcx, 64				; Pad 64 bytes
;	bump	rsi, 64				; Pad 64 bytes
;	bump	rbx, 64				; Pad 64 bytes
	dec	al				; Check loop counter
	jnz	usublp2				; Loop if necessary
	ad_epilog 0,0,rbx,rsi
gwzsubq1 ENDP

;;
;; Subtract two numbers with carry propagation (eight different versions)
;;

	; Irrational, not zero-padded
PROCFL	gwzsub1
	ad_prolog 0,0,rbx,rbp,rsi,rdi,r8,r9,r10,r12,r13,r14,zmm6,zmm7,zmm8,zmm9,zmm10,zmm11,zmm12,zmm13,zmm14,zmm15
	mov	rsi, SRCARG			; Address of first number
	mov	rcx, SRC2ARG			; Address of second number
	mov	rbx, DESTARG			; Address of destination
	mov	r13, pass1blkdst		; Distance between two sections (there are 8 sections)
	shl	r13, 1				; In znorm_op_1d macro, work on every other section
	lea	r14, [2*r13+r13]		; Distance to section #4
	mov	rdi, norm_biglit_array		; rdi = pointer to big/lit flags
	mov	r12, compressed_biglits		; Get pointer to compressed biglit table
	znorm_op_1d_preload exec		; Preload constants
	mov	al, 2				; Loop counter (even sections, then odd sections)
subsec:	mov	r10d, addcount1			; Count of double cache lines in a section
	mov	r9, rbx				; Save section start pointers for later zr1_add_carry and positioning to next set of sections
	mov	r8, rdi
	push	rsi
	push	rcx
sublp:	znorm_op_1d vsubpd, exec		; Add and normalize 4 double cache lines
	bump	rsi, 128			; Next source/dest pointers
	bump	rcx, 128
	bump	rbx, 128
	bump	rdi, 1				; Next big/little flags
	dec	r10d				; Are we done with the inner loop?
	jnz	sublp				; Loop til done
	call	zi1_add_carry			; Rotate and apply carries, create carries for next four sections
	; Position to next set of 4 sections and loop
	shr	r13, 1				; Compute section distance
	lea	rbx, [r9+r13]			; Calculate next source/dest ptrs
	pop	rcx
	pop	rsi
	add	rcx, r13
	add	rsi, r13
	shl	r13, 1				; Restore every other section distance
	dec	al				; Test outer loop counter
	jnz	subsec				; Iterate
	call	zi1_add_carry_cleanup		; Handle last eight carries out of the last four sections
	ad_epilog 0,0,rbx,rbp,rsi,rdi,r8,r9,r10,r12,r13,r14,zmm6,zmm7,zmm8,zmm9,zmm10,zmm11,zmm12,zmm13,zmm14,zmm15
gwzsub1 ENDP

	; Rational, not zero-padded
PROCFL	gwzsubr1
	ad_prolog 0,0,rbx,rbp,rsi,rdi,r8,r9,r10,r12,r13,r14,zmm6,zmm7,zmm8,zmm9,zmm10,zmm11,zmm12,zmm13,zmm14,zmm15
	mov	rsi, SRCARG			; Address of first number
	mov	rcx, SRC2ARG			; Address of second number
	mov	rbx, DESTARG			; Address of destination
	mov	r13, pass1blkdst		; Distance between two sections (there are 8 sections)
	shl	r13, 1				; In znorm_op_1d macro, work on every other section
	lea	r14, [2*r13+r13]		; Distance to section #4
	znorm_op_1d_preload noexec		; Preload constants
	mov	al, 2				; Loop counter (even sections, then odd sections)
rsubsec:mov	r10d, addcount1			; Count of double cache lines in a section
	mov	r9, rbx				; Save section start pointers for later zr1_add_carry and positioning to next set of sections
	push	rsi
	push	rcx
rsublp:	znorm_op_1d vsubpd, noexec		; Add and normalize 4 double cache lines
	bump	rsi, 128			; Next source/dest pointers
	bump	rcx, 128
	bump	rbx, 128
	dec	r10d				; Are we done with the inner loop?
	jnz	rsublp				; Loop til done
	call	zr1_add_carry			; Rotate and apply carries, create carries for next four sections
	; Position to next set of 4 sections and loop
	shr	r13, 1				; Compute section distance
	lea	rbx, [r9+r13]			; Calculate next source/dest ptrs
	pop	rcx
	pop	rsi
	add	rcx, r13
	add	rsi, r13
	shl	r13, 1				; Restore every other section distance
	dec	al				; Test outer loop counter
	jnz	rsubsec				; Iterate
	call	zr1_add_carry_cleanup		; Handle last eight carries out of the last four sections
	ad_epilog 0,0,rbx,rbp,rsi,rdi,r8,r9,r10,r12,r13,r14,zmm6,zmm7,zmm8,zmm9,zmm10,zmm11,zmm12,zmm13,zmm14,zmm15
gwzsubr1 ENDP

	; Irrational, zero-padded
PROCFL	gwzsubzp1
	ad_prolog 0,0,rbx,rbp,rsi,rdi,r8,r9,r10,r12,r13,r14,zmm6,zmm7,zmm8,zmm9,zmm10,zmm11,zmm12,zmm13,zmm14,zmm15
	mov	rsi, SRCARG			; Address of first number
	mov	rcx, SRC2ARG			; Address of second number
	mov	rbx, DESTARG			; Address of destination
	mov	r13, pass1blkdst		; Distance between two sections (there are 8 sections)
	shl	r13, 1				; In znorm_op_1d macro, work on every other section
	lea	r14, [2*r13+r13]		; Distance to section #4
	mov	rdi, norm_biglit_array		; rdi = pointer to big/lit flags
	mov	r12, compressed_biglits		; Get pointer to compressed biglit table
	znorm_op_1d_zpad_preload exec		; Preload constants
	mov	al, 2				; Loop counter (even sections, then odd sections)
zsubsec:mov	r10d, addcount1			; Count of double cache lines in a section
	mov	r9, rbx				; Save section start pointers for later zr1_add_carry and positioning to next set of sections
	mov	r8, rdi
	push	rsi
	push	rcx
	push	r12				; Save section start compressed biglit table pointer
zsublp:	znorm_op_1d_zpad vsubpd, exec		; Add and normalize 4 double cache lines
	bump	rsi, 128			; Next source/dest pointers
	bump	rcx, 128
	bump	rbx, 128
	mov	rdx, r12			; If r12 has been been bumped by 4, we'll need to unbump r12 and bump rdi
	shr	rdx, 3				; Set carry if r12 is in +4 state
	adc	rdi, 0				; Increment (or not) the biglit table pointer
	xor	r12, 4				; Bump or unbump pointer into the compressed biglit table
	dec	r10d				; Are we done with the inner loop?
	jnz	zsublp				; Loop til done
	pop	rbp				; Load section start compressed biglit table pointer
	call	zi1_add_carry_op_zpad		; Rotate and apply carries, create carries for next four sections
	; Position to next set of 4 sections and loop
	shr	r13, 1				; Compute section distance
	lea	rbx, [r9+r13]			; Calculate next source/dest ptrs
	pop	rcx
	pop	rsi
	add	rcx, r13
	add	rsi, r13
	shl	r13, 1				; Restore every other section distance
	dec	al				; Test outer loop counter
	jnz	zsubsec				; Iterate
	call	zi1_add_carry_op_zpad_cleanup	; Handle last eight carries out of the last four sections
	ad_epilog 0,0,rbx,rbp,rsi,rdi,r8,r9,r10,r12,r13,r14,zmm6,zmm7,zmm8,zmm9,zmm10,zmm11,zmm12,zmm13,zmm14,zmm15
gwzsubzp1 ENDP

	; Rational, zero-padded
PROCFL	gwzsubrzp1
	ad_prolog 0,0,rbx,rbp,rsi,rdi,r8,r9,r10,r12,r13,r14,zmm6,zmm7,zmm8,zmm9,zmm10,zmm11,zmm12,zmm13,zmm14,zmm15
	mov	rsi, SRCARG			; Address of first number
	mov	rcx, SRC2ARG			; Address of second number
	mov	rbx, DESTARG			; Address of destination
	mov	r13, pass1blkdst		; Distance between two sections (there are 8 sections)
	shl	r13, 1				; In znorm_op_1d macro, work on every other section
	lea	r14, [2*r13+r13]		; Distance to section #4
	znorm_op_1d_zpad_preload noexec		; Preload constants
	mov	al, 2				; Loop counter (even sections, then odd sections)
zrsubsec:mov	r10d, addcount1			; Count of double cache lines in a section
	mov	r9, rbx				; Save section start pointers for later zr1_add_carry and positioning to next set of sections
	mov	r8, rdi
	push	rsi
	push	rcx
zrsublp:znorm_op_1d_zpad vsubpd, noexec		; Add and normalize 4 double cache lines
	bump	rsi, 128			; Next source/dest pointers
	bump	rcx, 128
	bump	rbx, 128
	dec	r10d				; Are we done with the inner loop?
	jnz	zrsublp				; Loop til done
	call	zr1_add_carry_op_zpad		; Rotate and apply carries, create carries for next four sections
	; Position to next set of 4 sections and loop
	shr	r13, 1				; Compute section distance
	lea	rbx, [r9+r13]			; Calculate next source/dest ptrs
	pop	rcx
	pop	rsi
	add	rcx, r13
	add	rsi, r13
	shl	r13, 1				; Restore every other section distance
	dec	al				; Test outer loop counter
	jnz	zrsubsec			; Iterate
	call	zr1_add_carry_op_zpad_cleanup	; Handle last eight carries out of the last four sections
	ad_epilog 0,0,rbx,rbp,rsi,rdi,r8,r9,r10,r12,r13,r14,zmm6,zmm7,zmm8,zmm9,zmm10,zmm11,zmm12,zmm13,zmm14,zmm15
gwzsubrzp1 ENDP


;;
;; Add and subtract two numbers without carry propagation.
;;

PROCFL	gwzaddsubq1
	ad_prolog 0,0,rbx,rbp,rsi
	mov	rsi, SRCARG			; Address of first number
	mov	rcx, SRC2ARG			; Address of second number
	mov	rbx, DESTARG			; Address of destination
	mov	rbp, DEST2ARG			; Address of destination #2
	mov	al, 8				; Eight sections
uaslp2:	mov	edx, addcount1			; Count of double cache lines in a section
;	mov	r9, rsi				; Save section pointers for later carry rotate and add
uaslp:	vmovapd	zmm0, [rsi]			; Load first number
	vmovapd	zmm1, [rcx]			; Load second number
	vaddpd	zmm2, zmm0, zmm1		; Add in second number
	vsubpd	zmm3, zmm0, zmm1		; Subtract out second number
	vmovapd	zmm0, [rsi+64]			; Load first number
	vmovapd	zmm1, [rcx+64]			; Load second number
	vaddpd	zmm4, zmm0, zmm1		; Add in second number
	vsubpd	zmm5, zmm0, zmm1		; Subtract out second number
	zstore	[rbx], zmm2			; Save result
	zstore	[rbx+64], zmm4			; Save result
	zstore	[rbp], zmm3			; Save result
	zstore	[rbp+64], zmm5			; Save result
	bump	rcx, 128			; Next source
	bump	rsi, 128			; Next source
	bump	rbx, 128			; Next dest
	bump	rbp, 128			; Next dest
	dec	edx				; Test loop counter
	jnz	short uaslp			; Loop if necessary
;	bump	rcx, 64				; Pad 64 bytes
;	bump	rsi, 64				; Pad 64 bytes
;	bump	rbx, 64				; Pad 64 bytes
;	bump	rbp, 64				; Pad 64 bytes
	dec	al				; Check loop counter
	jnz	uaslp2				; Loop if necessary
	ad_epilog 0,0,rbx,rbp,rsi
gwzaddsubq1 ENDP

;;
;; Add and subtract two numbers with carry propagation (eight different versions)
;;

	; Irrational, not zero-padded
PROCFL	gwzaddsub1
	ad_prolog 0,0,rbx,rbp,rsi,rdi,r8,r9,r10,r12,r13,r14,r15,zmm6,zmm7,zmm8,zmm9,zmm10,zmm11,zmm12,zmm13,zmm14,zmm15
	mov	rsi, SRCARG			; Address of first number
	mov	rcx, SRC2ARG			; Address of second number
	mov	rbx, DESTARG			; Address of destination
	mov	rbp, DEST2ARG			; Address of destination #2
	mov	r13, pass1blkdst		; Distance between two sections (there are 8 sections)
	shl	r13, 1				; In znorm_op_1d macro, work on every other section
	lea	r14, [2*r13+r13]		; Distance to section #4
	mov	rdi, norm_biglit_array		; rdi = pointer to big/lit flags
	mov	r12, compressed_biglits		; Get pointer to compressed biglit table
	znorm_addsub_1d_preload exec		; Preload constants
	mov	al, 2				; Loop counter (even sections, then odd sections)
assec:	mov	r10d, addcount1			; Count of double cache lines in a section
	mov	r9, rbx				; Save section start pointers for later zr1_add_carry and positioning to next set of sections
	mov	r15, rbp
	mov	r8, rdi
	push	rsi
	push	rcx
aslp:	znorm_addsub_1d exec			; Add/sub and normalize 4 double cache lines
	bump	rsi, 128			; Next source/dest pointers
	bump	rcx, 128
	bump	rbx, 128
	bump	rbp, 128
	bump	rdi, 1				; Next big/little flags
	dec	r10d				; Are we done with the inner loop?
	jnz	aslp				; Loop til done
	znorm_addsub_save 1			; Save state of subtract carries
	call	zi1_add_carry			; Rotate and apply carries, create carries for next four sections
	znorm_addsub_save 2			; Save state of add carries, restore state of subtract carries
	call	zi1_add_carry			; Rotate and apply carries, create carries for next four sections
	znorm_addsub_save 3			; Restore state of add carries
	; Position to next set of 4 sections and loop
	shr	r13, 1				; Compute section distance
	lea	rbx, [r9+r13]			; Calculate next source/dest ptrs
	lea	rbp, [r15+r13]
	pop	rcx
	pop	rsi
	add	rcx, r13
	add	rsi, r13
	shl	r13, 1				; Restore every other section distance
	dec	al				; Test outer loop counter
	jnz	assec				; Iterate
	znorm_addsub_save 1			; Save state of subtract carries
	call	zi1_add_carry_cleanup		; Handle last eight carries out of the last four sections
	znorm_addsub_save 4			; Restore state of subtract carries
	call	zi1_add_carry_cleanup		; Handle last eight carries out of the last four sections
	ad_epilog 0,0,rbx,rbp,rsi,rdi,r8,r9,r10,r12,r13,r14,r15,zmm6,zmm7,zmm8,zmm9,zmm10,zmm11,zmm12,zmm13,zmm14,zmm15
gwzaddsub1 ENDP

	; Rational, not zero-padded
PROCFL	gwzaddsubr1
	ad_prolog 0,0,rbx,rbp,rsi,rdi,r8,r9,r10,r12,r13,r14,r15,zmm6,zmm7,zmm8,zmm9,zmm10,zmm11,zmm12,zmm13,zmm14,zmm15
	mov	rsi, SRCARG			; Address of first number
	mov	rcx, SRC2ARG			; Address of second number
	mov	rbx, DESTARG			; Address of destination
	mov	rbp, DEST2ARG			; Address of destination #2
	mov	r13, pass1blkdst		; Distance between two sections (there are 8 sections)
	shl	r13, 1				; In znorm_op_1d macro, work on every other section
	lea	r14, [2*r13+r13]		; Distance to section #4
	znorm_addsub_1d_preload noexec		; Preload constants
	mov	al, 2				; Loop counter (even sections, then odd sections)
rassec:	mov	r10d, addcount1			; Count of double cache lines in a section
	mov	r9, rbx				; Save section start pointers for later zr1_add_carry and positioning to next set of sections
	mov	r15, rbp
	push	rsi
	push	rcx
raslp:	znorm_addsub_1d noexec			; Add/sub and normalize 4 double cache lines
	bump	rsi, 128			; Next source/dest pointers
	bump	rcx, 128
	bump	rbx, 128
	bump	rbp, 128
	dec	r10d				; Are we done with the inner loop?
	jnz	raslp				; Loop til done
	znorm_addsub_save 1			; Save state of subtract carries
	call	zr1_add_carry			; Rotate and apply carries, create carries for next four sections
	znorm_addsub_save 2			; Save state of add carries, restore state of subtract carries
	call	zr1_add_carry			; Rotate and apply carries, create carries for next four sections
	znorm_addsub_save 3			; Restore state of add carries
	; Position to next set of 4 sections and loop
	shr	r13, 1				; Compute section distance
	lea	rbx, [r9+r13]			; Calculate next source/dest ptrs
	lea	rbp, [r15+r13]
	pop	rcx
	pop	rsi
	add	rcx, r13
	add	rsi, r13
	shl	r13, 1				; Restore every other section distance
	dec	al				; Test outer loop counter
	jnz	rassec				; Iterate
	znorm_addsub_save 1			; Save state of subtract carries
	call	zr1_add_carry_cleanup		; Handle last eight carries out of the last four sections
	znorm_addsub_save 4			; Restore state of subtract carries
	call	zr1_add_carry_cleanup		; Handle last eight carries out of the last four sections
	ad_epilog 0,0,rbx,rbp,rsi,rdi,r8,r9,r10,r12,r13,r14,r15,zmm6,zmm7,zmm8,zmm9,zmm10,zmm11,zmm12,zmm13,zmm14,zmm15
gwzaddsubr1 ENDP

	; Irrational, zero-padded
PROCFL	gwzaddsubzp1
	ad_prolog 0,0,rbx,rbp,rsi,rdi,r8,r9,r10,r12,r13,r14,r15,zmm6,zmm7,zmm8,zmm9,zmm10,zmm11,zmm12,zmm13,zmm14,zmm15
	mov	rsi, SRCARG			; Address of first number
	mov	rcx, SRC2ARG			; Address of second number
	mov	rbx, DESTARG			; Address of destination
	mov	rbp, DEST2ARG			; Address of destination #2
	mov	r13, pass1blkdst		; Distance between two sections (there are 8 sections)
	shl	r13, 1				; In znorm_op_1d macro, work on every other section
	lea	r14, [2*r13+r13]		; Distance to section #4
	mov	rdi, norm_biglit_array		; rdi = pointer to big/lit flags
	mov	r12, compressed_biglits		; Get pointer to compressed biglit table
	znorm_addsub_1d_zpad_preload exec	; Preload constants
	mov	al, 2				; Loop counter (even sections, then odd sections)
zassec:	mov	r10d, addcount1			; Count of double cache lines in a section
	mov	r9, rbx				; Save section start pointers for later zr1_add_carry and positioning to next set of sections
	mov	r15, rbp
	mov	r8, rdi
	push	rsi
	push	rcx
	push	r12				; Save section start compressed biglit table pointer
zaslp:	znorm_addsub_1d_zpad exec		; Add/sub and normalize 4 double cache lines
	bump	rsi, 128			; Next source/dest pointers
	bump	rcx, 128
	bump	rbx, 128
	bump	rbp, 128
	mov	rdx, r12			; If r12 has been been bumped by 4, we'll need to unbump r12 and bump rdi
	shr	rdx, 3				; Set carry if r12 is in +4 state
	adc	rdi, 0				; Increment (or not) the biglit table pointer
	xor	r12, 4				; Bump or unbump pointer into the compressed biglit table
	dec	r10d				; Are we done with the inner loop?
	jnz	zaslp				; Loop til done
	pop	rbp				; Load section start compressed biglit table pointer
	znorm_addsub_1d_zpad_save 1		; Save state of subtract carries
	call	zi1_add_carry_op_zpad		; Rotate and apply carries, create carries for next four sections
	znorm_addsub_1d_zpad_save 2		; Save state of add carries, restore state of subtract carries
	call	zi1_add_carry_op_zpad		; Rotate and apply carries, create carries for next four sections
	znorm_addsub_1d_zpad_save 3		; Restore state of add carries
	; Position to next set of 4 sections and loop
	shr	r13, 1				; Compute section distance
	lea	rbx, [r9+r13]			; Calculate next source/dest ptrs
	lea	rbp, [r15+r13]
	pop	rcx
	pop	rsi
	add	rcx, r13
	add	rsi, r13
	shl	r13, 1				; Restore every other section distance
	dec	al				; Test outer loop counter
	jnz	zassec				; Iterate
	znorm_addsub_1d_zpad_save 1		; Save state of subtract carries
	call	zi1_add_carry_op_zpad_cleanup	; Handle last eight carries out of the last four sections
	znorm_addsub_1d_zpad_save 4		; Restore state of subtract carries
	call	zi1_add_carry_op_zpad_cleanup	; Handle last eight carries out of the last four sections
	ad_epilog 0,0,rbx,rbp,rsi,rdi,r8,r9,r10,r12,r13,r14,r15,zmm6,zmm7,zmm8,zmm9,zmm10,zmm11,zmm12,zmm13,zmm14,zmm15
gwzaddsubzp1 ENDP

	; Rational, zero-padded
PROCFL	gwzaddsubrzp1
	ad_prolog 0,0,rbx,rbp,rsi,rdi,r8,r9,r10,r12,r13,r14,r15,zmm6,zmm7,zmm8,zmm9,zmm10,zmm11,zmm12,zmm13,zmm14,zmm15
	mov	rsi, SRCARG			; Address of first number
	mov	rcx, SRC2ARG			; Address of second number
	mov	rbx, DESTARG			; Address of destination
	mov	rbp, DEST2ARG			; Address of destination #2
	mov	r13, pass1blkdst		; Distance between two sections (there are 8 sections)
	shl	r13, 1				; In znorm_op_1d macro, work on every other section
	lea	r14, [2*r13+r13]		; Distance to section #4
	znorm_addsub_1d_zpad_preload noexec	; Preload constants
	mov	al, 2				; Loop counter (even sections, then odd sections)
zrassec:mov	r10d, addcount1			; Count of double cache lines in a section
	mov	r9, rbx				; Save section start pointers for later zr1_add_carry and positioning to next set of sections
	mov	r15, rbp
	push	rsi
	push	rcx
zraslp:	znorm_addsub_1d_zpad noexec		; Add/sub and normalize 4 double cache lines
	bump	rsi, 128			; Next source/dest pointers
	bump	rcx, 128
	bump	rbx, 128
	bump	rbp, 128
	dec	r10d				; Are we done with the inner loop?
	jnz	zraslp				; Loop til done
	znorm_addsub_1d_zpad_save 1		; Save state of subtract carries
	call	zr1_add_carry_op_zpad		; Rotate and apply carries, create carries for next four sections
	znorm_addsub_1d_zpad_save 2		; Save state of add carries, restore state of subtract carries
	call	zr1_add_carry_op_zpad		; Rotate and apply carries, create carries for next four sections
	znorm_addsub_1d_zpad_save 3		; Restore state of add carries
	; Position to next set of 4 sections and loop
	shr	r13, 1				; Compute section distance
	lea	rbx, [r9+r13]			; Calculate next source/dest ptrs
	lea	rbp, [r15+r13]
	pop	rcx
	pop	rsi
	add	rcx, r13
	add	rsi, r13
	shl	r13, 1				; Restore every other section distance
	dec	al				; Test outer loop counter
	jnz	zrassec				; Iterate
	znorm_addsub_1d_zpad_save 1		; Save state of subtract carries
	call	zr1_add_carry_op_zpad_cleanup	; Handle last eight carries out of the last four sections
	znorm_addsub_1d_zpad_save 4		; Restore state of subtract carries
	call	zr1_add_carry_op_zpad_cleanup	; Handle last eight carries out of the last four sections
	ad_epilog 0,0,rbx,rbp,rsi,rdi,r8,r9,r10,r12,r13,r14,r15,zmm6,zmm7,zmm8,zmm9,zmm10,zmm11,zmm12,zmm13,zmm14,zmm15
gwzaddsubrzp1 ENDP


;;
;; Copy 4KB a number with optional masking
;;

PROCFL	gwzcopy4kb
	ad_prolog 0,0,rsi,rdi
	mov	rsi, SRCARG			; Address of first number
	mov	rcx, SRC2ARG			; Address of mask
	mov	rdi, DESTARG			; Address of destination
	mov	al, 4096/256			; Count of 256 byte chunks in 4KB
z4klp:	vmovapd	zmm0, [rsi]
	vandpd	zmm0, zmm0, [rcx]
	vmovapd	zmm1, [rsi+64]
	vandpd	zmm1, zmm1, [rcx+64]
	vmovapd	zmm2, [rsi+128]
	vandpd	zmm2, zmm2, [rcx+128]
	vmovapd	zmm3, [rsi+192]
	vandpd	zmm3, zmm3, [rcx+192]
	zstore	[rdi], zmm0			; Store destination data
	zstore	[rdi+64], zmm1
	zstore	[rdi+128], zmm2
	zstore	[rdi+192], zmm3
	bump	rsi, 256			; Next src ptr
	bump	rcx, 256			; Next mask ptr
	bump	rdi, 256			; Next dest ptr
	dec	al				; Decrement count
	jnz	short z4klp
	ad_epilog 0,0,rsi,rdi
gwzcopy4kb ENDP

;;
;; Multiply a number by a small value (four versions)
;;

	; Irrational, not zero-padded
PROCFL	gwzmuls1
	ad_prolog 0,0,rbx,rbp,rsi,rdi,r8,r9,r10,r12,r13,r14,zmm6,zmm7,zmm8,zmm9,zmm10,zmm11,zmm12,zmm13,zmm14,zmm15
	mov	rsi, SRCARG			; Address of source number
	mov	rbx, DESTARG			; Address of destination
	vbroadcastsd zmm31, DBLARG		; Load small multiplier value
	mov	r13, pass1blkdst		; Distance between two sections (there are 8 sections)
	shl	r13, 1				; In znorm_op_1d macro, work on every other section
	lea	r14, [2*r13+r13]		; Distance to section #4
	mov	rdi, norm_biglit_array		; rdi = pointer to big/lit flags
	mov	r12, compressed_biglits		; Get pointer to compressed biglit table
	znorm_smallmul_1d_preload exec		; Preload constants
	mov	al, 2				; Loop counter (even sections, then odd sections)
mulsec:	mov	r10d, addcount1			; Count of double cache lines in a section
	mov	r9, rbx				; Save section start pointers for later zr1_add_carry and positioning to next set of sections
	mov	r8, rdi
	push	rsi
mullp:	znorm_smallmul_1d exec		; Add and normalize 4 double cache lines
	bump	rsi, 128			; Next source/dest pointers
	bump	rbx, 128
	bump	rdi, 1				; Next big/little flags
	dec	r10d				; Are we done with the inner loop?
	jnz	mullp				; Loop til done
	call	zi1_add_carry			; Rotate and apply carries, create carries for next four sections
	; Position to next set of 4 sections and loop
	shr	r13, 1				; Compute section distance
	lea	rbx, [r9+r13]			; Calculate next source/dest ptrs
	pop	rsi
	add	rsi, r13
	shl	r13, 1				; Restore every other section distance
	dec	al				; Test outer loop counter
	jnz	mulsec				; Iterate
	call	zi1_add_carry_cleanup		; Handle last eight carries out of the last four sections
	ad_epilog 0,0,rbx,rbp,rsi,rdi,r8,r9,r10,r12,r13,r14,zmm6,zmm7,zmm8,zmm9,zmm10,zmm11,zmm12,zmm13,zmm14,zmm15
gwzmuls1 ENDP

	; Rational, not zero-padded
PROCFL	gwzmulsr1
	ad_prolog 0,0,rbx,rbp,rsi,rdi,r8,r9,r10,r12,r13,r14,zmm6,zmm7,zmm8,zmm9,zmm10,zmm11,zmm12,zmm13,zmm14,zmm15
	mov	rsi, SRCARG			; Address of source number
	mov	rbx, DESTARG			; Address of destination
	vbroadcastsd zmm31, DBLARG		; Load small multiplier value
	mov	r13, pass1blkdst		; Distance between two sections (there are 8 sections)
	shl	r13, 1				; In znorm_op_1d macro, work on every other section
	lea	r14, [2*r13+r13]		; Distance to section #4
	znorm_smallmul_1d_preload noexec	; Preload constants
	mov	al, 2				; Loop counter (even sections, then odd sections)
rmulsec:mov	r10d, addcount1			; Count of double cache lines in a section
	mov	r9, rbx				; Save section start pointers for later zr1_add_carry and positioning to next set of sections
	push	rsi
rmullp:	znorm_smallmul_1d noexec		; Add and normalize 4 double cache lines
	bump	rsi, 128			; Next source/dest pointers
	bump	rbx, 128
	dec	r10d				; Are we done with the inner loop?
	jnz	rmullp				; Loop til done
	call	zr1_add_carry			; Rotate and apply carries, create carries for next four sections
	; Position to next set of 4 sections and loop
	shr	r13, 1				; Compute section distance
	lea	rbx, [r9+r13]			; Calculate next source/dest ptrs
	pop	rsi
	add	rsi, r13
	shl	r13, 1				; Restore every other section distance
	dec	al				; Test outer loop counter
	jnz	rmulsec				; Iterate
	call	zr1_add_carry_cleanup		; Handle last eight carries out of the last four sections
	ad_epilog 0,0,rbx,rbp,rsi,rdi,r8,r9,r10,r12,r13,r14,zmm6,zmm7,zmm8,zmm9,zmm10,zmm11,zmm12,zmm13,zmm14,zmm15
gwzmulsr1 ENDP

	; Irrational, zero-padded
PROCFL	gwzmulszp1
	ad_prolog 0,0,rbx,rbp,rsi,rdi,r8,r9,r10,r12,r13,r14,zmm6,zmm7,zmm8,zmm9,zmm10,zmm11,zmm12,zmm13,zmm14,zmm15
	mov	rsi, SRCARG			; Address of source number
	mov	rbx, DESTARG			; Address of destination
	vbroadcastsd zmm31, DBLARG		; Load small multiplier value
	mov	r13, pass1blkdst		; Distance between two sections (there are 8 sections)
	shl	r13, 1				; In znorm_op_1d macro, work on every other section
	lea	r14, [2*r13+r13]		; Distance to section #4
	mov	rdi, norm_biglit_array		; rdi = pointer to big/lit flags
	mov	r12, compressed_biglits		; Get pointer to compressed biglit table
	znorm_smallmul_1d_zpad_preload exec	; Preload constants
	mov	al, 2				; Loop counter (even sections, then odd sections)
zmulsec:mov	r10d, addcount1			; Count of double cache lines in a section
	mov	r9, rbx				; Save section start pointers for later zr1_add_carry and positioning to next set of sections
	mov	r8, rdi
	push	rsi
	push	r12				; Save section start compressed biglit table pointer
zmullp:	znorm_smallmul_1d_zpad exec		; Add and normalize 4 double cache lines
	bump	rsi, 128			; Next source/dest pointers
	bump	rbx, 128
	mov	rdx, r12			; If r12 has been been bumped by 4, we'll need to unbump r12 and bump rdi
	shr	rdx, 3				; Set carry if r12 is in +4 state
	adc	rdi, 0				; Increment (or not) the biglit table pointer
	xor	r12, 4				; Bump or unbump pointer into the compressed biglit table
	dec	r10d				; Are we done with the inner loop?
	jnz	zmullp				; Loop til done
	pop	rbp				; Load section start compressed biglit table pointer
	call	zi1_add_carry_op_zpad		; Rotate and apply carries, create carries for next four sections
	; Position to next set of 4 sections and loop
	shr	r13, 1				; Compute section distance
	lea	rbx, [r9+r13]			; Calculate next source/dest ptrs
	pop	rsi
	add	rsi, r13
	shl	r13, 1				; Restore every other section distance
	dec	al				; Test outer loop counter
	jnz	zmulsec				; Iterate
	call	zi1_add_carry_op_zpad_cleanup	; Handle last eight carries out of the last four sections
	ad_epilog 0,0,rbx,rbp,rsi,rdi,r8,r9,r10,r12,r13,r14,zmm6,zmm7,zmm8,zmm9,zmm10,zmm11,zmm12,zmm13,zmm14,zmm15
gwzmulszp1 ENDP

	; Rational, zero-padded
PROCFL	gwzmulsrzp1
	ad_prolog 0,0,rbx,rbp,rsi,rdi,r8,r9,r10,r12,r13,r14,zmm6,zmm7,zmm8,zmm9,zmm10,zmm11,zmm12,zmm13,zmm14,zmm15
	mov	rsi, SRCARG			; Address of source number
	mov	rbx, DESTARG			; Address of destination
	vbroadcastsd zmm31, DBLARG		; Load small multiplier value
	mov	r13, pass1blkdst		; Distance between two sections (there are 8 sections)
	shl	r13, 1				; In znorm_op_1d macro, work on every other section
	lea	r14, [2*r13+r13]		; Distance to section #4
	znorm_smallmul_1d_zpad_preload noexec	; Preload constants
	mov	al, 2				; Loop counter (even sections, then odd sections)
zrmulsec:mov	r10d, addcount1			; Count of double cache lines in a section
	mov	r9, rbx				; Save section start pointers for later zr1_add_carry and positioning to next set of sections
	push	rsi
zrmullp:znorm_smallmul_1d_zpad noexec		; Add and normalize 4 double cache lines
	bump	rsi, 128			; Next source/dest pointers
	bump	rbx, 128
	dec	r10d				; Are we done with the inner loop?
	jnz	zrmullp				; Loop til done
	call	zr1_add_carry_op_zpad		; Rotate and apply carries, create carries for next four sections
	; Position to next set of 4 sections and loop
	shr	r13, 1				; Compute section distance
	lea	rbx, [r9+r13]			; Calculate next source/dest ptrs
	pop	rsi
	add	rsi, r13
	shl	r13, 1				; Restore every other section distance
	dec	al				; Test outer loop counter
	jnz	zrmulsec			; Iterate
	call	zr1_add_carry_op_zpad_cleanup	; Handle last eight carries out of the last four sections
	ad_epilog 0,0,rbx,rbp,rsi,rdi,r8,r9,r10,r12,r13,r14,zmm6,zmm7,zmm8,zmm9,zmm10,zmm11,zmm12,zmm13,zmm14,zmm15
gwzmulsrzp1 ENDP

;;
;; Routines to do the normalization after a multiply
;;

;; Six internal routines to process carries after one set of four sections.
;; Rotate and apply carries, create carries for next set of four sections.
;; See zadd_carry_1d macro for register usage.
PROCFL	zr1_add_carry
	int_prolog 0,0,0
	zadd_carry_1d noexec
	int_epilog 0,0,0
zr1_add_carry ENDP
PROCFL	zi1_add_carry
	int_prolog 0,0,0
	zadd_carry_1d exec
	int_epilog 0,0,0
zi1_add_carry ENDP
PROCFL	zr1_add_carry_zpad
	int_prolog 0,0,0
	zadd_carry_1d_zpad noexec
	int_epilog 0,0,0
zr1_add_carry_zpad ENDP
PROCFL	zi1_add_carry_zpad
	int_prolog 0,0,0
	zadd_carry_1d_zpad exec
	int_epilog 0,0,0
zi1_add_carry_zpad ENDP
PROCFL	zr1_add_carry_op_zpad
	int_prolog 0,0,0
	zadd_carry_op_1d_zpad noexec
	int_epilog 0,0,0
zr1_add_carry_op_zpad ENDP
PROCFL	zi1_add_carry_op_zpad
	int_prolog 0,0,0
	zadd_carry_op_1d_zpad exec
	int_epilog 0,0,0
zi1_add_carry_op_zpad ENDP

;; Six internal routines to handle last eight carries out of the last four sections
PROCFL	zr1_add_carry_cleanup
	int_prolog 0,0,0
	zadd_carry_1d_cleanup noexec	;; Process last 8 carries
	int_epilog 0,0,0
zr1_add_carry_cleanup ENDP
PROCFL	zi1_add_carry_cleanup
	int_prolog 0,0,0
	znorm_top_carry_1d		;; Adjust carry in xmm7 if k > 1
	zadd_carry_1d_cleanup exec	;; Process last 8 carries
	int_epilog 0,0,0
zi1_add_carry_cleanup ENDP
PROCFL	zr1_add_carry_zpad_cleanup
	int_prolog 0,0,0
	zadd_carry_1d_zpad_cleanup noexec ;; Process last 8 carries
	zpad_prep_final_div_by_k_after_multiply noexec
	zpad_final_div_by_k noexec
	int_epilog 0,0,0
zr1_add_carry_zpad_cleanup ENDP
PROCFL	zi1_add_carry_zpad_cleanup
	int_prolog 0,0,0
	zadd_carry_1d_zpad_cleanup exec	;; Process last 8 carries
	zpad_prep_final_div_by_k_after_multiply exec
	zpad_final_div_by_k exec
	int_epilog 0,0,0
zi1_add_carry_zpad_cleanup ENDP
PROCFL	zr1_add_carry_op_zpad_cleanup
	int_prolog 0,0,0
	zadd_carry_op_1d_zpad_cleanup noexec ;; Process last 8 carries
	zpad_prep_final_div_by_k_after_op noexec
	zpad_final_div_by_k noexec
	int_epilog 0,0,0
zr1_add_carry_op_zpad_cleanup ENDP
PROCFL	zi1_add_carry_op_zpad_cleanup
	int_prolog 0,0,0
	zadd_carry_op_1d_zpad_cleanup exec ;; Process last 8 carries
	zpad_prep_final_div_by_k_after_op exec
	zpad_final_div_by_k exec
	int_epilog 0,0,0
zi1_add_carry_op_zpad_cleanup ENDP

ENDIF 

; Macro to loop through all the FFT values and apply the proper normalization routine.

;; On input (see pass1_normalize macro in zmult.mac):
;; rsi = FFT data (DESTARG)

;; During computation registers used thusly:
;; rax = outer loop counter
;; rcx = inner loop counter
;; rdx = register for loading compressed biglit index
;; rsi = ptr to first FFT section data
;; r13 = distance to second FFT section data
;; r14 = distance to fourth FFT section data
;; rdi = big/lit/fudge flags array
;; r12 = compressed biglit table pointer
;; r10 = pointer to inverse weight multipliers
;; r9 = saved dest ptr at section start
;; r8 = saved big/lit ptr at section start
;; rbx = available
;; rbp = available

inorm	MACRO	lab, ttp, echk, const
	LOCAL	ilp1, ilp2, ilp2dq, ilp2dn
	PROCFLP	lab
	int_prolog 0,0,0

	mov	rsi, DESTARG		;; Destination pointer
	mov	edi, ADDIN_OFFSET	;; Get address to add value into
	vmovsd	xmm0, ADDIN_VALUE	;; Get the addin value
	vaddsd	xmm0, xmm0, Q [rsi][rdi] ;; Add in the FFT value
	vmovsd	Q [rsi][rdi], xmm0	;; Save the new value

echk	vbroadcastsd zmm31, MAXERR	;; Load maximum error
	vbroadcastsd zmm30, ZMM_RNDVAL	;; Load rounding constant
	znorm_1d_preload ttp, echk, const

	mov	r13, pass1blkdst	;; Distance between two sections (there are 8 sections)
	mov	r10, norm_grp_mults	;; Inverse weights
ttp	mov	rdi, norm_biglit_array	;; rdi = pointer to big/lit flags
ttp	mov	r12, compressed_biglits	;; Get pointer to compressed biglit table
	shl	r13, 1			;; In znorm_1d macro, work on every other section
	lea	r14, [2*r13+r13]	;; Calc distance to 4th section

	vmovapd zmm0, zmm30		;; Init carries
	vmovapd zmm1, zmm30
	vmovapd zmm2, zmm30
	vmovapd zmm3, zmm30
	vmovapd zmm4, zmm30
	vmovapd zmm5, zmm30
	vmovapd zmm6, zmm30
	vmovapd zmm7, zmm30

	mov	al, 2			;; Loop counter (even sections, then odd sections)
ilp1:	mov	ecx, addcount1		;; Count of double cache lines in a section

	mov	r9, rsi			;; Save section pointers for later carry rotate and add
ttp	mov	r8, rdi

const ttp vsubpd zmm0, zmm0, zmm30	;; Remove RNDVAL from irrational mul-by-small-const normalization
const ttp vsubpd zmm1, zmm1, zmm30
const ttp vsubpd zmm2, zmm2, zmm30
const ttp vsubpd zmm3, zmm3, zmm30
const ttp vsubpd zmm4, zmm4, zmm30
const ttp vsubpd zmm5, zmm5, zmm30
const ttp vsubpd zmm6, zmm6, zmm30
const ttp vsubpd zmm7, zmm7, zmm30

;; Error checking uses AVX512DQ instructions available on most CPUs.  This section is for the rare CPU that does not support AVX512DQ.
	AVX512DQ = 0			;; Disable use of AVX512DQ instructions
echk	test	CPU_FLAGS, 0400000h	;; Test the AVX512DQ flag
echk	jnz	ilp2			;; Go use the AVX512DQ instructions
ilp2dq:					;; Copy the ilp2 loop
echk		znorm_1d ttp, echk, const	;; Normalize 64 values (4 simultaneous sections)
echk		bump	rsi, 128		;; Next source pointer
echk ttp	bump	r10, 4*128		;; Next set of inverse weights
echk no ttp	bump	r10, 4*64		;; Rational negacyclics also need inverse weights due to delayed mul-by-sine
echk ttp	bump	rdi, 1			;; Next big/little flags
echk		dec	ecx			;; Are we done with the inner loop?
echk		jnz	ilp2dq			;; Loop til done
echk	jmp	ilp2dn

	AVX512DQ = 1			;; Enable use of AVX512DQ instructions
ilp2:	znorm_1d ttp, echk, const	;; Normalize 64 values (4 simultaneous sections)
	bump	rsi, 128		;; Next source pointer
ttp	bump	r10, 4*128		;; Next set of inverse weights
no ttp	bump	r10, 4*64		;; Rational negacyclics also need inverse weights due to delayed mul-by-sine
ttp	bump	rdi, 1			;; Next big/little flags
	dec	ecx			;; Are we done with the inner loop?
	jnz	ilp2			;; Loop til done
ilp2dn:

const ttp vaddpd zmm0, zmm0, zmm30	;; Reapply RNDVAL for irrational mul-by-small-const normalization
const ttp vaddpd zmm1, zmm1, zmm30
const ttp vaddpd zmm2, zmm2, zmm30
const ttp vaddpd zmm3, zmm3, zmm30
const ttp vaddpd zmm4, zmm4, zmm30
const ttp vaddpd zmm5, zmm5, zmm30
const ttp vaddpd zmm6, zmm6, zmm30
const ttp vaddpd zmm7, zmm7, zmm30

	;; Rotate and apply four section carries, create carries for next four sections
no ttp	call	zr1_add_carry
ttp	call	zi1_add_carry

	;; Position to next 4 sections and loop
	shr	r13, 1			;; Compute section distance
	lea	rsi, [r9+r13]		;; Calculate next source ptr
	shl	r13, 1			;; Restore every other section distance
	dec	al			;; Test outer loop counter
	jnz	ilp1			;; Iterate

	;; Handle last eight carries out of the last four sections
no ttp	call	zr1_add_carry_cleanup
ttp	call	zi1_add_carry_cleanup

	;; Apply post-addin
	mov	rsi, DESTARG		;; Restore destination address
	mov	edi, ADDIN_OFFSET	;; Get address to add value into
	vmovsd	xmm8, POSTADDIN_VALUE	;; Get the addin value
	vaddsd	xmm8, xmm8, Q [rsi][rdi] ;; Add in the FFT value
	vmovsd	Q [rsi][rdi], xmm8	;; Save the new value

	;; Compress 8 MAXERR values into one
echk	zcollapse_maxerr		;; Collapse and save MAXERR

	int_epilog 0,0,0
	ENDPP lab
	ENDM


zpnorm	MACRO	lab, ttp, echk, const, khi
	LOCAL	ilp1, ilp2
	PROCFLP	lab
	int_prolog 0,0,4

	c_call	ZPAD_SUB7		;; Subtract 7 ZPAD words from lowest FFT words

	mov	rsi, DESTARG		;; Source pointer
echk	vbroadcastsd zmm31, MAXERR	;; Load maximum error
	znorm_1d_zpad_preload ttp, echk, const, khi

	mov	r13, pass1blkdst	;; Distance between two sections (there are 8 sections)
	mov	r10, norm_grp_mults	;; Inverse weights
ttp	mov	rdi, norm_biglit_array	;; rdi = pointer to big/lit flags
ttp	mov	r12, compressed_biglits	;; Get pointer to compressed biglit table
	shl	r13, 1			;; In znorm_1d_zpad macro, work on every other section
	lea	r14, [2*r13+r13]	;; Calc distance to 4th section

	vmovapd zmm0, zmm30		;; Init carries
	vmovapd zmm1, zmm30
	vmovapd zmm2, zmm30
	vmovapd zmm3, zmm30
no ttp	vmovapd zmm4, zmm30
ttp	vxorpd	zmm4, zmm4, zmm4	;; Remove RNDVAL from high carries for irrational normalization
	vmovapd zmm5, zmm4
	vmovapd zmm6, zmm4
	vmovapd zmm7, zmm4

	mov	al, 2			;; Loop counter (even sections, then odd sections)
ilp1:	mov	ecx, addcount1		;; Count of double cache lines in a section

	mov	r9, rsi			;; Save section pointers for later carry rotate and add
ttp	mov	r8, rdi
ttp	mov	rbp, r12

ilp2:	znorm_1d_zpad ttp, echk, const, khi ;; Normalize 64 values (4 simultaneous sections)
	bump	rsi, 128		;; Next source pointer
	bump	r10, 4*64		;; Rational negacyclics also need inverse weights due to delayed mul-by-sine
ttp	mov	rdx, r12		;; If r12 has been been bumped by 4, we'll need to unbump r12 and bump rdi
ttp	shr	rdx, 3			;; Set carry if r12 is in the +4 state
ttp	adc	rdi, 0			;; Increment (or not) the biglit table pointer
ttp	xor	r12, 4			;; Bump or unbump pointer into the compressed biglit table
	dec	ecx			;; Are we done with the inner loop?
	jnz	ilp2			;; Loop til done

	;; Rotate and apply four section carries, create carries for next four sections
no ttp	call	zr1_add_carry_zpad
ttp	call	zi1_add_carry_zpad

	;; Position to next 4 sections and loop
	shr	r13, 1			;; Compute section distance
	lea	rsi, [r9+r13]		;; Calculate next source ptr
	shl	r13, 1			;; Restore every other section distance
	dec	al			;; Test outer loop counter
	jnz	ilp1			;; Iterate

	;; Handle last eight carries out of the last four sections
no ttp	call	zr1_add_carry_zpad_cleanup
ttp	call	zi1_add_carry_zpad_cleanup

	;; Compress 8 MAXERR values into one
echk	zcollapse_maxerr		;; Collapse and save MAXERR

	int_epilog 0,0,4
	ENDPP lab
	ENDM


; The many different normalization routines.  One for each valid combination of rational/irrational, error check/no error check, mul by const/no mul by const.

	inorm	zr1, noexec, noexec, noexec
	inorm	zr1e, noexec, exec, noexec
	inorm	zr1c, noexec, noexec, exec
	inorm	zr1ec, noexec, exec, exec
	inorm	zi1, exec, noexec, noexec
	inorm	zi1e, exec, exec, noexec
	inorm	zi1c, exec, noexec, exec
	inorm	zi1ec, exec, exec, exec
	zpnorm	zr1zp, noexec, noexec, noexec, exec
	zpnorm	zr1zpe, noexec, exec, noexec, exec
	zpnorm	zr1zpc, noexec, noexec, exec, exec
	zpnorm	zr1zpec, noexec, exec, exec, exec
	zpnorm	zi1zp, exec, noexec, noexec, exec
	zpnorm	zi1zpe, exec, exec, noexec, exec
	zpnorm	zi1zpc, exec, noexec, exec, exec
	zpnorm	zi1zpec, exec, exec, exec, exec
	zpnorm	zr1zpk, noexec, noexec, noexec, noexec
	zpnorm	zr1zpek, noexec, exec, noexec, noexec
	zpnorm	zr1zpck, noexec, noexec, exec, noexec
	zpnorm	zr1zpeck, noexec, exec, exec, noexec
	zpnorm	zi1zpk, exec, noexec, noexec, noexec
	zpnorm	zi1zpek, exec, exec, noexec, noexec
	zpnorm	zi1zpck, exec, noexec, exec, noexec
	zpnorm	zi1zpeck, exec, exec, exec, noexec

_TEXT	ENDS
END
