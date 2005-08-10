; Copyright 1995-2005 Just For Fun Software, Inc., all rights reserved
; Author:  George Woltman
; Email: woltman@alum.mit.edu
;
; This file implements helper routines for the ECM code
;

	TITLE   setup

IFNDEF X86_64
	.386
	.MODEL	FLAT
ENDIF

INCLUDE unravel.mac

EXTRN	_SRCARG:EXTPTR
EXTRN	_SRC2ARG:EXTPTR
EXTRN	_DESTARG:EXTPTR
EXTRN	_NUMARG:DWORD
EXTRN	_NUM2ARG:DWORD
EXTRN	_CARRYH:DWORD
EXTRN	_CARRYL:DWORD
EXTRN	_RES:DWORD
EXTRN	_EGCD_A:DWORD
EXTRN	_EGCD_B:DWORD
EXTRN	_EGCD_C:DWORD
EXTRN	_EGCD_D:DWORD
EXTRN	_EGCD_ODD:DWORD

;
; Global variables
;

_DATA SEGMENT
fltval		DD	0
fltvalhi	DD	0
fltval1		DD	0
fltval1hi	DD	0
TWOPOW32	DD	4294967296.0		; 2^32
BIGVAL		DD	5F000000h		; 2^63
FPCNTRL		DW	0
_DATA	ENDS

_TEXT	SEGMENT

; Add integer to carryh,carryl,res

	PUBLIC	_eaddhlp
_eaddhlp PROC
	mov	eax, _NUMARG		; Integer1
	add	_RES, eax		; Add two integers
	adc	_CARRYL, 0		; Add the carry
	adc	_CARRYH, 0		; Add the carry
	ret
_eaddhlp ENDP

; Subtract integer from carryh,carryl,res

	PUBLIC	_esubhlp
_esubhlp PROC
	mov	eax, _NUMARG		; Integer1
	sub	_RES, eax		; Subtract two integers
	sbb	_CARRYL, 0		; Subtract the carry
	sbb	_CARRYH, 0		; Subtract the carry
	ret
_esubhlp ENDP

; Multiply integer1 and integer2 adding result to carryh,carryl,res

	PUBLIC	_emuladdhlp
_emuladdhlp PROC
	mov	eax, _NUMARG		; Integer1
	mul	_NUM2ARG		; Integer2
	add	_RES, eax		; Add result to 3 word accumulator
	adc	_CARRYL, edx
	adc	_CARRYH, 0
	ret
_emuladdhlp ENDP

; Multiply integer1 and integer2 adding twice the result to carryh,carryl,res

	PUBLIC	_emuladd2hlp
_emuladd2hlp PROC
	mov	eax, _NUMARG		; Integer1
	mul	_NUM2ARG		; Integer2
	add	eax, eax		; Add 2*result to 3 word accumulator
	adc	edx, edx
	adc	_CARRYH, 0
	add	_RES, eax
	adc	_CARRYL, edx
	adc	_CARRYH, 0
	ret
_emuladd2hlp ENDP

; Multiply integer1 and integer2 subtracting result from carryh,carryl,res

	PUBLIC	_emulsubhlp
_emulsubhlp PROC
	mov	eax, _NUMARG		; Integer1
	mul	_NUM2ARG		; Integer2
	sub	_RES, eax		; Sub result from 3 word accumulator
	sbb	_CARRYL, edx
	sbb	_CARRYH, 0
	ret
_emulsubhlp ENDP

; Routine to help in computing extended GCD quickly
;
; Do several single-precision steps for the extended GCD code
; U is larger than V.  Both are the same length or U is one larger
; than V.  Returns A,B,C,D as defined in Knuth vol 2. description of
; extended GCD for large numbers.  This was implemented in assembly
; language to better deal with 64-bit integers.

	PUBLIC	_egcdhlp
_egcdhlp PROC
IFNDEF X86_64
	fninit
	push	ebp
	push	ebx
	push	edi
	push	esi

; Load up to 64 bits of U and V
;	U will be in edx:eax with bits shifted in from edi
;	V will be in esi:ebp with bits shifted in from ebx

	mov	edi, _SRCARG		; Giant U
	mov	ebx, _NUMARG		; U->sign
	mov	ecx, _SRC2ARG		; Giant V
	mov	ebp, _NUM2ARG		; V->sign

	mov	edx, [edi-4][ebx*4]	; U[Ulen-1]

	xor	esi, esi
	cmp	ebx, ebp
	jne	short noload
	mov	esi, [ecx-4][ebx*4]	; V[Ulen-1]

noload:	cmp	ebx, 1			; Are there more words to shift
	jg	short multi		; bits from?
	xor	eax, eax		; No, zero out MSWs
	xor	ebp, ebp
	xchg	eax, edx
	xchg	esi, ebp
	jmp	short noshft

multi:	mov	eax, [edi-8][ebx*4]	; U[Ulen-2]
	mov	ebp, [ecx-8][ebx*4]	; V[Ulen-2]
	
	cmp	ebx, 2			; Are there more words to shift
	je	short noshft		; bits from?

	mov	edi, [edi-12][ebx*4]	; U[Ulen-3]
	mov	ebx, [ecx-12][ebx*4]	; V[Ulen-3]

	bsr	ecx, edx		; Count bits to shift U
	xor	ecx, 31			; Turn bit # into a shift count

	shld	edx, eax, cl		; Shift U
	shld	eax, edi, cl
	shld	esi, ebp, cl		; Shift V
	shld	ebp, ebx, cl

; Init extended GCD information
;	EGCD_A = 1;
;	EGCD_C = 0;
;	EGCD_B = 0;
;	EGCD_D = 1;
;	EGCD_ODD = 0;

noshft:	fld1				; A
	fldz				; C
	fldz				; B
	fld1				; D

; Load the values into the FPU
;	U is in edx:eax
;	V is in esi:ebp

	mov	fltval, edx
	mov	fltval1, esi
	fild	QWORD PTR fltval	; uhi
	fild	QWORD PTR fltval1	; vhi, uhi
	mov	fltval, eax
	mov	fltval1, ebp
	fmul	TWOPOW32
	fxch	st(1)			; uhi, vhi
	fmul	TWOPOW32
	fild	QWORD PTR fltval	; ulo, uhi, vhi
	fild	QWORD PTR fltval1	; vlo, ulo, uhi, vhi
	faddp	st(3), st		; ulo, uhi, v
	faddp	st(1), st		; u, v

; Turn on truncating mode

	fstcw	FPCNTRL
	or	FPCNTRL, 0C00h
	fldcw	FPCNTRL

; Check if we are doing an exact GCD

	cmp	_NUMARG, 2
	jle	simple

; Do as many single precision operations as we can
;	FPU contains U, V, D, B, C, A
; As Knuth suggests:
;	Compute (U-B)/(V+D), the smaller quotient
;	Compute (U+A)/(V-C), the larger quotient, break if not equal
;	Set newB = D, newD = B+Q*D, break if newD won't fit in 32 bits
;	Set newA = C, newC = A+Q*C
;	Set newU = V, newV = U-Q*V

dloop:	fld	st(0)			; U, U, V, D, B, C, A
	fsub	st(0), st(4)		; U-B, U, V, D, B, C, A
	fld	st(2)			; V, U-B, U, V, D, B, C, A
	fadd	st(0), st(4)		; V+D, U-B, U, V, D, B, C, A
	fdivp	st(1), st(0)		; Q, U, V, D, B, C, A
	fadd	BIGVAL
	fld	st(2)			; V, Q, U, V, D, B, C, A
	fsub	st(0), st(6)		; V-C, Q, U, V, D, B, C, A
	fxch	st(1)			; Q, V-C, U, V, D, B, C, A
	fsub	BIGVAL
	fst	QWORD PTR fltval
	fmul	st(0), st(1)		; Q*(V-C), V-C, U, V, D, B, C, A
	faddp	st(1), st(0)		; (Q+1)*(V-C), U, V, D, B, C, A
	fld	QWORD PTR fltval	; Q, (Q+1)*(V-C), U, V, D, B, C, A
	fmul	st(0), st(4)		; Q*D, (Q+1)*(V-C), U, V, D, B, C, A
	fxch	st(1)			; (Q+1)*(V-C), Q*D, U, V, D, B, C, A
	fsub	st(0), st(2)		; (Q+1)*(V-C)-U, Q*D, U, V, D, B, C, A
	fxch	st(1)			; Q*D, (Q+1)*(V-C)-U, U, V, D, B, C, A
	fadd	st(0), st(5)		; B+Q*D,(Q+1)*(V-C)-U, U, V, D, B, C, A
	fxch	st(1)			; (Q+1)*(V-C)-U,B+Q*D, U, V, D, B, C, A
	fcomp	st(7)			; B+Q*D, U, V, D, B, C, A
	fstsw	ax			; Copy comparison results
	and	eax, 4100h		; Isolate C0 & C3 bits
	jnz	lpdone			; Break if less than or equal
	fcom	TWOPOW32		; Will new D be too large?
	fstsw	ax			; Copy comparison results
	and	eax, 0100h		; Isolate C0 bit
	jz	lpdone			; Break if not less than
	fxch	st(4)			; oldB, U, V, newB, newD, C, A
	fcomp	st(1)			; U, V, newB, newD, C, A
	fld	QWORD PTR fltval	; Q, U, V, newB, newD, C, A
	fmul	st, st(5)		; Q*C, U, V, newB, newD, C, A
	fld	QWORD PTR fltval	; Q, Q*C, U, V, newB, newD, C, A
	fmul	st, st(3)		; Q*V, Q*C, U, V, newB, newD, C, A
	fxch	st(1)			; Q*C, Q*V, U, V, newB, newD, C, A
	faddp	st(7), st		; Q*V, U, V, newB, newD, newA, newC
	fsubp	st(1), st		; newV, newU, newB, newD, newA, newC

;	Compute (U-A)/(V+C), the smaller quotient
;	Compute (U+B)/(V-D), the larger quotient, break if not equal
;	Set newB = D, newD = B+Q*D, break if newD won't fit in 32 bits
;	Set newA = C, newC = A+Q*C
;	Set newU = V, newV = U-Q*V

	fld	st(1)			; U, V, U, B, D, A, C
	fsub	st(0), st(5)		; U-A, V, U, B, D, A, C
	fld	st(1)			; V, U-A, V, U, B, D, A, C
	fadd	st(0), st(7)		; V+C, U-A, V, U, B, D, A, C
	fdivp	st(1), st(0)		; Q, V, U, B, D, A, C
	fadd	BIGVAL
	fld	st(1)			; V, Q, V, U, B, D, A, C
	fsub	st(0), st(5)		; V-D, Q, V, U, B, D, A, C
	fxch	st(1)			; Q, V-D, V, U, B, D, A, C
	fsub	BIGVAL
	fst	QWORD PTR fltval
	fmul	st(0), st(1)		; Q*(V-D), V-D, V, U, B, D, A, C
	faddp	st(1), st(0)		; (Q+1)*(V-D), V, U, B, D, A, C
	fld	QWORD PTR fltval	; Q, (Q+1)*(V-D), V, U, B, D, A, C
	fmul	st(0), st(5)		; Q*D, (Q+1)*(V-D), V, U, B, D, A, C
	fxch	st(1)			; (Q+1)*(V-D), Q*D, V, U, B, D, A, C
	fsub	st(0), st(3)		; (Q+1)*(V-D)-U, Q*D, V, U, B, D, A, C
	fxch	st(1)			; Q*D, (Q+1)*(V-D)-U, V, U, B, D, A, C
	fadd	st(0), st(4)		; B+Q*D,(Q+1)*(V-D)-U, V, U, B, D, A, C
	fxch	st(1)			; (Q+1)*(V-D)-U,B+Q*D, V, U, B, D, A, C
	fcomp	st(4)			; B+Q*D, V, U, B, D, A, C
	fstsw	ax			; Copy comparison results
	and	eax, 4100h		; Isolate C0 & C3 bits
	jnz	lpdone1			; Break if less than or equal
	fcom	TWOPOW32		; Will new D be too large?
	fstsw	ax			; Copy comparison results
	and	eax, 0100h		; Isolate C0 bit
	jz	lpdone1			; Break if not less than
	fxch	st(3)			; oldB, V, U, newD, newB, A, C
	fcomp	st(1)			; V, U, newD, newB, A, C
	fld	QWORD PTR fltval	; Q, V, U, newD, newB, A, C
	fmul	st, st(6)		; Q*C, V, U, newD, newB, A, C
	fld	QWORD PTR fltval	; Q, Q*C, V, U, newD, newB, A, C
	fmul	st, st(2)		; Q*V, Q*C, V, U, newD, newB, A, C
	fxch	st(1)			; Q*C, Q*V, V, U, newD, newB, A, C
	faddp	st(6), st		; Q*V, V, U, newD, newB, newC, newA
	fsubp	st(2), st		; newU, newV, newD, newB, newC, newA
	jmp	dloop

; The single precision case:
;	Compute U/V, the quotient
;	Set newB = D, newD = B+Q*D, break if newD won't fit in 32 bits
;	Set newA = C, newC = A+Q*C
;	Set newU = V, newV = U-Q*V

simple:	xor	ecx, ecx		; Clear EGCD_ODD flag
sloop:	fld	st(0)			; U, U, V, D, B, C, A
	fdiv	st(0), st(2)		; Q, U, V, D, B, C, A
	fadd	BIGVAL
	fsub	BIGVAL
	fst	QWORD PTR fltval
	fmul	st(0), st(3)		; Q*D, U, V, D, B, C, A
	fadd	st(0), st(4)		; B+Q*D, U, V, D, B, C, A
	fcom	TWOPOW32		; Will new D be too large?
	fstsw	ax			; Copy comparison results
	and	eax, 0100h		; Isolate C0 bit
	jz	short lpdone2		; Break if not less than
	fxch	st(3)			; newB, U, V, newD, oldB, C, A
	fxch	st(4)			; oldB, U, V, newD, newB, C, A
	fcomp	st(1)			; U, V, newD, newB, C, A
	fld	QWORD PTR fltval	; Q, U, V, newD, newB, C, A
	fmul	st, st(5)		; Q*C, U, V, newD, newB, C, A
	fld	QWORD PTR fltval	; Q, Q*C, U, V, newD, newB, C, A
	fmul	st, st(3)		; Q*V, Q*C, U, V, newD, newB, C, A
	fxch	st(1)			; Q*C, Q*V, U, V, newD, newB, C, A
	faddp	st(7), st		; Q*V, U, V, newD, newB, newA, newC
	fxch	st(5)
	fxch	st(6)
	fxch	st(5)			; Q*V, U, V, newD, newB, newC, newA
	fsubp	st(1), st		; newV, newU, newD, newB, newC, newA
	inc	ecx			; Toggle EGCD_ODD flag
	ftst
	fxch	st(1)			; newU, newV, newD, newB, newC, newA
	fstsw	ax			; Copy comparison results
	and	eax, 4000h		; Isolate C3 bit
	jz	short sloop		; Loop if V is not zero

; Copy extended GCD info to globals
; FPU has U, V, D, B, C, A

	fld	st(0)			; Push value for following code to pop
lpdone2:and	ecx, 1
	mov	_EGCD_ODD, ecx
	jmp	short lpd1

; Copy extended GCD info to globals
; FPU has B+Q*D, U, V, D, B, C, A

lpdone:	mov	_EGCD_ODD, 0
lpd1:	fcompp				; Pop 3 values
	fcomp	st(3)
	fistp	QWORD PTR fltval	; D
	fistp	QWORD PTR fltval1	; B
	mov	eax, fltval
	mov	edx, fltval1
	cmp	edx, 0			; Check for special case
	je	short lpd2
	mov	_EGCD_D, eax
	mov	_EGCD_B, edx
	fistp	QWORD PTR fltval	; C
	fistp	QWORD PTR fltval1	; A
	mov	eax, fltval
	mov	edx, fltval1
	mov	_EGCD_C, eax
	mov	_EGCD_A, edx
	jmp	egdone

; In this case the main loop couldn't determine even a single quotient.
; Return FALSE, so caller can handle by more sophisticated means.

lpd2:	fcompp				; Pop 2 values
	xor	eax, eax
	jmp	short egdone1

; Just like lpdone, but handles the second chunk of dloop code

lpdone1:fcompp				; Pop 3 values
	fcomp	st(3)
	fistp	QWORD PTR fltval	; B
	fistp	QWORD PTR fltval1	; D
	mov	eax, fltval
	mov	edx, fltval1
	mov	_EGCD_B, eax
	mov	_EGCD_D, edx
	fistp	QWORD PTR fltval	; A
	fistp	QWORD PTR fltval1	; C
	mov	eax, fltval
	mov	edx, fltval1
	mov	_EGCD_A, eax
	mov	_EGCD_C, edx
	mov	_EGCD_ODD, 1
	jmp	short egdone

; Restore FPU control word, return TRUE

egdone:	mov	eax, 1
egdone1:fstcw	FPCNTRL
	and	FPCNTRL, 0F3FFh
	fldcw	FPCNTRL
	pop	esi
	pop	edi
	pop	ebx
	pop	ebp


; The x86-64 version uses integer instructions


ELSE
	push	rbp
	push	rbx
	push	rdi
	push	rsi
	push	r12
	push	r14
	push	r15

; Init extended GCD information
;	EGCD_A = 1;			; R8
;	EGCD_B = 0;			; R9
;	EGCD_C = 0;			; R10
;	EGCD_D = 1;			; R11
;	EGCD_ODD = 0;			; R12

	mov	r8, 1			; A
	xor	r9, r9			; B
	xor	r10, r10		; C
	mov	r11, 1			; D
	xor	r12, r12		; ODD

; Load up to 64 bits of U and V
;	U will be in r14
;	V will be in r15

	mov	rdi, _SRCARG		; Giant U
	mov	ebx, _NUMARG		; U->sign
	mov	r14d, [rdi-4][rbx*4]	; U[Ulen-1]

	mov	rsi, _SRC2ARG		; Giant V
	mov	ebp, _NUM2ARG		; V->sign
	xor	r15, r15		; Zero V in case U and V not same len
	cmp	rbx, rbp		; Load top V word if U and V same len
	jne	short noload
	mov	r15d, [rsi-4][rbx*4]	; V[Ulen-1]

	cmp	rbx, 1			; Are there more words to shift
	je	simple			; bits from?

noload:	shl	r14, 32
	shl	r15, 32
	mov	eax, [rdi-8][rbx*4]	; U[Ulen-2]
	mov	edx, [rsi-8][rbx*4]	; V[Ulen-2]
	add	r14, rax
	add	r15, rdx

	cmp	ebx, 2			; Are there more words to shift
	je	simple			; bits from?

	mov	eax, [rdi-12][rbx*4]	; U[Ulen-3]
	mov	edx, [rsi-12][rbx*4]	; V[Ulen-3]

	bsr	rcx, r14		; Count bits to shift U
	xor	cl, 63			; Turn bit # into a shift count

	shl	rax, 32
	shl	rdx, 32
	shld	r14, rax, cl		; Shift U
	shld	r15, rdx, cl		; Shift V

; Do as many operations as we can constrained by a maximum 32-bit result
;	FPU contains U, V, D, B, C, A
; As Knuth suggests:
;	Compute (U-B)/(V+D), the smaller quotient
;	Compute (U+A)/(V-C), the larger quotient, break if not equal
;	Set newB = D, newD = B+Q*D, break if newD won't fit in 32 bits
;	Set newA = C, newC = A+Q*C
;	Set newU = V, newV = U-Q*V

	mov	rsi, 1			; Compute 2^32
	shl	rsi, 32
	shr	r14, 1			; Lose a bit so we can use signed
	shr	r15, 1			; multiply and divide instructions
dloop:	mov	rax, r14		; U
	sub	rax, r9			; U-B
	mov	rcx, r15		; V
	add	rcx, r11		; V+D
	xor	rdx, rdx
	idiv	rcx			; Q = (U-B)/(V+D)

	mov	rcx, r15		; V
	sub	rcx, r10		; V-C
	mov	rbp, rax		; Copy Q
	imul	rbp, rcx		; Q*(V-C)
	add	rbp, rcx		; (Q+1)*(V-C)
	lea	rbx, [r14+r8]		; U+A
	cmp	rbx, rbp		; Compare U+A to (Q+1)*(V-C)
	jae	lpdone			; Break if above or equal

	mov	rbx, r11		; D
	imul	rbx, rax		; Q*D
	add	rbx, r9			; B+Q*D (newD)
	cmp	rbx, rsi		; Will new D be too large?
	jae	lpdone			; Break if newD is 32-bits or more
	mov	r9, r11			; newB = D
	mov	r11, rbx		; newD = B+Q*D

	mov	rcx, r10		; C
	imul	rcx, rax		; Q*C
	add	rcx, r8			; A+Q*C (newC)
	mov	r8, r10			; newA = C
	mov	r10, rcx		; newC = A+Q*C

	imul	rax, r15		; Q*V
	sub	r14, rax		; U-Q*V (newV)
	xchg	r14, r15		; newU = V, newV = U-Q*V

	xor	r12, 1			; Flip ODD

;	Compute (U-A)/(V+C), the smaller quotient
;	Compute (U+B)/(V-D), the larger quotient, break if not equal
;	Set newB = D, newD = B+Q*D, break if newD won't fit in 32 bits
;	Set newA = C, newC = A+Q*C
;	Set newU = V, newV = U-Q*V

	mov	rax, r14		; U
	sub	rax, r8			; U-A
	mov	rcx, r15		; V
	add	rcx, r10		; V+C
	xor	rdx, rdx
	idiv	rcx			; Q = (U-A)/(V+C)

	mov	rcx, r15		; V
	sub	rcx, r11		; V-D
	mov	rbp, rax		; Copy Q
	imul	rbp, rcx		; Q*(V-D)
	add	rbp, rcx		; (Q+1)*(V-D)
	lea	rbx, [r14+r9]		; U+B
	cmp	rbx, rbp		; Compare U+B to (Q+1)*(V-D)
	jae	lpdone			; Break if above or equal

	mov	rbx, r11		; D
	imul	rbx, rax		; Q*D
	add	rbx, r9			; B+Q*D (newD)
	cmp	rbx, rsi		; Will new D be too large?
	jae	lpdone			; Break if newD is 32-bits or more
	mov	r9, r11			; newB = D
	mov	r11, rbx		; newD = B+Q*D

	mov	rcx, r10		; C
	imul	rcx, rax		; Q*C
	add	rcx, r8			; A+Q*C (newC)
	mov	r8, r10			; newA = C
	mov	r10, rcx		; newC = A+Q*C

	imul	rax, r15		; Q*V
	sub	r14, rax		; U-Q*V (newV)
	xchg	r14, r15		; newU = V, newV = U-Q*V

	xor	r12, 1			; Flip ODD
	jmp	dloop

; The single precision case:
;	Compute U/V, the quotient
;	Set newB = D, newD = B+Q*D, break if newD won't fit in 32 bits
;	Set newA = C, newC = A+Q*C
;	Set newU = V, newV = U-Q*V

simple:	mov	rsi, 1			; Compute 2^32
	shl	rsi, 32
sloop:	mov	rax, r14		; U
	xor	rdx, rdx
	div	r15			; Q = U/V
	mov	rbp, rax		; Copy Q

	mul	r11			; Q*D
	add	rax, r9			; B+Q*D
	cmp	rax, rsi		; Will new D be too large?
	jae	short lpdone		; Break if D >= 2^32
	mov	r9, r11			; newB = D
	mov	r11, rax		; newD = B+Q*D

	mov	rax, r10		; C
	mul	rbp			; Q*C
	add	rax, r8			; A+Q*C (newC)
	mov	r8, r10			; newA = C
	mov	r10, rax		; newC = A+Q*C

	mov	rax, r15		; V
	mul	rbp			; Q*V
	sub	r14, rax		; U-Q*V (newV)
	xchg	r14, r15		; newU = V, newV = U-Q*V

	xor	r12, 1			; Flip ODD

	and	r15, r15		; Loop if V is not zero
	jnz	short sloop

; Copy extended GCD info to globals

lpdone:	mov	_EGCD_A, r8d
	mov	_EGCD_B, r9d
	mov	_EGCD_C, r10d
	mov	_EGCD_D, r11d
	mov	_EGCD_ODD, r12d
	and	r9, r9			; Check for special case (B = 0)
	jnz	short egdone

; In this case the main loop couldn't determine even a single quotient.
; Return FALSE, so caller can handle by more sophisticated means.

	xor	rax, rax
	jmp	short egdone1

egdone:	mov	rax, 1
egdone1:pop	r15
	pop	r14
	pop	r12
	pop	rsi
	pop	rdi
	pop	rbx
	pop	rbp
ENDIF
	ret
_egcdhlp ENDP

_TEXT	ENDS
END
