; Copyright 1995-2000 Just For Fun Software, Inc., all rights reserved
; Author:  George Woltman
; Email: woltman@alum.mit.edu
;
; This file implements helper routines for the ECM code
;

	TITLE   setup

	.386

_TEXT32 SEGMENT PAGE USE32 PUBLIC 'DATA'

EXTRN	_PARG:DWORD
EXTRN	_FFTLEN:DWORD
EXTRN	_SRCARG:DWORD
EXTRN	_SRC2ARG:DWORD
EXTRN	_DESTARG:DWORD
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

	align	32
fltval		DD	0
fltvalhi	DD	0
fltval1		DD	0
fltval1hi	DD	0
TWOPOW32	DD	4294967296.0		; 2^32
BIGVAL		DD	5F000000h		; 2^63
FPCNTRL		DW	0
_TEXT32	ENDS

	ASSUME  CS: _TEXT32, DS: _TEXT32, SS: _TEXT32, ES: _TEXT32
_TEXT32	SEGMENT

INCLUDE unravel.mac

; Add integer to carryh,carryl,res

	PUBLIC	_eaddhlp
_eaddhlp PROC NEAR
	mov	eax, _SRCARG		; Integer1
	add	_RES, eax		; Add two integers
	adc	_CARRYL, 0		; Add the carry
	adc	_CARRYH, 0		; Add the carry
	ret
_eaddhlp ENDP

; Subtract integer from carryh,carryl,res

	PUBLIC	_esubhlp
_esubhlp PROC NEAR
	mov	eax, _SRCARG		; Integer1
	sub	_RES, eax		; Subtract two integers
	sbb	_CARRYL, 0		; Subtract the carry
	sbb	_CARRYH, 0		; Subtract the carry
	ret
_esubhlp ENDP

; Multiply integer1 and integer2 adding result to carryh,carryl,res

	PUBLIC	_emuladdhlp
_emuladdhlp PROC NEAR
	mov	eax, _SRCARG		; Integer1
	mul	_SRC2ARG		; Integer2
	add	_RES, eax		; Add result to 3 word accumulator
	adc	_CARRYL, edx
	adc	_CARRYH, 0
	ret
_emuladdhlp ENDP

; Multiply integer1 and integer2 adding twice the result to carryh,carryl,res

	PUBLIC	_emuladd2hlp
_emuladd2hlp PROC NEAR
	mov	eax, _SRCARG		; Integer1
	mul	_SRC2ARG		; Integer2
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
_emulsubhlp PROC NEAR
	mov	eax, _SRCARG		; Integer1
	mul	_SRC2ARG		; Integer2
	sub	_RES, eax		; Sub result from 3 word accumulator
	sbb	_CARRYL, edx
	sbb	_CARRYH, 0
	ret
_emulsubhlp ENDP

; Routine to help in computing extended GCD quickly
;
; Do several single-precision steps for the extended GCD code
; U is larger than V - both are the same length.  Returns A,B,C,D
; as defined in Knuth vol 2. description of extended GCD for
; large numbers.
; This was implemented in assembly language to deal with 64-bit integers.

	PUBLIC	_egcdhlp
_egcdhlp PROC NEAR
	fninit
	push	ebp
	push	ebx
	push	edi
	push	esi

; Load up to 64 bits of U and V
;	U will be in edx:eax with bits shifted in from edi
;	V will be in esi:ebp with bits shifted in from ebx

	mov	edi, _SRCARG		; Giant U
	mov	ebx, [edi]		; U->sign
	mov	edi, [edi+4]		; U->n
	mov	ecx, _SRC2ARG		; Giant V
	mov	ebp, [ecx]		; V->sign
	mov	ecx, [ecx+4]		; V->n

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

	mov	ebx, _SRCARG		; Giant U
	mov	ebx, [ebx]		; U->sign
	cmp	ebx, 2
	JLE_X	simple

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
	JNZ_X	lpdone			; Break if less than or equal
	fcom	TWOPOW32		; Will new D be too large?
	fstsw	ax			; Copy comparison results
	and	eax, 0100h		; Isolate C0 bit
	JZ_X	lpdone			; Break if not less than
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
	JNZ_X	lpdone1			; Break if less than or equal
	fcom	TWOPOW32		; Will new D be too large?
	fstsw	ax			; Copy comparison results
	and	eax, 0100h		; Isolate C0 bit
	JZ_X	lpdone1			; Break if not less than
	fxch	st(3)			; oldB, V, U, newD, newB, A, C
	fcomp	st(1)			; V, U, newD, newB, A, C
	fld	QWORD PTR fltval	; Q, V, U, newD, newB, A, C
	fmul	st, st(6)		; Q*C, V, U, newD, newB, A, C
	fld	QWORD PTR fltval	; Q, Q*C, V, U, newD, newB, A, C
	fmul	st, st(2)		; Q*V, Q*C, V, U, newD, newB, A, C
	fxch	st(1)			; Q*C, Q*V, V, U, newD, newB, A, C
	faddp	st(6), st		; Q*V, V, U, newD, newB, newC, newA
	fsubp	st(2), st		; newU, newV, newD, newB, newC, newA
	JMP_X	dloop

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
	JMP_X	egdone

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
	ret
_egcdhlp ENDP


_TEXT32	ENDS
END
