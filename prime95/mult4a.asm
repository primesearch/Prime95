; Copyright 1995-2003 Just For Fun Software, Inc., all rights reserved
; Author:  George Woltman
; Email: woltman@alum.mit.edu
;
; This code handles the middle 6 or 7 levels of three-pass FFTs
; that use the very convoluted memory model.
;

	TITLE   setup

	.686
	.XMM

_TEXT32 SEGMENT PARA USE32 PUBLIC 'DATA'

	ASSUME  CS: _TEXT32, DS: _TEXT32, SS: _TEXT32, ES: _TEXT32

INCLUDE extrn.mac
INCLUDE	unravel.mac
INCLUDE	lucas.mac
INCLUDE mult.mac
INCLUDE pass1.mac
INCLUDE pass2.mac
INCLUDE memory.mac

IFDEF PPRO
INCLUDE	lucasp.mac
ENDIF

	very_convoluted_distances

	PUBLICP	pass1_6_levels_fft
	PUBLICP	pass1_6_levels_unfft
	PUBLICP	pass1_6a_levels_fft
	PUBLICP	pass1_6a_levels_unfft
	PUBLICP	pass1_6_levels_fftp
	PUBLICP	pass1_6_levels_unfftp
	PUBLICP	pass1_7_levels_fft
	PUBLICP	pass1_7_levels_unfft
	PUBLICP	pass1_7a_levels_fft
	PUBLICP	pass1_7a_levels_unfft
	PUBLICP	pass1_7_levels_fftp
	PUBLICP	pass1_7_levels_unfftp

PROCP	pass1_procs4

;; Do 6 levels of pass 1 of the forward FFT
;; Caller must set ecx for the proper number of complex iterations

LABELP	pass1_6a_levels_fft
	pass1_six_levels_real64_fft	;; The all real sub-section
	JMP_X	b7c			;; V - Jump to complex sections

LABELP	pass1_6_levels_fft
	pass1_six_levels_real128_fft	;; The all real sub-sections

LABELP	pass1_6_levels_fftp
b7c:	mov	edx, pass1_premults	;; V - Address of the group multipliers
b8b:	pass1_six_levels_complex_fft	;; Do an all complex sub-section
	lea	edx, [edx+7*PMD]	;; V - Next group pre-multipliers
	dec	cl			;; U - Test if we are done
	jz	short b9b		;; V - Jump if done
	add	ch, 256/4		;; U - Test loop counter
	JNC_X	b8b			;; V - Iterate if necessary
	lea	esi, [esi-8*dist16K+dist128K];; U - Next source pointer
	dec	ch			;; V - Test loop counter
	JNZ_X	b8b			;;*V - Iterate if necessary
	lea	esi, [esi-8*dist128K+dist1M];; U - Next source pointer
	mov	ch, 8			;; V - Restore loop counter
	JMP_X	b8b			;;*V - Iterate
b9b:	ret

;; Do 6 levels of pass 1 of the inverse FFT
;; Caller must set ecx for the proper number of complex iterations

LABELP	pass1_6_levels_unfftp
	IFDEF	PFETCH
	mov	norm_ptr1, esi
	ENDIF
	JMP_X	c7c
LABELP	pass1_6a_levels_unfft
	IFDEF	PFETCH
	lea	ebx, [esi+2*dist16K]
	mov	norm_ptr1, esi
	ENDIF
c7a:	pass1_six_levels_real64_unfft	;; The all real sub-section
	add	cl, 256/2		;; U - Test loop counter
	JNC_X	c7a			;; V - Iterate if necessary
	lea	esi, [esi-2*dist128+dist16K];; U - Next source pointer
	JMP_X	c7c			;; V - Jump to complex sections
LABELP	pass1_6_levels_unfft
	IFDEF	PFETCH
	lea	ebx, [esi+2*dist16K]
	mov	norm_ptr1, esi
	ENDIF
c7b:	pass1_six_levels_real128_unfft	;; The all real sub-section
	add	cl, 256/2		;; U - Test loop counter
	JNC_X	c7b			;; V - Iterate if necessary
	lea	esi, [esi-2*dist128+2*dist16K];; U - Next source pointer
c7c:	mov	edx, pass1_premults	;; V - Address of the group multipliers
c8b:	pass1_six_levels_complex_unfft	;; Do an all complex sub-section
	add	cl, 256/2		;; U - Test inner loop counter
	JNC_X	c8b			;; V - Iterate if necessary
	lea	esi, [esi-2*dist128+2*dist16K];; U - Next source pointer
	lea	edx, [edx+7*PMD]	;; V - Next group pre-multipliers
	dec	cl			;; U - Test if we are done
	jz	short c9b		;; V - Jump if done
	add	ch, 256/4		;; U - Test loop counter
	JNC_X	c8b			;; V - Iterate if necessary
	lea	esi, [esi-8*dist16K+dist128K];; U - Next source pointer
	dec	ch			;; V - Test loop counter
	JNZ_X	c8b			;;*V - Iterate if necessary
	lea	esi, [esi-8*dist128K+dist1M];; U - Next source pointer
	mov	ch, 8			;; V - Restore loop counter
	JMP_X	c8b			;;*V - Iterate
c9b:	ret

;; Do 7 levels of pass 1 of the forward FFT
;; Caller must set ecx for the proper number of complex iterations

LABELP	pass1_7_levels_fft
	pass1_seven_levels_real256_fft	;; The all real sub-section
LABELP	pass1_7_levels_fftp
	mov	edx, pass1_premults	;; V - Address of the group multipliers
b2b:	pass1_seven_levels_complex_fft	;; Do an all complex sub-section
	lea	edx, [edx+11*PMD]	;; V - Next group pre-multipliers
	dec	cl			;; U - Test if we are done
	jz	short b3b		;; V - Jump if done
	add	ch, 256/2		;; U - Test loop counter
	JNC_X	b2b			;; V - Iterate if necessary
	lea	esi, [esi-8*dist16K+dist128K];; U - Next source pointer
	dec	ch			;; V - Test loop counter
	JNZ_X	b2b			;;*V - Iterate if necessary
	lea	esi, [esi-8*dist128K+dist1M];; U - Next source pointer
	mov	ch, 8			;; V - Restore loop counter
	JMP_X	b2b			;;*V - Iterate
b3b:	ret

;; Do 7 levels of pass 1 of the inverse FFT
;; Caller must set ecx for the proper number of complex iterations

LABELP	pass1_7_levels_unfftp
	IFDEF	PFETCH
	mov	norm_ptr1, esi
	ENDIF
	JMP_X	c33b
LABELP	pass1_7_levels_unfft
	IFDEF	PFETCH
	lea	ebx, [esi+4*dist16K]
	mov	norm_ptr1, esi
	ENDIF
c1b:	pass1_seven_levels_real256_unfft;; The all real sub-section
	add	cl, 256/2		;; U - Test loop counter
	JNC_X	c1b			;; V - Iterate if necessary
	lea	esi, [esi-2*dist128+4*dist16K];; U - Next source pointer
c33b:	mov	edx, pass1_premults	;; V - Address of the group multipliers
c2b:	pass1_seven_levels_complex_unfft;; Do an all complex sub-section
	add	cl, 256/2		;; U - Test inner loop counter
	JNC_X	c2b			;; V - Iterate if necessary
	lea	esi, [esi-2*dist128+4*dist16K];; U - Next source pointer
	lea	edx, [edx+11*PMD]	;; V - Next group pre-multipliers
	dec	cl			;; U - Test if we are done
	jz	short c3b		;; V - Jump if done
	add	ch, 256/2		;; U - Test loop counter
	JNC_X	c2b			;; V - Iterate if necessary
	lea	esi, [esi-8*dist16K+dist128K];; U - Next source pointer
	dec	ch			;; V - Test loop counter
	JNZ_X	c2b			;;*V - Iterate if necessary
	lea	esi, [esi-8*dist128K+dist1M];; U - Next source pointer
	mov	ch, 8			;; V - Restore loop counter
	JMP_X	c2b			;;*V - Iterate
c3b:	ret

;; Do 7 levels of pass 1 of the forward FFT
;; Caller must set ecx for the proper number of complex iterations
;; This oddball version is for the case where we have just done a
;; five_reals_fft or seven_reals_fft.

LABELP	pass1_7a_levels_fft
	pass1_seven_levels_real128_fft	;; The all real sub-section
	mov	edx, pass1_premults	;; V - Address of the group multipliers
	push	ebx			;; U - Save ebx for caller
	mov	ebx, 2*dist16K		;; V - Put dist32K in ebx
d2b:	pass1_seven_levels_a_complex_fft;; Do an all complex sub-section
	add	cl, 256/2		;; U - Test inner loop counter
	JNC_X	d2b			;; V - Iterate if necessary
	lea	esi, [esi-2*dist128][ebx+2*dist16K];; U - Next source pointer
	lea	edx, [edx+11*PMD]	;; V - Next group pre-multipliers
	dec	cl			;; U - Test if we are done
	jz	short d3b		;; V - Jump if done
	xor	ebx, (2*dist16K XOR (dist128K-6*dist16K))
					;; U - Dist32K = (dist128K-6*dist16K)
	JMP_X	d2b			;; V - Iterate
d3b:	pop	ebx			;; U - Restore ebx for caller
	ret

;; Do 7 levels of pass 1 of the inverse FFT
;; Caller must set ecx for the proper number of complex iterations
;; This oddball version is for the case where we have just done a
;; five_reals_fft or seven_reals_fft.

LABELP	pass1_7a_levels_unfft
e1b:	pass1_seven_levels_real128_unfft;; The all real sub-section
	add	cl, 256/2		;; U - Test loop counter
	JNC_X	e1b			;; V - Iterate if necessary
	lea	esi, [esi-2*dist128+2*dist16K];; U - Next source pointer
	mov	edx, pass1_premults	;; V - Address of the group multipliers
	push	ebx			;; U - Save ebx for caller
	mov	ebx, 2*dist16K		;; V - Put dist32K in ebx
e2b:	pass1_seven_levels_a_complex_unfft;; Do an all complex sub-section
	add	cl, 256/2		;; U - Test inner loop counter
	JNC_X	e2b			;; V - Iterate if necessary
	lea	esi, [esi-2*dist128][ebx+2*dist16K];; U - Next source pointer
	lea	edx, [edx+11*PMD]	;; V - Next group pre-multipliers
	dec	cl			;; U - Test if we are done
	jz	short e3b		;; V - Jump if done
	xor	ebx, (2*dist16K XOR (dist128K-6*dist16K))
					;; U - Dist32K = (dist128K-6*dist16K)
	JMP_X	e2b			;; V - Iterate
e3b:	pop	ebx			;; U - Restore ebx for caller
	ret

ENDPP	pass1_procs4

_TEXT32	ENDS
END
