; Copyright 1995-2002 Just For Fun Software, Inc., all rights reserved
; Author:  George Woltman
; Email: woltman@alum.mit.edu
;
; This file implements helper routines for the CPU identification code
; as well as other helper functions
;

	TITLE   setup

	.586

_TEXT32 SEGMENT PAGE USE32 PUBLIC 'DATA'

EXTRN	_CPUID_EAX:DWORD
EXTRN	_CPUID_EBX:DWORD
EXTRN	_CPUID_ECX:DWORD
EXTRN	_CPUID_EDX:DWORD

;
; Global variables
;

_TEXT32	ENDS

	ASSUME  CS: _TEXT32, DS: _TEXT32, SS: _TEXT32, ES: _TEXT32
_TEXT32	SEGMENT

;
; Utility routine to initialize the FPU
;

	PUBLIC	_fpu_init
_fpu_init PROC NEAR
	fninit
	ret
_fpu_init ENDP

;
; Utility routine to read the time stamp counter
;

	PUBLIC	_erdtsc
_erdtsc PROC NEAR
	rdtsc
	mov	_CPUID_EAX, eax		; low 32 bits
	mov	_CPUID_EDX, edx		; high 32 bits
	ret
_erdtsc ENDP


;
; Utility routine to see if CPUID is supported
; Returns non-zero if CPUID is supported
;

	PUBLIC	_ecpuidsupport
_ecpuidsupport PROC NEAR
        pushfd				; Get original EFLAGS
	pop	eax
	mov 	ecx, eax
        xor     eax, 200000h		; Flip ID bit in EFLAGS
        push    eax			; Save new EFLAGS value on stack
        popfd				; Replace current EFLAGS value
        pushfd				; Get new EFLAGS
        pop     eax			; Store new EFLAGS in EAX
        xor     eax, ecx		; Test toggled ID bit
	ret				; Processor supports CPUID if eax != 0
_ecpuidsupport ENDP

;
; Utility routine to execute CPUID
;

	PUBLIC	_ecpuid
_ecpuid PROC NEAR
	push	ebx
	mov	eax, _CPUID_EAX		; Argument to CPUID
	cpuid				; Perform the CPUID
	mov	_CPUID_EAX, eax		; Return the results
	mov	_CPUID_EBX, ebx		; Return the results
	mov	_CPUID_ECX, ecx		; Return the results
	mov	_CPUID_EDX, edx		; Return the results
	pop	ebx
	ret
_ecpuid ENDP


_TEXT32	ENDS
END
