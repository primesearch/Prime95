; Copyright 1995-2003 Just For Fun Software, Inc., all rights reserved
; Author:  George Woltman
; Email: woltman@alum.mit.edu
;

INCLUDE pfa.mac
INCLUDE lucas.mac
INCLUDE lucasp.mac
INCLUDE xmult.mac
INCLUDE xlucas.mac
INCLUDE xnormal.mac
INCLUDE xpass2.mac

X87_CASES	EQU	13
SSE2_CASES	EQU	22

loopent	MACRO	y,z		; Create a entry in the loop entry table
	DD	&y&z
	ENDM
looptab	MACRO	y, cnt		; Create the loop entry table
	x = 0
	REPT	cnt
	loopent	y, %x
	x = x + 1
	ENDM
	ENDM

x87table: looptab case, X87_CASES
sse2table: looptab sscase, SSE2_CASES

;; Macros from p4notes.doc

read1	MACRO	mem, c			;; Bytes to read
	LOCAL loop1, loop2
	cnt = mem/64		        ;; Read 64 bytes per iteration
	mov	edx, c
loop2:	mov     ecx, cnt
loop1:	mov     eax, [esi]		;; Read one cache line
	mov     eax, [esi+4]		;; 4 bytes at a time
	mov     eax, [esi+8]
	mov     eax, [esi+12]
	mov     eax, [esi+16]
	mov     eax, [esi+20]
	mov     eax, [esi+24]
	mov     eax, [esi+28]
	mov     eax, [esi+32]
	mov     eax, [esi+36]
	mov     eax, [esi+40]
	mov     eax, [esi+44]
	mov     eax, [esi+48]
	mov     eax, [esi+52]
	mov     eax, [esi+56]
	mov     eax, [esi+60]
	lea     esi, [esi+64]			; Next cache line
	sub     ecx, 1
	jnz     loop1
	lea     esi, [esi-mem]		; Restore esi
	dec	edx
	jnz	loop2
	ENDM

read2	MACRO	mem, c			;; Bytes to read
	LOCAL loop1, loop2
        cnt = mem/64                    ; 64 bytes per iteration
	mov	edx, c
loop2:	mov     ecx, cnt
loop1:  movapd  xmm1, [esi]			; Read one cache line
        movapd  xmm1, [esi+16]		; 16 bytes at a time
        movapd  xmm1, [esi+32]
        movapd  xmm1, [esi+48]
        lea     esi, [esi+64]			; Next cache line
        sub     ecx, 1
        jnz     loop1
        lea     esi, [esi-mem]
	dec	edx
	jnz	loop2
	ENDM

write1	MACRO	mem, c			;; Bytes to write
	LOCAL loop1, loop2
        cnt = mem/64                    ; 64 bytes per iteration
	mov	edx, c
loop2:	mov     ecx, cnt
loop1:  movapd  [esi], xmm1			; Write one cache line
        movapd  [esi+16], xmm1		; 16 bytes at a time
        movapd  [esi+32], xmm1
        movapd  [esi+48], xmm1
        lea     esi, [esi+64]			; Next cache line
        sub     ecx, 1
        jnz     loop1
        lea     esi, [esi-mem]
	dec	edx
	jnz	loop2
	ENDM

write2	MACRO	mem, c			;; Bytes to write
	LOCAL loop1, loop2
        cnt = mem/(4*128)               ; 128 bytes per iteration
        dist = 64
	mov	edx, c
	sub     ebx, ebx
loop2:	mov     ecx, cnt
loop1:  movapd  [esi+0*dist], xmm1      ; Write 8 cache lines
        movapd  [esi+1*dist], xmm1
        movapd  [esi+2*dist], xmm1
        movapd  [esi+3*dist], xmm1
        movapd  [esi+4*dist], xmm1
        movapd  [esi+5*dist], xmm1
        movapd  [esi+6*dist], xmm1
        movapd  [esi+7*dist], xmm1
        lea     esi, [esi+16]           ; Same cache lines
        add     bl, 256/4			; 4 inner loop iterations
        jnc     loop1
        lea     esi, [esi-4*16+8*dist]  ; Next set of 8 cache lines
        sub     ecx, 1
        jnz     loop1
        lea     esi, [esi-mem]
	dec	edx
	jnz	loop2
	ENDM

readwrite1 MACRO mem, c			;; Bytes to write
	LOCAL loop1, loop2
        cnt = mem/64
	mov	edx, c
loop2:	mov     ecx, cnt
loop1:  movapd  xmm0, [esi]             ; Read one cache line
        movapd  xmm1, [esi+16]
        movapd  xmm2, [esi+32]
        movapd  xmm3, [esi+48]
        subpd   xmm0, xmm0              ; Operate on the data
        pxor    xmm1, xmm1
        subpd   xmm2, xmm2
        pxor    xmm3, xmm3
        movapd  [esi], xmm0             ; Write the cache line
        movapd  [esi+16], xmm1
        movapd  [esi+32], xmm2
        movapd  [esi+48], xmm3
        lea     esi, [esi+64]           ; Next cache line
        sub     ecx, 1
        jnz     loop1
        lea     esi, [esi-mem]
	dec	edx
	jnz	loop2
	ENDM

;; Time one of the basic FFT building blocks

x87mac	MACRO	memused, memarea, ops:vararg
	LOCAL	ss0a, ss0b
	inner_iters = memarea / memused
	outer_iters = 10000 / inner_iters
	mov	eax, outer_iters
ss0a:	mov	edi, _SRCARG
	lea	esi, [edi+4096]
	mov	ecx, inner_iters
ss0b:	disp	&ops
	lea	esi, [esi+memused]	;; Next source pointer
	lea	edi, [edi+SCD]		;; Next sine/cosine pointer
	dec	ecx
	jnz	ss0b
	dec	eax
	jnz	ss0a
	ENDM

sse2mac	MACRO	memused, memarea, ops:vararg
	LOCAL	ss0a, ss0b
	inner_iters = memarea / memused
	outer_iters = 10000 / inner_iters
	mov	eax, outer_iters
ss0a:	mov	edi, _SRCARG
	lea	esi, [edi+4096]
	mov	ecx, inner_iters
ss0b:	&ops
	lea	edi, [edi+2*XMM_SCD]	;; Next sine/cosine pointer
	dec	ecx
	jnz	ss0b
	dec	eax
	jnz	ss0a
	ENDM


PUBLIC _timeit
_timeit	PROC NEAR
	push	esi
	push	edi
	push	ebp
	push	ebx

	clear_timers

	sub	ebx, ebx
	sub	ecx, ecx
	sub	ebp, ebp
	mov	edi, _SRCARG
	mov	esi, _SRCARG
	
	mov	edx, _DESTARG		; Get test number

	mov	eax, X87_CASES
	cmp	edx, -1			; -1 = get num x87 cases
	je	exit
	mov	eax, SSE2_CASES
	cmp	edx, -2			; -2 = get num sse2 cases
	je	exit

	cmp	edx, 1000		; Tests above 1000 are SSE2 code
	jl	short x87
	subpd	xmm0, xmm0
	subpd	xmm1, xmm1
	subpd	xmm2, xmm2
	subpd	xmm3, xmm3
	subpd	xmm4, xmm4
	subpd	xmm5, xmm5
	subpd	xmm6, xmm6
	subpd	xmm7, xmm7

; Jump to desired test case

	sub	edx, 1000
	mov	eax, sse2table[edx*4]	; Get address of test to execute
	jmp	eax
x87:	mov	eax, x87table[edx*4]	; Get address of test to execute
	jmp	eax

; Time the "do-nothing" case.

case0:
	jmp	exit

; Time the loop of "do-nothing" case.  1000 iterations.

case1:
	mov	ecx, 1000
c1a:	dec	ecx
	jnz	c1a
	jmp	exit

; This code reads a contiguous block of memory.
; Timings are done on 3 memory size.  4KB will read from the L1 cache
; only, 96KB will read from the L2 cache only, and 2MB will test reading
; from main memory.

case2:	read1	4096, 1000	; Read 4KB
	jmp	exit

case3:	read1	96*1024, 100	; Read 96KB
	jmp	exit

case4:	read1	2048*1024, 2	; Read 2MB
	jmp	exit

case5:	x87mac	64, 4096, eight_reals_fft, 8, 16, 32
	jmp	exit

case6:	x87mac	64, 100000, eight_reals_fft, 8, 16, 32
	jmp	exit

case7:	x87mac	64, 4096, eight_reals_unfft, 8, 16, 32
	jmp	exit

case8:	x87mac	64, 100000, eight_reals_unfft, 8, 16, 32
	jmp	exit

case9:	x87mac	64, 4096, four_complex_fft, 8, 16, 32
	jmp	exit

case10:	x87mac	64, 100000, four_complex_fft, 8, 16, 32
	jmp	exit

case11:	x87mac	64, 4096, four_complex_unfft, 8, 16, 32
	jmp	exit

case12:	x87mac	64, 100000, four_complex_unfft, 8, 16, 32
	jmp	exit

; This code reads a contiguous block of memory.
; Timings are done on 3 memory size.  4KB will read from the L1 cache
; only, 96KB will read from the L2 cache only, and 2MB will test reading
; from main memory.

sscase0:
	read2	4096, 1000	; Read 4KB
	jmp	exit

sscase1:
	read2	96*1024, 100	; Read 96KB
	jmp	exit

sscase2:
	read2	2048*1024, 2	; Read 2MB
	jmp	exit

; This code writes a contiguous block of memory.
; Timings are done on 3 memory size.  4KB will write to the L1 cache
; only, 96KB will write to L2 cache only, and 2MB will test writing
; to main memory.

sscase3:
	write1	4096, 1000	; Read 4KB
	jmp	exit

sscase4:
	write1	96*1024, 100	; Read 96KB
	jmp	exit

sscase5:
	write1	2048*1024, 2	; Read 2MB
	jmp	exit

; This code writes a block of memory non-contiguously.
; Timings are done on 3 memory size.  4KB will write to the L1 cache
; only, 96KB will write to L2 cache only, and 2MB will test writing
; to main memory.

sscase6:
	write2	4096, 1000	; Read 4KB
	jmp	exit

sscase7:
	write2	96*1024, 100	; Read 96KB
	jmp	exit

sscase8:
	write2	2048*1024, 2	; Read 2MB
	jmp	exit

; This code reads & writes a block of memory.
; Timings are done on 3 memory size.  4KB will write to the L1 cache
; only, 96KB will write to L2 cache only, and 2MB will test writing
; to main memory.

sscase9:
	readwrite1	4096, 1000	; Read 4KB
	jmp	exit

sscase10:
	readwrite1	96*1024, 100	; Read 96KB
	jmp	exit

sscase11:
	readwrite1	2048*1024, 2	; Read 2MB
	jmp	exit

; Time ~10000 iterations of the SSE2 macros in L1 and L2 caches

sscase12:
	sse2mac	128, 4096, x2cl_eight_reals_fft esi, 2*64, 64
	jmp	exit
sscase13:
	sse2mac	128, 100000, x2cl_eight_reals_fft esi, 2*64, 64
	jmp	exit
sscase14:
	sse2mac	256, 4096, x4cl_eight_reals_last_unfft esi, 4*64, 64, 2*64
	jmp	exit
sscase15:
	sse2mac	256, 100000, x4cl_eight_reals_last_unfft esi, 4*64, 64, 2*64
	jmp	exit

sscase16:
	sse2mac	256, 4096, x4cl_four_complex_fft esi, 4*64, 64, 2*64
	jmp	exit
sscase17:
	sse2mac	256, 100000, x4cl_four_complex_fft esi, 4*64, 64, 2*64
	jmp	exit
sscase18:
	sse2mac	256, 4096, x4cl_four_complex_unfft esi, 4*64, 64, 2*64
	jmp	exit
sscase19:
	sse2mac	256, 100000, x4cl_four_complex_unfft esi, 4*64, 64, 2*64
	jmp	exit

sscase20:
	sse2mac	128, 4096, x2cl_four_complex_fft esi, 2*64, 64
	jmp	exit
sscase21:
	sse2mac	128, 100000, x2cl_four_complex_fft esi, 2*64, 64
	jmp	exit

; Exit the timing code

exit:	pop	ebx
	pop	ebp
	pop	edi
	pop	esi
	ret

_timeit	ENDP
