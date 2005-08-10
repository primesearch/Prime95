/*----------------------------------------------------------------------
| Copyright 1995-2005 Just For Fun Software, Inc., all rights reserved
| Author:  George Woltman
| Email: woltman@alum.mit.edu
|
| This file contains routines to determine the CPU type and speed.
| Plus, as a bonus, you get 3 routines to portably access the high
| resolution timer!
+---------------------------------------------------------------------*/

#ifndef _CPUID_H
#define _CPUID_H

/* This is a C library.  If used in a C++ program, don't let the C++ */
/* compiler mangle names. */

#ifdef __cplusplus
extern "C" {
#endif

/* Include common definitions */

#include "common.h"

/* Routines that the user should call */

void guessCpuType (void);
void guessCpuSpeed (void);

/* Routines to access the high resolution timer */

int isHighResTimerAvailable (void);
double getHighResTimer (void);
double getHighResTimerFrequency (void);

/* Handle the difference between the naming conventions in */
/* C compilers.  We need to to this for global variables that are */
/* referenced by the assembly routines.  Most non-Windows systems */
/* should #define ADD_UNDERSCORES before including this file. */

#ifdef ADD_UNDERSCORES
#define CPU_FLAGS	_CPU_FLAGS
#define CPU_TYPE	_CPU_TYPE
#define CPUID_EAX	_CPUID_EAX
#define CPUID_EBX	_CPUID_EBX
#define CPUID_ECX	_CPUID_ECX
#define CPUID_EDX	_CPUID_EDX
#define fpu_init	_fpu_init
#define erdtsc		_erdtsc
#define ecpuidsupport	_ecpuidsupport
#define ecpuid		_ecpuid
#endif

/* Global variables describing the CPU we are running on */

extern char CPU_BRAND[49];		/* Text description of CPU type */
extern double CPU_SPEED;		/* Actual CPU Speed in MHz */
#define CPU_RDTSC	0x0001
#define CPU_CMOV	0x0002
#define CPU_PREFETCH	0x0004
#define CPU_SSE		0x0008
#define CPU_SSE2	0x0010
#define CPU_MMX		0x0020
#define CPU_3DNOW	0x0040
extern unsigned int CPU_FLAGS;		/* Cpu capabilities */
extern int CPU_L1_CACHE_SIZE;
extern int CPU_L2_CACHE_SIZE;
extern int CPU_L1_CACHE_LINE_SIZE;
extern int CPU_L2_CACHE_LINE_SIZE;
extern int CPU_L1_DATA_TLBS;
extern int CPU_L2_DATA_TLBS;
extern int CPU_L1_SET_ASSOCIATIVE;
extern int CPU_L2_SET_ASSOCIATIVE;

/* Internal global variables */

extern unsigned long CPUID_EAX;	/* For communicating to asm routines */
extern unsigned long CPUID_EBX;
extern unsigned long CPUID_ECX;
extern unsigned long CPUID_EDX;

/* Internal assembly routines */

void fpu_init (void);
unsigned long ecpuidsupport (void);
void ecpuid (void);
#define isCpuidSupported()	ecpuidsupport ()
#define Cpuid(a)		{ CPUID_EAX = a; ecpuid (); }

/* Routine used to time code chunks */
#define rdtsc(hi,lo)	erdtsc(),*(hi)=CPUID_EDX,*(lo)=CPUID_EAX
void erdtsc(void);

#ifdef __cplusplus
}
#endif

#endif
