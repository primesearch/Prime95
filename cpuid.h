/*----------------------------------------------------------------------
| Copyright 1995-2003 Just For Fun Software, Inc., all rights reserved
| Author:  George Woltman
| Email: woltman@alum.mit.edu
|
| This file contains routines to determine the CPU type and speed.
+---------------------------------------------------------------------*/

void guessCpuType (void);
void guessCpuSpeed (void);

/* Global variables describing the CPU we are running on */

extern char CPU_BRAND[49];		/* Text description of CPU type */
extern double CPU_SPEED;		/* Actual CPU Speed in MHz */
#define CPU_RDTSC	0x0001
#define CPU_CMOV	0x0002
#define CPU_PREFETCH	0x0004
#define CPU_SSE		0x0008
#define CPU_SSE2	0x0010
#define CPU_MMX		0x0020
EXTERNC unsigned int CPU_FLAGS;		/* Cpu capabilities */
extern int CPU_L1_CACHE_SIZE;
EXTERNC int CPU_L2_CACHE_SIZE;
extern int CPU_L1_CACHE_LINE_SIZE;
EXTERNC int CPU_L2_CACHE_LINE_SIZE;
extern int CPU_L1_DATA_TLBS;
extern int CPU_L2_DATA_TLBS;

/* Internal global variables */

EXTERNC unsigned long CPUID_EAX;	/* For communicating to asm routines */
EXTERNC unsigned long CPUID_EBX;
EXTERNC unsigned long CPUID_ECX;
EXTERNC unsigned long CPUID_EDX;

/* Internal assembly routines */

EXTERNC void fpu_init (void);
EXTERNC unsigned long ecpuidsupport (void);
EXTERNC void ecpuid (void);
#define isCpuidSupported()	ecpuidsupport ()
#define Cpuid(a)		{ CPUID_EAX = a; ecpuid (); }

/* Routine used to time code chunks */
#define rdtsc(hi,lo)	erdtsc(),*(hi)=CPUID_EDX,*(lo)=CPUID_EAX
EXTERNC void erdtsc(void);

