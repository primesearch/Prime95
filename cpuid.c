/*----------------------------------------------------------------------
| Copyright 1995-2002 Just For Fun Software, Inc., all rights reserved
| Author:  George Woltman
| Email: woltman@alum.mit.edu
|
| This file contains routines to determine the CPU type and speed.
+---------------------------------------------------------------------*/

/* Global variables describing the CPU we are running on */

char	CPU_BRAND[49] = "";
double	CPU_SPEED = 100.0;
EXTERNC unsigned int CPU_FLAGS = 0;
int	CPU_L1_CACHE_SIZE = -1;
int	CPU_L2_CACHE_SIZE = -1;
int	CPU_L1_CACHE_LINE_SIZE = -1;
int	CPU_L2_CACHE_LINE_SIZE = -1;
int	CPU_L1_DATA_TLBS = -1;
int	CPU_L2_DATA_TLBS = -1;

/* Internal global variables */

EXTERNC unsigned long CPUID_EAX = 0;	/* For communicating to asm routines */
EXTERNC unsigned long CPUID_EBX = 0;
EXTERNC unsigned long CPUID_ECX = 0;
EXTERNC unsigned long CPUID_EDX = 0;

/* Internal routines to see if CPU-specific instructions (RDTSC, CMOV */
/* SSE, SSE2) are supported.  CPUID could report them as supported yet */
/* the OS might not support them. */

#if defined (__linux__) || defined (__FreeBSD__) || defined (__EMX__)

#include <setjmp.h>
int	boom;
jmp_buf	env;
void sigboom_handler (int i)
{
	boom = TRUE;
	longjmp (env, 1);
}
int canExecInstruction (
	unsigned long cpu_flag)
{
	boom = FALSE;
	(void) signal (SIGILL, sigboom_handler);
	if (setjmp (env) == 0) {
		switch (cpu_flag) {
		case CPU_RDTSC:		/* RDTSC */
			__asm__ __volatile__ (".byte 0x0F\n .byte 0x31\n");
			break;
		case CPU_CMOV:		/* CMOV */
			__asm__ __volatile__ (".byte 0x0F\n .byte 0x42\n .byte 0xC0\n");
			break;
		case CPU_MMX:		/* PADDB */
			__asm__ __volatile__ (".byte 0x0F\n .byte 0xFC\n .byte 0xC0\n");
			break;
		case CPU_SSE:		/* ORPS */
			__asm__ __volatile__ (".byte 0x0F\n .byte 0x56\n .byte 0xC0\n");
			break;
		case CPU_SSE2:		/* ADDPD */
			__asm__ __volatile__ (".byte 0x66\n .byte 0x0F\n .byte 0x58\n .byte 0xC0\n");
			break;
		case CPU_PREFETCH:	/* PREFETCHT1 */
			__asm__ __volatile__ (".byte 0x0F\n .byte 0x18\n .byte 0x16\n");
			break;
		}
	}
	(void) signal (SIGILL, SIG_DFL);
	fpu_init ();
	return (!boom);
}

#elif defined (__IBMC__)

int canExecInstruction (
	unsigned long cpu_flag)
{
	return TRUE;
}

#else

int canExecInstruction (
	unsigned long cpu_flag)
{
	int	succeeded;
	__try {
		switch (cpu_flag) {
		case CPU_RDTSC:		/* RDTSC */
			__asm __emit 0x0F
			__asm __emit 0x31
			break;
		case CPU_CMOV:		/* CMOV */
			__asm __emit 0x0F
			__asm __emit 0x42
			__asm __emit 0xC0
			break;
		case CPU_MMX:		/* PADDB */
			__asm __emit 0x0F
			__asm __emit 0xFC
			__asm __emit 0xC0
			break;
		case CPU_SSE:		/* ORPS */
			__asm __emit 0x0F
			__asm __emit 0x56
			__asm __emit 0xC0
			break;
		case CPU_SSE2:		/* ADDPD */
			__asm __emit 0x66
			__asm __emit 0x0F
			__asm __emit 0x58
			__asm __emit 0xC0
			break;
		case CPU_PREFETCH:	/* PREFETCHT1 */
			__asm __emit 0x0F
			__asm __emit 0x18
			__asm __emit 0x16
			break;
		}
		succeeded = TRUE;
	}
	__except (EXCEPTION_EXECUTE_HANDLER) {
		succeeded = FALSE;
	}
	fpu_init ();
	return (succeeded);
}
#endif


/* Work with CPUID instruction to guess the cpu type and features */
/* See Intel's document AP-485 for using CPUID on Intel processors */
/* AMD and VIA have similar documents */

void guessCpuType (void)
{
	unsigned long max_cpuid_value;
	unsigned long max_extended_cpuid_value;
	unsigned long extended_family, extended_model, type, family_code;
	unsigned long model_number, stepping_id, brand_index;
	char	vendor_id[13];

/* Set up default values for features we cannot determine with CPUID */

	CPU_BRAND[0] = 0;
	CPU_SPEED = 100.0;
	CPU_FLAGS = 0;
	CPU_L1_CACHE_SIZE = -1;
	CPU_L2_CACHE_SIZE = -1;
	CPU_L1_CACHE_LINE_SIZE = -1;
	CPU_L2_CACHE_LINE_SIZE = -1;
	CPU_L1_DATA_TLBS = -1;
	CPU_L2_DATA_TLBS = -1;

/* If CPUID instruction is not supported, assume we have a 486 (not all */
/* 486 chips supported CPUID.  The CPU might be something else, but that */
/* isn't particularly important. */

	if (! isCpuidSupported ()) {
		strcpy (CPU_BRAND, "CPUID not supported - 486 CPU assumed");
		return;
	}

/* Call CPUID with 0 argument.  It returns how highest argument CPUID */
/* can accept as well as the vendor string */

	Cpuid (0);
	max_cpuid_value = CPUID_EAX;
	memcpy (vendor_id, &CPUID_EBX, 4);
	memcpy (vendor_id+4, &CPUID_EDX, 4);
	memcpy (vendor_id+8, &CPUID_ECX, 4);
	vendor_id[12] = 0;

/* So far all vendors have adopted Intel's definition of CPUID with 1 as an */
/* argument.  Let's assume future vendors will do the same.  CPUID returns */
/* the processor family, stepping, etc.  It also returns the feature flags. */

	if (max_cpuid_value >= 1) {
		Cpuid (1);
		extended_family = (CPUID_EAX >> 20) & 0xFF;
		extended_model = (CPUID_EAX >> 16) & 0xF;
		type = (CPUID_EAX >> 12) & 0x3;
		family_code = (CPUID_EAX >> 8) & 0xF;
		model_number = (CPUID_EAX >> 4) & 0xF;
		stepping_id = CPUID_EAX & 0xF;
		brand_index = CPUID_EBX & 0xFF;
		if ((CPUID_EDX >> 4) & 0x1 && canExecInstruction (CPU_RDTSC))
			CPU_FLAGS |= CPU_RDTSC;
		if ((CPUID_EDX >> 15) & 0x1 && canExecInstruction (CPU_CMOV))
			CPU_FLAGS |= CPU_CMOV;
		if ((CPUID_EDX >> 23) & 0x1 && canExecInstruction (CPU_MMX))
			CPU_FLAGS |= CPU_MMX;
		if ((CPUID_EDX >> 25) & 0x1 && canExecInstruction (CPU_PREFETCH))
			CPU_FLAGS |= CPU_PREFETCH;
		if ((CPUID_EDX >> 25) & 0x1 && canExecInstruction (CPU_SSE))
			CPU_FLAGS |= CPU_SSE;
		if ((CPUID_EDX >> 26) & 0x1 && canExecInstruction (CPU_SSE2))
			CPU_FLAGS |= CPU_SSE2;
	}

/* Call CPUID with 0x80000000 argument.  It tells us how many extended CPU */
/* functions are supported. */

	Cpuid (0x80000000);
	max_extended_cpuid_value = CPUID_EAX;

/* Although not guaranteed, all vendors have standardized on putting the */
/* brand string (if supported) at cpuid calls 0x8000002, 0x80000003, and */
/* 0x80000004.  We'll assume future vendors will do the same. */

	if (max_extended_cpuid_value >= 80000004) {
		Cpuid (0x80000002);
		memcpy (CPU_BRAND, &CPUID_EAX, 4);
		memcpy (CPU_BRAND+4, &CPUID_EBX, 4);
		memcpy (CPU_BRAND+8, &CPUID_ECX, 4);
		memcpy (CPU_BRAND+12, &CPUID_EDX, 4);
		Cpuid (0x80000003);
		memcpy (CPU_BRAND+16, &CPUID_EAX, 4);
		memcpy (CPU_BRAND+20, &CPUID_EBX, 4);
		memcpy (CPU_BRAND+24, &CPUID_ECX, 4);
		memcpy (CPU_BRAND+28, &CPUID_EDX, 4);
		Cpuid (0x80000004);
		memcpy (CPU_BRAND+32, &CPUID_EAX, 4);
		memcpy (CPU_BRAND+36, &CPUID_EBX, 4);
		memcpy (CPU_BRAND+40, &CPUID_ECX, 4);
		memcpy (CPU_BRAND+44, &CPUID_EDX, 4);
		CPU_BRAND[48] = 0;
		while (CPU_BRAND[0] == ' ') strcpy (CPU_BRAND, CPU_BRAND+1);
	}

/*-------------------------------------------------------------------+
| Check for INTEL vendor string.  Perform INTEL-specific operations. |
+-------------------------------------------------------------------*/

	if (strcmp ((const char *) vendor_id, "GenuineIntel") == 0) {

/* If we haven't figured out the brand string, create one based on the */
/* brand_index as recommended in the AP-485 document. */

		if (CPU_BRAND[0] == 0) {
			switch (brand_index) {
			case 1:
				strcpy (CPU_BRAND, "Intel(R) Celeron(R) processor");
				break;
			case 2:
			case 4:
				strcpy (CPU_BRAND, "Intel(R) Pentium(R) III processor");
				break;
			case 3:
				if (family_code != 6 ||
				    model_number != 0xB ||
				    stepping_id != 1)
					strcpy (CPU_BRAND, "Intel(R) Pentium(R) III Xeon processor");
				else
					strcpy (CPU_BRAND, "Intel(R) Celeron(R) processor");
				break;
			case 6:
				strcpy (CPU_BRAND, "Mobile Intel(R) Pentium(R) III processor");
				break;
			case 7:
				strcpy (CPU_BRAND, "Mobile Intel(R) Celeron(R) processor");
				break;
			case 8:
				if ((family_code << 8) +
				    (model_number << 4) +
				    stepping_id < 0xF20)
					strcpy (CPU_BRAND, "Intel(R) Pentium(R) 4 processor");
				else
					strcpy (CPU_BRAND, "Intel(R) Genuine processor");
				break;
			case 9:
				strcpy (CPU_BRAND, "Intel(R) Pentium(R) 4 processor");
				break;
			case 0xB:
			case 0xE:
				strcpy (CPU_BRAND, "Intel(R) Xeon processor");
				break;
			}
		}

/* Call CPUID with 2 argument.  It returns the cache size and structure */
/* in a series of 8-bit descriptors */

		if (max_cpuid_value >= 2) {
			Cpuid (2);
			if ((CPUID_EAX & 0xFF) > 0) {
				unsigned int descriptors[15];
				int i, count;
				count = 0;
				if (! (CPUID_EAX & 0x80000000)) {
					descriptors[count++] = (CPUID_EAX >> 24) & 0xFF;
					descriptors[count++] = (CPUID_EAX >> 16) & 0xFF;
					descriptors[count++] = (CPUID_EAX >> 8) & 0xFF;
				}
				if (! (CPUID_EBX & 0x80000000)) {
					descriptors[count++] = (CPUID_EBX >> 24) & 0xFF;
					descriptors[count++] = (CPUID_EBX >> 16) & 0xFF;
					descriptors[count++] = (CPUID_EBX >> 8) & 0xFF;
					descriptors[count++] = CPUID_EBX & 0xFF;
				}
				if (! (CPUID_ECX & 0x80000000)) {
					descriptors[count++] = (CPUID_ECX >> 24) & 0xFF;
					descriptors[count++] = (CPUID_ECX >> 16) & 0xFF;
					descriptors[count++] = (CPUID_ECX >> 8) & 0xFF;
					descriptors[count++] = CPUID_ECX & 0xFF;
				}
				if (! (CPUID_EDX & 0x80000000)) {
					descriptors[count++] = (CPUID_EDX >> 24) & 0xFF;
					descriptors[count++] = (CPUID_EDX >> 16) & 0xFF;
					descriptors[count++] = (CPUID_EDX >> 8) & 0xFF;
					descriptors[count++] = CPUID_EDX & 0xFF;
				}
				for (i = 0; i < count; i++) {
					switch (descriptors[i]) {
					case 0x03:
						CPU_L2_DATA_TLBS = 64;
						break;
					case 0x0A:
						CPU_L1_CACHE_SIZE = 8;
						CPU_L1_CACHE_LINE_SIZE = 32;
						break;
					case 0x0C:
						CPU_L1_CACHE_SIZE = 16;
						CPU_L1_CACHE_LINE_SIZE = 32;
						break;
					case 0x40:
						if (family_code == 15) {
							/* no L3 cache */
						} else {
							CPU_L2_CACHE_SIZE = 0;
						}
						break;
					case 0x41:
						CPU_L2_CACHE_SIZE = 128;
						CPU_L2_CACHE_LINE_SIZE = 32;
						break;
					case 0x42:
						CPU_L2_CACHE_SIZE = 256;
						CPU_L2_CACHE_LINE_SIZE = 32;
						break;
					case 0x43:
						CPU_L2_CACHE_SIZE = 512;
						CPU_L2_CACHE_LINE_SIZE = 32;
						break;
					case 0x44:
						CPU_L2_CACHE_SIZE = 1024;
						CPU_L2_CACHE_LINE_SIZE = 32;
						break;
					case 0x45:
						CPU_L2_CACHE_SIZE = 2048;
						CPU_L2_CACHE_LINE_SIZE = 32;
						break;
					case 0x5B:
						CPU_L2_DATA_TLBS = 64;
						break;
					case 0x5C:
						CPU_L2_DATA_TLBS = 128;
						break;
					case 0x5D:
						CPU_L2_DATA_TLBS = 256;
						break;
					case 0x66:
						CPU_L1_CACHE_SIZE = 8;
						CPU_L1_CACHE_LINE_SIZE = 64;
						break;
					case 0x67:
						CPU_L1_CACHE_SIZE = 16;
						CPU_L1_CACHE_LINE_SIZE = 64;
						break;
					case 0x68:
						CPU_L1_CACHE_SIZE = 32;
						CPU_L1_CACHE_LINE_SIZE = 64;
						break;
					case 0x39:
					case 0x79:
						CPU_L2_CACHE_SIZE = 128;
						CPU_L2_CACHE_LINE_SIZE = 64;
						break;
					case 0x3C:
					case 0x7A:
						CPU_L2_CACHE_SIZE = 256;
						CPU_L2_CACHE_LINE_SIZE = 64;
						break;
					case 0x7B:
						CPU_L2_CACHE_SIZE = 512;
						CPU_L2_CACHE_LINE_SIZE = 64;
						break;
					case 0x7C:
						CPU_L2_CACHE_SIZE = 1024;
						CPU_L2_CACHE_LINE_SIZE = 64;
						break;
					case 0x82:
						CPU_L2_CACHE_SIZE = 256;
						CPU_L2_CACHE_LINE_SIZE = 32;
						break;
					case 0x83:
						CPU_L2_CACHE_SIZE = 512;
						CPU_L2_CACHE_LINE_SIZE = 32;
						break;
					case 0x84:
						CPU_L2_CACHE_SIZE = 1024;
						CPU_L2_CACHE_LINE_SIZE = 32;
						break;
					case 0x85:
						CPU_L2_CACHE_SIZE = 2048;
						CPU_L2_CACHE_LINE_SIZE = 32;
						break;
					}
				}
			}
		}

/* Deduce the cpu type given the family, model, stepping, etc.  If we */
/* haven't figured out the brand string, create one based on the cpu type. */

		if (family_code == 4) {
			strcpy (CPU_BRAND, "Intel 486 processor");
		}
		if (family_code == 5) {
			if (type == 0 && model_number <= 2)
				strcpy (CPU_BRAND, "Intel Pentium processor");
			if (type == 1 && model_number <= 3)
				strcpy (CPU_BRAND, "Intel Pentium OverDrive processor");
			if (type == 0 && model_number >= 4)
				strcpy (CPU_BRAND, "Intel Pentium MMX processor");
			if (type == 1 && model_number >= 4)
				strcpy (CPU_BRAND, "Intel Pentium MMX OverDrive processor");
		}
		if (family_code == 6 && model_number == 1) {
			strcpy (CPU_BRAND, "Intel Pentium Pro processor");
		}
		if (family_code == 6 && model_number <= 6) {
			if (type == 0 && model_number == 3)
				strcpy (CPU_BRAND, "Intel Pentium II processor");
			if (model_number == 5 && CPU_L2_CACHE_SIZE == 0 ||
			    model_number == 6) {
				strcpy (CPU_BRAND, "Intel Celeron processor");
			}
			if (model_number == 5 && CPU_L2_CACHE_SIZE == 512)
				strcpy (CPU_BRAND, "Intel Pentium II or Pentium II Xeon processor");
			if (model_number == 5 && CPU_L2_CACHE_SIZE >= 1024)
				strcpy (CPU_BRAND, "Intel Pentium II Xeon processor");
			if (type == 1 && model_number == 3)
				strcpy (CPU_BRAND, "Intel Pentium II OverDrive processor");
		}
		if (family_code == 6 && model_number >= 7) {
			if (model_number == 7 && CPU_L2_CACHE_SIZE >= 1024)
				strcpy (CPU_BRAND, "Intel Pentium III Xeon processor");
			if (model_number == 7 && CPU_L2_CACHE_SIZE == 512)
				strcpy (CPU_BRAND, "Intel Pentium III or Pentium III Xeon processor");
		}

/* If we've failed to figure out the brand string, create a default. */

		if (CPU_BRAND[0] == 0) {
			strcpy (CPU_BRAND, "Unknown Intel CPU");
		}
	}

/*---------------------------------------------------------------+
| Check for AMD vendor string.  Perform AMD-specific operations. |
+---------------------------------------------------------------*/

	else if (strcmp ((const char *) vendor_id, "AuthenticAMD") == 0) {

/* Deduce the cpu type given the family, model, stepping, etc.  If we */
/* haven't figured out the brand string, create one based on the cpu type. */

		if (family_code == 4) {
			strcpy (CPU_BRAND, "AMD Am486 or Am5x86 processor");
		}
		if (family_code == 5 && model_number <= 3) {
			if (CPU_BRAND[0] == 0)
				strcpy (CPU_BRAND, "AMD K5 processor");
		}

/* Early Athlon CPUs support the SSE prefetch instructions even though */
/* they do not support the full SSE instruction set.  I think testing for */
/* the AMD MMX extensions capability will detect this case. */

		if (max_extended_cpuid_value >= 80000001 &&
		    ! (CPU_FLAGS & CPU_PREFETCH)) {
			Cpuid (0x80000001);
			if ((CPUID_EDX >> 22) & 0x1 &&
			    canExecInstruction (CPU_PREFETCH))
				CPU_FLAGS |= CPU_PREFETCH;
		}

/* Get the L1 cache size and number of data TLBs */

		if (max_extended_cpuid_value >= 80000005) {
			Cpuid (0x80000005);
			CPU_L1_DATA_TLBS = (CPUID_EBX >> 16) & 0xFF;
			CPU_L1_CACHE_SIZE = (CPUID_ECX >> 24) & 0xFF;
			CPU_L1_CACHE_LINE_SIZE = CPUID_ECX & 0xFF;
		}

/* Get the L2 cache size */

		if (max_extended_cpuid_value >= 80000006) {
			Cpuid (0x80000006);
			CPU_L2_DATA_TLBS = (CPUID_EBX >> 16) & 0xFFF;
			CPU_L2_CACHE_SIZE = (CPUID_ECX >> 16) & 0xFFFF;
			CPU_L2_CACHE_LINE_SIZE = CPUID_ECX & 0xFF;
		}

/* If we haven't figured out the brand string, create a default one */

		if (CPU_BRAND[0] == 0)
			strcpy (CPU_BRAND, "Unknown AMD CPU");
	}

/*---------------------------------------------------------------+
| Check for VIA vendor string.  Perform VIA-specific operations. |
+---------------------------------------------------------------*/

	else if (strcmp ((const char *) vendor_id, "CentaurHauls") == 0) {

/* Get the L1 cache size and number of data TLBs */

		if (max_extended_cpuid_value >= 80000005) {
			Cpuid (0x80000005);
			CPU_L2_DATA_TLBS = (CPUID_EBX >> 16) & 0xFF;
			CPU_L1_CACHE_SIZE = (CPUID_ECX >> 24) & 0xFF;
			CPU_L1_CACHE_LINE_SIZE = CPUID_ECX & 0xFF;
		}

/* Get the L2 cache size */

		if (max_extended_cpuid_value >= 80000006) {
			Cpuid (0x80000006);
			CPU_L2_CACHE_SIZE = (CPUID_ECX >> 24) & 0xFF;
			CPU_L2_CACHE_LINE_SIZE = CPUID_ECX & 0xFF;
		}

/* If we haven't figured out the brand string, create a default one */

		if (CPU_BRAND[0] == 0)
			strcpy (CPU_BRAND, "Unknown VIA/CYRIX CPU");
	}

/*--------------------------------------------------------+
| An unknown CPU vendor.  Fill in defaults as best we can |
+--------------------------------------------------------*/

	else {
		if (CPU_BRAND[0] == 0) {
			strcpy (CPU_BRAND, "Unrecognized CPU vendor: ");
			strcat (CPU_BRAND, vendor_id);
		}
	}
}

void guessCpuSpeed (void)
{

/* If RDTSC is not supported, then measuring the CPU speed is real hard */
/* so we just assume the CPU speed is 100 MHz.  This isn't a big deal */
/* since CPUs not supporting RDTSC aren't powerful enough to run prime95. */

	if (! (CPU_FLAGS & CPU_RDTSC)) {
		CPU_SPEED = 100.0;
		return;
	}

/* If this machine supports a high resolution counter, use that for timing */

	if (isHighResTimerAvailable ()) {
		unsigned long start_hi, start_lo, end_hi, end_lo;
		double	frequency, temp, start_time, end_time;
		unsigned long iterations;
		double	speed1, speed2, speed3, avg_speed;
		int	tries;

/* Compute the number of high resolution ticks in one millisecond */
/* This should give us good accuracy while hopefully avoiding time slices */

		frequency = getHighResTimerFrequency ();
		iterations = (unsigned long) (frequency / 1000.0);

/* Do up to 20 iterations (idea lifted from Intel's code) until 3 straight */
/* speed calculations are within 1 MHz of each other.  This is good since */
/* outside forces can interfere with this calculation. */

		tries = 0; speed1 = 0.0; speed2 = 0.0;
		do {

/* Shuffle the last calculations, bump counter */

			speed3 = speed2;
			speed2 = speed1;
			tries++;

/* Loop waiting for high resolution timer to change */

			temp = getHighResTimer ();
			while ((start_time = getHighResTimer ()) == temp);
			rdtsc (&start_hi, &start_lo);

/* Now loop waiting for timer to tick off about 1 millisecond */

			temp = start_time + (double) iterations;
			while ((end_time = getHighResTimer ()) < temp);
			rdtsc (&end_hi, &end_lo);

/* Compute speed based on number of clocks in the time interval */

			speed1 = (end_hi * 4294967296.0 + end_lo -
				  start_hi * 4294967296.0 - start_lo) *
				 frequency /
				 (end_time - start_time) / 1000000.0;

/* Caclulate average of last 3 speeds.  Loop if this average isn't */
/* very close to all of the last three speed calculations. */

			avg_speed = (speed1 + speed2 + speed3) / 3.0;
		} while (tries < 3 ||
		         (tries < 20 &&
		          (fabs (speed1 - avg_speed) > 1.0 ||
		           fabs (speed2 - avg_speed) > 1.0 ||
		           fabs (speed3 - avg_speed) > 1.0)));

/* Final result is average speed of last three calculations */

		CPU_SPEED = avg_speed;
	}

/* Otherwise use the low resolution timer to measure CPU speed */

	else {
		struct _timeb temp, start_time, end_time;
		unsigned long start_hi, start_lo, end_hi, end_lo;
		double	speed1, speed2, speed3, avg_speed, elapsed_time;
		int	tries;

/* Do up to 10 iterations until 3 straight speed calculations are within */
/* 1 MHz of each other. */

		tries = 0; speed1 = 0.0; speed2 = 0.0;
		do {

/* Shuffle the last calculations, bump counter */

			speed3 = speed2;
			speed2 = speed1;
			tries++;

/* Loop waiting for low resolution timer to change */

			_ftime (&temp);
			do
				_ftime (&start_time);
			while (temp.millitm == start_time.millitm);
			rdtsc (&start_hi, &start_lo);

/* Now loop waiting for timer to change again */

			do
				_ftime (&end_time);
			while (start_time.millitm == end_time.millitm);
			rdtsc (&end_hi, &end_lo);

/* Compute elapsed time.  Since most PCs have a low resolution clock */
/* that ticks every 18.20648193 seconds, then if elapsed time is close */
/* to 1 / 18.2 seconds, then assume the elapsed time is */
/* 1 / 18.206... = 0.054925493 seconds. */

			elapsed_time = (end_time.time - start_time.time) +
				  ((int) end_time.millitm -
				   (int) start_time.millitm) / 1000.0;
			if (elapsed_time >= 0.049 && elapsed_time <= 0.061)
				elapsed_time = 0.054925493;

/* Compute speed based on number of clocks in the time interval */

			speed1 = (end_hi * 4294967296.0 + end_lo -
				  start_hi * 4294967296.0 - start_lo) /
				 elapsed_time / 1000000.0;

/* Caclulate average of last 3 speeds.  Loop if this average isn't */
/* very close to all of the last three speed calculations. */

			avg_speed = (speed1 + speed2 + speed3) / 3.0;
		} while (tries < 3 ||
		         (tries < 10 &&
		          (fabs (speed1 - avg_speed) > 1.0 ||
		           fabs (speed2 - avg_speed) > 1.0 ||
		           fabs (speed3 - avg_speed) > 1.0)));

/* Final result is average speed of last three calculations */

		CPU_SPEED = avg_speed;
	}
}
