/***************************************************************
*
* C file:   speed.h... for cpuinf32 DLL
*
*       This program has been developed by Intel Corporation.  
*		You have Intel's permission to incorporate this code 
*       into your product, royalty free.  Intel has various 
*	    intellectual property rights which it may assert under
*       certain circumstances, such as if another manufacturer's
*       processor mis-identifies itself as being "GenuineIntel"
*		when the CPUID instruction is executed.
*
*       Intel specifically disclaims all warranties, express or
*       implied, and all liability, including consequential and
*		other indirect damages, for the use of this code, 
*		including liability for infringement of any proprietary
*		rights, and including the warranties of merchantability
*		and fitness for a particular purpose.  Intel does not 
*		assume any responsibility for any errors which may 
*		appear in this code nor any responsibility to update it.
*
*  * Other brands and names are the property of their respective
*    owners.
*
*  Copyright (c) 1995, Intel Corporation.  All rights reserved.
***************************************************************/
 
#ifndef speed_h
#define speed_h



// CONSTANT DEFINITIONS ////////////////////////////////////////

	// GENERIC CONSTANTS ///////////////////////////////////////
#define MAXCLOCKS		150		// Maximum number of cycles per
								//   BSF instruction

#define EOA				0		// End of Array variable.

#define CLONE_MASK		0x8000	// Mask to be 'OR'ed with proc-
								//   essor family type

	// ACCURACY AFFECTING CONSTANTS ////////////////////////////
#define ITERATIONS		4000	// Number of times to repeat BSF
								//   instruction in samplings.
								//   Initially set to 4000.

#define INITIAL_DELAY	3		// Number of ticks to wait 
								//   through before starting 
								//   test sampling. Initially
								//   set to 3.
							
#define SAMPLING_DELAY	60		// Number of ticks to allow to
								//   to elapse during sampling.
								//   Initially set to 60.

#define TOLERANCE		1		// Number of MHz to allow
								//   samplings to deviate from
								//   average of samplings.
								//   Initially set to 2.
								
#define MAX_TRIES		20		// Maximum number of samplings
								//   to allow before giving up
								//   and returning current 
								//   average. Initially set to
								//   20.																

#define	SAMPLINGS		10		// Number of BSF sequence 
								//   samplings to make.
								//   Initially set to 10.
								
#define TOL386			2		// Number of MHz above normal-
								//   ised value to normalise
								//   down to the current normal-
								//   ised speed. Initially set
								//   to 2.

#define TOL486			4		// Number of MHz above normal-
								//   ised value to normalise
								//   down to the current normal-
								//   ised speed. Initially set
								//   to 4.

#define TOLP5			5		// Number of MHz above normal-
								//   ised value to normalise
								//   down to the current normal-
								//   ised speed. Initially set
								//   to 5.

#define TOLP6			5		// Number of MHz above normal-
								//   ised value to normalise
								//   down to the current normal-
								//   ised speed. Initially set
								//   to 5.



// VARIABLE STRUCTURE DEFINITIONS //////////////////////////////
struct FREQ_INFO
{
	unsigned long in_cycles;	// Internal clock cycles during
								//   test
								
	unsigned long ex_ticks;		// Microseconds elapsed during 
								//   test
								
	unsigned long raw_freq;		// Raw frequency of CPU in MHz
	
	unsigned long norm_freq;	// Normalized frequency of CPU
								//   in MHz.
};


typedef unsigned short ushort;
typedef unsigned long  ulong;



/***************************************************************
* CpuSpeed() -- Return the raw clock rate of the host CPU.
*
* Inputs:
*	clocks:		NULL: Use default value for number of cycles
*				   per BSF instruction.
*   			Positive Integer: Use clocks value for number
*				   of cycles per BSF instruction.
*
* Returns:
*		If error then return all zeroes in FREQ_INFO structure
*		Else return FREQ_INFO structure containing calculated 
*       clock frequency, normalized clock frequency, number of 
*       clock cycles during test sampling, and the number of 
*       microseconds elapsed during the sampling.
***************************************************************/
struct FREQ_INFO cpuspeed(int clocks);


/***************************************************************
* NormFreq() -- Given an approximate  clock frequency of the
*				host processor, return nearest actual frequency
*				value, e.g. 64 MHz on i486(TM) would return 66 
*				MHz.
*  
* Inputs:
*	processor	Processor type
*	freq		Raw processor clock frequency 
*
* Returns:
*	Normalized processor frequency
***************************************************************/
static ulong NormFreq(ushort processor,ulong freq);


#endif 
