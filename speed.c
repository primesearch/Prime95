/***************************************************************
* C file:  Speed.c... for cpuinf32 DLL
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
 
#include <stdio.h>
#include <math.h>
#include <limits.h>


// Number of cycles needed to execute a single BSF instruction.
//    Note that processors below i386(tm) are not supported.
static ulong processor_cycles[] = {
	00,  00,  00, 115, 47, 43, 
	38,  38,  38, 38,  38, 38, 
};



/***************************************************************
* CpuSpeed() -- Return the raw clock rate of the host CPU.
*
* Inputs:
*	clocks:		0: Use default value for number of cycles
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

struct FREQ_INFO cpuspeed(int clocks) {
	
	ulong  ticks;			// Microseconds elapsed 
					//   during test
	
	ulong  cycles;			// Clock cycles elapsed 
					//   during test
	
	ushort processor = wincpuid();	// Family of processor

	DWORD features = wincpufeatures();	// Features of Processor

	ulong  stamp0, stamp1;			// Time Stamp Variable 
						//   for beginning and end 
						//   of test

	LARGE_INTEGER t0,t1;			// Variables for High-
						//   Resolution Performance
						//   Counter reads
	
	ulong freq  =0;			// Most current frequ. calculation
	ulong freq2 =0;			// 2nd most current frequ. calc.
	ulong freq3 =0;			// 3rd most current frequ. calc.
	
	ulong total;			// Sum of previous three frequency
					//   calculations

	int manual=0;			// Specifies whether the user 
					//   manually entered the number of
					//   cycles for the BSF instruction.

	int tries=0;			// Number of times a calculation has
					//   been made on this call to 
					//   cpuspeed

	struct FREQ_INFO cpu_speed;	// Return structure for cpuspeed
					
	LARGE_INTEGER count_freq;	// High Resolution 
					//   Performance Counter 
					//   frequency
	
	cpu_speed.in_cycles = 0;	// Initialize return structure
	cpu_speed.ex_ticks  = 0;
	cpu_speed.raw_freq  = 0;
	cpu_speed.norm_freq = 0;
	
	if ( processor & CLONE_MASK )
		return cpu_speed;
	  
	// Check for manual BSF instruction clock count
	if (0 == clocks) {
		cycles = ITERATIONS * processor_cycles[processor];
	}
	else if (0 < clocks && clocks <= MAXCLOCKS)  {
		cycles = ITERATIONS * clocks;
		manual = 1;		// Toggle manual control flag.
					//   Note that this mode will not
					// 	 work properly with processors
					//   which can process multiple
					//   BSF instructions at a time.
					//   For example, manual mode
					//   will not work on a 
					//   PentiumPro(R)
	}
	else {				// If invalid clock value
		return cpu_speed;	//    was passed, return
					//    error value
	} 
  	

	if ( !QueryPerformanceFrequency ( &count_freq ) )
					// Checks whether the high-
					//   resolution counter exists
					//   and returns an error if
					//   it does not exist.
		return cpu_speed;

	if ( ( features&0x00000010 ) && !(manual) ) {
					// On processors supporting the Read 
					// Time Stamp opcode, compare elapsed
					// time on the High-Resolution Counter
					// with elapsed cycles on the Time 
					//   Stamp Register.
	
	    do {			// This do loop runs up to 20 times or
	    				//  until the average of the previous 
	    				//  three calculated frequencies is 
	    				//  within 1 MHz of each of the 
	    				//  individual calculated frequencies. 
					//  This resampling increases the 
					//  accuracy of the results since
					//  outside factors could affect this
					//  calculation
			
			tries++;	// Increment number of times sampled
					//   on this call to cpuspeed
			
			freq3 = freq2;	// Shift frequencies back to make
			freq2 = freq;	//   room for new frequency 
					//   measurement

    		QueryPerformanceCounter(&t0);	
    					// Get high-resolution performance 
    					//   counter time
			
			t1.LowPart = t0.LowPart; // Set Initial time
	  		t1.HighPart = t0.HighPart;

    		while ( (ulong)t1.LowPart - (ulong)t0.LowPart<50) {	  
    					// Loop until 50 ticks have 
    					//   passed	since last read of hi-
					//   res counter. This accounts for
					//   overhead later.

				QueryPerformanceCounter(&t1);

			}
			
			
			RDTSC;		// Read Time Stamp

   			_asm
        		{
        		MOV stamp0, EAX
   				}

			t0.LowPart = t1.LowPart;	// Reset Initial 
			t0.HighPart = t1.HighPart;	//   Time

			while ((ulong)t1.LowPart-(ulong)t0.LowPart<1000 ) {
    					// Loop until 1000 ticks have 
    					//   passed since last read of hi-
    					//   res counter. This allows for
					//   elapsed time for sampling.
				QueryPerformanceCounter(&t1);
      			}

			
			RDTSC;		// Read Time Stamp
			
			__asm
	        	{
    	    	MOV stamp1, EAX
        		}

        	cycles = stamp1 - stamp0;	// Number of internal 
        					//   clock cycles is 
        					//   difference between 
        					//   two time stamp 
        					//   readings.

    		ticks = (ulong) t1.LowPart - (ulong) t0.LowPart;	
    						// Number of external ticks is
						//   difference between two
						//   hi-res counter reads.
	
			// Note INTEL's code has been replaced here.

			freq = (unsigned long)
				((double) cycles *
				 (double) count_freq.LowPart /
				 (double) ticks /
					  1000000.0 + 0.5);
          		
			total = ( freq + freq2 + freq3 );
						// Total last three frequency 
						//   calculations
				
		} while ( (tries < 3 ) || 		
		          (tries < 20)&&
		          ((abs(3 * freq -total) > 3*TOLERANCE )||
		           (abs(3 * freq2-total) > 3*TOLERANCE )||
		           (abs(3 * freq3-total) > 3*TOLERANCE )));	
					// Compare last three calculations to 
		          		// average of last three calculations.
		
		
		if ( total / 3  !=  ( total + 1 ) / 3 )
			total ++; 		// Round up if necessary

		freq = total / 3;		// Average last three 
						//   calculations
    
    }

	else if ( processor >= 3 ) {		
					// If processor does not support time 
					//  stamp reading, but is at least a 
					//  386 or above, utilize method of 
					//  timing a loop of BSF instructions 
					//  which take a known number of cycles
					//  to run on i386(tm), i486(tm), and
					//  Pentium(R) processors.
						
		int i;			// Temporary Variable

		ulong current = 0;      // Variable to store time
					//   elapsed during loop of
					//   of BSF instructions

		ulong lowest  = ULONG_MAX;	// Since algorithm finds 
						//   the lowest value out of
						//   a set of samplings, 
						//   this variable is set 
						//   intially to the max 
						//   unsigned long value). 
						//   This guarantees that 
						//   the initialized value 
						//   is not later used as 
						//   the least time through 
						//   the loop.

		for ( i = 0; i < SAMPLINGS; i++ ) { // Sample Ten times. Can
						 //   be increased or 
						 //   decreased depending
						 //   on accuracy vs. time
						 //   requirements

			QueryPerformanceCounter(&t0);	// Get start time

				_asm 
				{
				
					mov eax, 80000000h	
					mov bx, ITERATIONS		
						// Number of consecutive BSF 
						//   instructions to execute. 
						//   Set identical to 
						//   nIterations constant in
						//   speed.h
	           
				loop1:	bsf ecx,eax
               				dec	bx
					jnz	loop1
				}
							 
			QueryPerformanceCounter(&t1);	// Get end time

			current = (ulong) t1.LowPart - (ulong) t0.LowPart;	
						// Number of external ticks is
						//   difference between two
						//   hi-res counter reads.

			if ( current < lowest )		// Take lowest elapsed
				lowest = current;	// time to account
		}					// for some samplings
							// being interrupted
							// by other operations 
		 
		ticks = lowest;


			// Note that some seemingly arbitrary mulitplies and
			//   divides are done below. This is to maintain a 
			//   high level of precision without truncating the 
			//   most significant data. According to what value 
			//   ITERATIIONS is set to, these multiplies and
			//   divides might need to be shifted for optimal
			//   precision.

		ticks = ticks * 100000;	
						// Convert ticks to hundred
						//   thousandths of a tick
			
		ticks = ticks / ( count_freq.LowPart/10 );		
						// Hundred Thousandths of a 
						// Ticks / ( 10 ticks/second )
						//   = microseconds (us)
		
		if ( ticks%count_freq.LowPart > count_freq.LowPart/2 )	
			ticks++;		// Round up if necessary
			
		freq = cycles/ticks;		// Cycles / us  = MHz
		
	    if ( cycles%ticks > ticks/2 )
    		freq++;				// Round up if necessary	
          	
	}

    else {				// In case of 286 or below class
    					// processors, return all zeroes. This
    	cycles = 0;		//   code does not support those
    	ticks  = 0;		//   processors.
    	freq   = 0;
    
    }

	cpu_speed.in_cycles = cycles;	// Return variable structure
	cpu_speed.ex_ticks  = ticks;	//   determined by one of 
	cpu_speed.raw_freq  = freq;	//   the algorithms above
	cpu_speed.norm_freq = NormFreq(processor, freq);

	return cpu_speed;
   	
} // cpuspeed()



 
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

static ulong NormFreq(ushort processor,ulong freq) {
	ushort *speeds;
	
	ushort i386Speeds[] = { 16, 20, 25, 33, 40, EOA };
	ushort i486Speeds[] = { 25, 33, 50, 66, 75, 100, EOA };
	ushort iPentiumSpeeds[]  = { 60, 66, 75, 90, 100, 120, 133,
			150, 166, 185, 200, 233, 266, 300, EOA };  
	ushort iPentiumProSpeeds[]  = { 133, 150, 167, 185, 200,
			220, 233, 240, 266, 300, 333, 350, 366, 400, 433,
			450, 466, 500, 533, 550, 600, 650, 700, 733, 750,
			800, 866, 933, 1000,
			1050, 1066, 1100, 1133, 1150, 1200, 1250, 1266,
			1300, 1333, 1350, 1400, 1466, 1500, EOA };

	int ptr = 0;

	if (3 == processor) {                
		
		speeds = i386Speeds;
		
		while (speeds[ptr] != EOA) {
			if ( (int)freq <= (int)(speeds[ptr] + TOL386))
				return  speeds[ptr];
					// Scan each speed in array until
					//  current calculated frequency
					//  is less than or equal to 
					//  normalized value plus TOL386.
		ptr++;
		}
		
		return freq;		// If raw speed is higher than 
					//   highest normalized speed plus
					//   TOL386, return raw frequency.
	}
	else if (4 == processor) {
		
		speeds = i486Speeds;
		
		while (speeds[ptr] != EOA) {
			if ( (int)freq <= (int) (speeds[ptr] + TOL486))
				return  speeds[ptr];
					// Scan each speed in array until
					//  current calculated frequency
					//  is less than or equal to 
					//  normalized value plus TOL486.
		ptr++;
		}
		
		return freq;		// If raw speed is higher than 
					//   highest normalized speed plus
					//   TOL486, return raw frequency.
	}
	else if (5 == processor) {

		speeds = iPentiumSpeeds;
		
		while (speeds[ptr] != EOA) {
			if ( (int)freq <= (int) (speeds[ptr] + TOLP5))
				return  speeds[ptr];
					// Scan each speed in array until
					//  current calculated frequency
					//  is less than or equal to 
					//  normalized value plus TOLP5.
		ptr++;
		}
		
		return freq;		// If raw speed is higher than 
					//   highest normalized speed plus
					//   TOLP5, return raw frequency.
	}
	else if (6 == processor || 8 == processor ||
		 9 == processor || 10 == processor) {

		speeds = iPentiumProSpeeds;
		
		while (speeds[ptr] != EOA) {
			if ( (int)freq <= (int) (speeds[ptr] + TOLP6))
				return  speeds[ptr];
					// Scan each speed in array until
					//  current calculated frequency
					//  is less than or equal to 
					//  normalized value plus TOLP6.
		ptr++;
		}
		
		return freq;		// If raw speed is higher than 
					//   highest normalized speed plus
					//   TOLP6, return raw frequency.
	}
	else {
		return freq;		// return raw frequency 
	}
} // NormFreq()
