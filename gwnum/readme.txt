Welcome to the gwnum library.

-> how to use
	Before calling one of the gwsetup routines, you must call the
	two cpuid functions:
		#include "cpuid.h"
		guessCpuType ();
		guessCpuSpeed ();
	Alternatively, you can set the global variables defined in cpuid.h

	Next, study gwnum.h for a list of functions.  Also, study prp.c
		for an example of how to use the gwnum library.

-> how to compile / how to port

Windows:
	To make the libraries, your path must already be set to the proper
	C compiler (32-bit or 64-bit).  Then execute
		nmake /f compile
	to build the 32-bit libraries and 
		nmake /f compil64
	to build the 64-bit libraries.

Linux:
	makefile is used to convert the COFF object files built under Windows
	into ELF object files for Linux.  The makefile also compiles the
	necessary C and C++ source files.

	Align the _GWDATA data segment on a 32 byte boundary.  You can link
	with the gwnum.ld file or one of the dummy*.obj files to do this.

	A C application must link with -lstdc++

FreeBSD:
	make.bsd is used to compile the	necessary C and C++ source files.
	For the assembly code, you must copy the ELF .o files from a
	Linux port.

	Align the _GWDATA data segment on a 32 byte boundary.  You can link
	with one of the dummy*.obj files to do this.

	A C application must link with -lcompat and -lstdc++

-> Legal stuff

Copyright (c) 1996-2005, George Woltman

All rights reserved. 

Redistribution and use in source and binary forms, with or without
modification, are permitted provided that the following conditions are 
met: 

(1) Redistributing source code must contain this copyright notice,
limitations, and disclaimer. 
(2) If this software is used to find Mersenne Prime numbers, then
GIMPS will be considered the discoverer of any prime numbers found
and the prize rules at http://mersenne.org/prize.htm will apply.

THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
"AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR
A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT
OWNER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT
LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE,
DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY
THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
(INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE. 
  

