/* Handy definitions */

#define FALSE	0
#define TRUE	1

/* Port number used in version numbers and result reporting. */

#ifdef __linux__
#define PORT	2
#endif
#ifdef __FreeBSD__
#define PORT	6
#endif
#if defined (__EMX__) || defined (__IBMC__) || defined (__OS2__)
#define PORT	7
#endif

/* This controls whether we want to pause computation if the load average */
/* becomes too great.  This does not apply to OS/2. */

#if defined (__linux__) || defined (__FreeBSD__)
#define MPRIME_LOADAVG
#endif

/* Handle the difference between the naming conventions in the two */
/* C compilers.  We only need to to this for global variables that */
/* referenced by the assembly routines */

#define ERRCHK	 _ERRCHK
#define CPU_TYPE _CPU_TYPE
#define CPU_L2_CACHE_SIZE _CPU_L2_CACHE_SIZE
#define CPU_L2_CACHE_LINE_SIZE _CPU_L2_CACHE_LINE_SIZE
#define MAXERR	 _MAXERR
#define GWERROR	 _GWERROR
#define PARG	 _PARG
#define FACLSW	 _FACLSW
#define FACMSW	 _FACMSW
#define FACHSW	 _FACHSW
#define FACPASS	 _FACPASS

#define FFTLEN	_FFTLEN
#define SRCARG	_SRCARG
#define SRC2ARG	_SRC2ARG
#define DESTARG	_DESTARG
#define DEST2ARG _DEST2ARG
#define MAXDIFF _MAXDIFF
#define INFP	_INFP
#define INFF	_INFF
#define INFT	_INFT
#define GWPROCPTRS _GWPROCPTRS
#define UU	_UU
#define VV	_VV
#define Ulen	_Ulen
#define Vlen	_Vlen
#define EGCD_A	_EGCD_A
#define EGCD_B	_EGCD_B
#define EGCD_C	_EGCD_C
#define EGCD_D	_EGCD_D
#define EGCD_ODD _EGCD_ODD
#define NORM_ARRAY1 _NORM_ARRAY1
#define NORM_ARRAY2 _NORM_ARRAY2
#define RES	_RES
#define CARRYL	_CARRYL
#define CARRYH	_CARRYH
#define PLUS1	_PLUS1

#define NORMRTN		_NORMRTN
#define FFTZERO		_FFTZERO
#define CPU_FLAGS	_CPU_FLAGS
#define NUMLIT		_NUMLIT
#define NUMBIG		_NUMBIG
#define BITS_PER_WORD	_BITS_PER_WORD
#define FFTLEN_INV	_FFTLEN_INV
#define PROTHVALS 	_PROTHVALS
#define ADDIN_ROW	_ADDIN_ROW
#define ADDIN_OFFSET	_ADDIN_OFFSET
#define ADDIN_VALUE	_ADDIN_VALUE
#define COPYZERO	_COPYZERO
#define POSTFFT		_POSTFFT

#define CPUID_EAX	_CPUID_EAX
#define CPUID_EBX	_CPUID_EBX
#define CPUID_ECX	_CPUID_ECX
#define CPUID_EDX	_CPUID_EDX

/* Handle the difference in the way the two C compilers name routines */

#define timeit _timeit
#define setupf	 _setupf
#define factor64 _factor64
#define gw_setup1 _gw_setup1
#define gw_setup2 _gw_setup2
#define egetval _egetval
#define egcdhlp _egcdhlp
#define eaddhlp	_eaddhlp
#define emuladdhlp _emuladdhlp
#define emulsubhlp _emulsubhlp
#define emuladd2hlp _emuladd2hlp
#define gwinfo1	_gwinfo1
#define gwsetup1 _gwsetup1
#define gwsetup2 _gwsetup2
#define esubhlp _esubhlp
#define emulmod	_emulmod
#define eisvaliddouble _eisvaliddouble
#define fpu_init _fpu_init
#define erdtsc	_erdtsc
#define etwo_to_pow _etwo_to_pow
#define esincos	_esincos
#define esincos3 _esincos3
#define etwo_to_pow_over_fftlen _etwo_to_pow_over_fftlen
#define eset_mul_const _eset_mul_const
#define ecpuidsupport _ecpuidsupport
#define ecpuid _ecpuid

/* Handle differences between Windows and Linux runtime libraries */

#define _commit(f)	fsync(f)
#define _open		open
#define _close		close
#define _read		read
#define _write		write
#define _lseek		lseek
#define _unlink		unlink
#define _creat		creat
#define _chdir		chdir
#define _ftime		ftime
#define closesocket	close
#define IsCharAlphaNumeric(c) isalnum(c)

#ifndef __WATCOMC__
#define stricmp(x,y)	strcasecmp(x,y)
#define _timeb		timeb
#define _O_APPEND	O_APPEND
#define _O_RDONLY	O_RDONLY
#define _O_WRONLY	O_WRONLY
#define _O_RDWR		O_RDWR
#define _O_CREAT	O_CREAT
#define _O_TRUNC	O_TRUNC
#define _O_BINARY 	0
#define _O_TEXT		0

/* Handle differences between Windows and OS/2 runtime libraries */

#ifdef __IBMC__
#define stricmp(x,y)  stricmp(x,y)
#define _commit(f)    /* no commit/fsync on OS/2 */
#define _ftime        _ftime
#endif
#endif

/* The common include files */

#include <time.h>
/*#define SERVER_TESTING*/
#ifdef __WATCOMC__
#define EXTERNC extern
#else
#define EXTERNC
#endif
extern int NO_GUI;
#include "cpuid.h"
#include "giants.h"
#include "gwnum.h"
#include "commona.h"
#include "commonc.h"
#include "commonb.h"
#include "primenet.h"

/* Global variables */

extern int volatile THREAD_STOP;	/* TRUE if thread should stop */
extern int volatile THREAD_KILL;	/* TRUE if program should terminate */
extern int MENUING;			/* TRUE when main menu active */

/* Internal routines */

void main_menu (void);
void linuxContinue (char *);
void Sleep (long);
void test_user(void);
void test_welcome(void);

/* Assembly routines */

void setupf (void);
int factor64 (void);

