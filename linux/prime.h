/* Handy definitions */

#define FALSE	0
#define TRUE	1

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

/* Handle the difference in the way the two C compilers name routines */

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
#define fpu_init _fpu_init
#define erdtsc	_erdtsc
#define etwo_to_pow _etwo_to_pow
#define esincos	_esincos
#define esincos3 _esincos3
#define etwo_to_pow_over_fftlen _etwo_to_pow_over_fftlen
#define eset_mul_const _eset_mul_const

/* Handle differences between Windows and Linux runtime libraries */

#define stricmp(x,y)	strcasecmp(x,y)
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
#define _timeb		timeb
#define IsCharAlphaNumeric(c) isalnum(c)
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

/* The common include files */

#include <time.h>
/*#define SERVER_TESTING*/
#define EXTERNC
extern int NO_GUI;
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
extern int (*PRIMENET)(short, void *);

/* Internal routines */

void main_menu (void);
void linuxContinue (char *);
void Sleep (long);
void test_user(void);
void test_welcome(void);

/* Assembly routines */

void setupf (void);
int factor64 (void);

