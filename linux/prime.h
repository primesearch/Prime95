/* Handy definitions */

#include "common.h"

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

/* Handle the difference between the naming conventions in */
/* C compilers.  We need to do this for global variables that are */
/* referenced by the assembly routines.  Most non-Windows systems */
/* should #define ADD_UNDERSCORES before including files. */

#define ADD_UNDERSCORES

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
#define closesocket	close
#define IsCharAlphaNumeric(c) isalnum(c)

#ifndef __WATCOMC__
#define stricmp(x,y)	strcasecmp(x,y)
#define _timeb		timeb
#define _ftime		ftime
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
extern int NO_GUI;
#include "common.h"
#include "cpuid.h"
#include "gwnum.h"
#include "giants.h"
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

/* Routine definitions */

void rangeStatus (void);
void options_cpu (void);
