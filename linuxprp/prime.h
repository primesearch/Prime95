/* Handy definitions */

#include "common.h"

/* This controls whether we want to pause computation if the load average */
/* becomes too great.  This does not apply to OS/2. */

#if defined (__linux__) || defined (__FreeBSD__)
#define MPRIME_LOADAVG
#endif

/* Handle the difference between the naming conventions in the two */
/* C compilers.  We only need to do this for global variables that */
/* referenced by the assembly routines */

#define ADD_UNDERSCORES

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

/* The common include files */

#include <time.h>
/*#define SERVER_TESTING*/
extern int NO_GUI;
#include "cpuid.h"
#include "giants.h"
#include "gwnum.h"
#include "prp.h"

/* Global variables */

extern int volatile THREAD_STOP;	/* TRUE if thread should stop */
extern int volatile THREAD_KILL;	/* TRUE if program should terminate */
extern int MENUING;			/* TRUE when main menu active */

/* Internal routines */

void main_menu ();
void linuxContinue (char *);
void Sleep (long);

