/* These macros produce the security codes used for sending in results */
/* The security.h file is not required to compile a trial version of */
/* this program. */

#include "security.h"
#ifndef SEC1
#define SEC1(p)			0
#define SEC2(p,hi,lo,u,e)	0
#define SEC3(p)			0
#define SEC4(p)			0
#endif

/* Common routines */

int isKnownMersennePrime (unsigned long);
unsigned int avail_mem ();
void clear_memory (unsigned long, unsigned long);
void makestr (unsigned long, unsigned long, unsigned long, char *buf);
void updateWorkToDo (unsigned long, int, unsigned long);

int stopCheck ();
#define STOP_ESCAPE		1
#define STOP_TIMEOUT		2
#define STOP_MEM_CHANGED	3
#define STOP_NOT_ENOUGH_MEM	4
extern int STOP_REASON;		/* Reason stopCheck stopped processing */
void memSettingsChanged ();

void primeContinue ();
int prime (unsigned long, unsigned long, unsigned long, unsigned long);
int selfTest (unsigned long);
int selfTestInternal (unsigned long, int);
void primeTime (unsigned long, unsigned long);
int primeFactor (unsigned long, unsigned int, int *, int, int);
int primeSieve (unsigned long, unsigned long, unsigned short,
		unsigned short, int);
void primeSieveTest ();
int pfactor (struct work_unit *);

/* Routines used to time code chunks */
extern double __timers[10];		/* 10 timers are available */

/* Routines called by common routines */

void title (char *);
void flashWindowAndBeep ();
void SetPriority ();
int ecm (unsigned long, unsigned long, unsigned long, unsigned long,
	 unsigned long, unsigned long, double, int);
int pminus1 (unsigned long, unsigned long, unsigned long, unsigned long,
	     int, int);
