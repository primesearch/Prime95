/* These macros produce the security codes used for sending in results */
/* The security.h file is not required to compile a trial version of */
/* this program. */

#include "security.h"
#ifndef SEC1
#define SEC1(p)			0
#define SEC2(p,hi,lo,u,e)	0
#define SEC3(p)			0
#define SEC4(p)			0
#define SEC5(p,b1,b2)		0
#endif

/* Common routines */

int isKnownMersennePrime (unsigned long);
unsigned int avail_mem (void);
void clear_memory (unsigned long, unsigned long);
void makestr (unsigned long, unsigned long, unsigned long, char *buf);
void updateWorkToDo (unsigned long, int, unsigned long);

int stopCheck (void);
#define STOP_ESCAPE		1
#define STOP_TIMEOUT		2
#define STOP_MEM_CHANGED	3
#define STOP_NOT_ENOUGH_MEM	4
extern int STOP_REASON;		/* Reason stopCheck stopped processing */
void memSettingsChanged (void);

void primeContinue (void);
int prime (unsigned long, unsigned long, unsigned long, unsigned long);
int selfTest (unsigned long);
int selfTestInternal (unsigned long, unsigned int, int);
void primeTime (unsigned long, unsigned long);
void primeBench (void);
int primeFactor (unsigned long, unsigned int, int *, int, int);
int primeSieve (unsigned long, unsigned long, unsigned short,
		unsigned short, int);
void primeSieveTest (void);
int pfactor (struct work_unit *);

/* Routines used to time code chunks */
extern double __timers[10];		/* 10 timers are available */

/* Routines called by common routines */

void title (char *);
void flashWindowAndBeep (void);
void SetPriority (void);
int ecm (unsigned long, unsigned long, unsigned long, unsigned long,
	 unsigned long, unsigned long, double, int);
int pminus1 (unsigned long, unsigned long, unsigned long, unsigned long,
	     int, int);

/* Messages */

#define BENCH_SPEED  "The CPU speed in Options/CPU may not be correct.  An incorrect value will result in inaccurate timings.  Are you sure this is the correct CPU speed value?"
