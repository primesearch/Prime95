//#define SERVER_TESTING

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

int LaunchWorkerThreads (int);
int LaunchTortureTest (unsigned long, int);
int LaunchBench (void);
int LaunchAdvancedTime (unsigned long, unsigned long);

/* Callback routines */

#define LD_CONTINUE	0		/* Call primeContinue */
#define LD_TIME		1		/* Call primeTime */
#define LD_BENCH	2		/* Call primeBench */
#define LD_TORTURE	3		/* Call selfTest */
void PreLaunchCallback (int launch_type);	/* Implemented for each OS */
void PostLaunchCallback (int launch_type);	/* Implemented for each OS */
void stopCheckCallback (int thread_num);
int checkPauseListCallback (void);

/* Routines and structures to set thread priority.  Since we may have */
/* different rules doing normal work vs. benchmarking vs. torture testing, */
/* we send that info to the SetPriority routine to use as it sees fit. */

#define SET_PRIORITY_NORMAL_WORK	1
#define SET_PRIORITY_BENCHMARKING	2
#define SET_PRIORITY_TORTURE		3
#define SET_PRIORITY_QA			4
struct PriorityInfo {
 	int	type;		/* Type defined above */
	int	thread_num;	/* Worker thread number */
	int	aux_thread_num;	/* Set when gwnum launches auxillary threads */
	int	num_threads;	/* Total number of torture test threads */
};
void SetPriority (struct PriorityInfo *);

/* Internal routines that do the real work */

void create_worker_windows (int num_threads);
void Launcher (void *arg);
void LauncherDispatch (void *arg);
int primeContinue (int);
int tortureTest (int, int);
int selfTest (int, struct PriorityInfo *, struct work_unit *);
int selfTestInternal (int, struct PriorityInfo *, unsigned long, unsigned int,
		      int *, unsigned int, void *, int *, int *);
int primeTime (int, unsigned long, unsigned long);
int primeBench (int);
int primeFactor (int, struct PriorityInfo *, struct work_unit *, unsigned int);
int prime (int, struct PriorityInfo *, struct work_unit *, int);
int prp (int, struct PriorityInfo *, struct work_unit *, int);
int ecm (int, struct PriorityInfo *, struct work_unit *, int);
int pminus1 (int, struct PriorityInfo *, struct work_unit *, int);
int pfactor (int, struct PriorityInfo *, struct work_unit *, int);

/* Utility routines */

int isKnownMersennePrime (unsigned long);
void makestr (unsigned long, unsigned long, unsigned long, char *);

/* Stop routines */

/* Reasons returned for why we stopped processing a work unit */
#define STOP_ESCAPE		1	/* User hit escape key */
#define STOP_OUT_OF_MEM		2	/* Fatal out-of-memory error */
#define STOP_FILE_IO_ERROR	3	/* Fatal file I/O error */
#define STOP_FATAL_ERROR	4	/* Other fatal error */
#define STOP_ABORT		5	/* Abort current work unit */
#define STOP_WORK_UNIT_COMPLETE	50	/* Work unit is done! */
#define STOP_PRIORITY_WORK	51	/* Priority work, restart thread */
#define STOP_BATTERY		52	/* On battery - pause */
#define STOP_TIMEOUT		100	/* Time= prime.ini change */
#define STOP_RESTART		101	/* Important INI option changed */
#define STOP_MEM_CHANGED	102	/* Day/night memory change */
#define STOP_NOT_ENOUGH_MEM	103	/* Not enough memory for P-1 stage 2 */
#define STOP_NOT_DESIRED_MEM	104	/* Rather not share mem P-1 stage 2 */

EXTERNC int stopCheck (int);
void stop_workers_for_escape (void);
void stop_workers_for_restart (void);
void stop_workers_for_add_files (void);
void restart_waiting_workers (void);
void stop_worker_for_abort (int);

/* Routines dealing with day/night memory settings */

void start_mem_changed_timer (void);
void mem_settings_have_changed (void);
unsigned int max_mem (void);
unsigned int avail_mem (int thread_num);

/* battery routines */

void run_on_battery_changed (void);
int OnBattery (void);			/* Implemented for each OS */
void test_battery (void);

/* priority work routines */

void check_for_priority_work (void);
void stop_worker_for_advanced_test (int);

/* "pause while running" routines */

void checkPauseWhileRunning (void);
void implementPause (int thread_num);

/* throttle routines */

void start_throttle_timer (void);
void stop_throttle_timer (void);
int handleThrottleTimerEvent (void);
void implementThrottle (int thread_num);

/* Routines called by common routines */

void clearThreadHandleArray (void);
void setThreadPriorityAndAffinity (int, int);
void registerThreadTermination (void);
void raiseAllWorkerThreadPriority (void);
void flashWindowAndBeep (void);
int primeSieveTest (int);
int setN (gwhandle *, int, struct work_unit *, giant *);
int ecm_QA (int, struct PriorityInfo *);
int pminus1_QA (int, struct PriorityInfo *);
int test_randomly (int, struct PriorityInfo *);
int test_all_impl (int, struct PriorityInfo *);

/* Messages */

#define BENCH_SPEED  "The CPU speed in Options/CPU may not be correct.  An incorrect value will result in inaccurate timings.  Are you sure this is the correct CPU speed value?"
