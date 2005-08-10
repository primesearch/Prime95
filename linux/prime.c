/* Copyright 1995-2005 Just For Fun Software, Inc. */
/* Author:  George Woltman */
/* Email: woltman@alum.mit.edu */

/* Include files */

#include "prime.h"
#ifdef __FreeBSD__
/* FreeBSD needs to process sys/types.h before it can understand either
/* sys/time.h or sys/resource.h */
#include <sys/types.h>
#endif
#include <ctype.h>
#ifndef __IBMC__
#include <dirent.h>
#endif
#include <fcntl.h>
#include <math.h>
#include <memory.h>
#include <signal.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <sys/stat.h>
#ifndef __IBMC__
#include <unistd.h>
#endif
#if defined (__linux__) || defined (__FreeBSD__) || defined(__WATCOMC__)
#include <sys/time.h>
#include <sys/timeb.h>
#endif
#ifndef __IBMC__
#include <sys/resource.h>
#endif
/* Required OS/2 header files */
#ifdef __IBMC__
#define INCL_DOS
#define INCL_DOSPROFILE
#include <os2.h>
#include <direct.h>
#include <io.h>
#include <process.h>
#include <sys/timeb.h>
typedef int pid_t;
#include "dosqss.h"
#endif

/* Globals */

#ifndef __WATCOMC__
#define OPEN_MAX 20
#endif

#ifdef MPRIME_LOADAVG
#define LINUX_LDAV_FILE "/proc/loadavg"
int volatile SLEEP_STOP = 0;
long LOAD_CHECK_TIME = 0;
double HI_LOAD = 0.0;
double LO_LOAD = 0.0;
#endif

int volatile THREAD_STOP = 0;
int volatile THREAD_KILL = 0;
int NO_GUI = 1;
int VERBOSE = 0;
int MENUING = 0;

/* Common code */

#include "gwutil.h"
#include "commona.c"
#include "commonb.c"
#include "commonc.c"
#include "ecm.c"
#include "primenet.c"

/* Signal handlers */

void sigterm_handler(int signo)
{
	THREAD_STOP = TRUE;
	if (signo != SIGINT) THREAD_KILL = TRUE;
	(void)signal(signo, sigterm_handler);
}

#ifdef MPRIME_LOADAVG

/* Routine to get the current load average */
double get_load_average (void)
{
#ifdef __linux__
	char	ldavgbuf[40];
	double	load_avg;
	int	fd, count;

	fd = open (LINUX_LDAV_FILE, O_RDONLY);
	if (fd == -1) return (-1.0);
	count = read (fd, ldavgbuf, 40);
	(void) close (fd);
	if (count <= 0) return (-1.0);
	count = sscanf (ldavgbuf, "%lf", &load_avg);
	if (count < 1) return (-1.0);
	return (load_avg);
#endif
#ifdef __FreeBSD__
	double load[3];

	if (getloadavg (load, sizeof(load)/sizeof(load[0])) < 0) return (-1.0);
	return (load[0]);
#endif
}

/* load_handler: call by signal routine,
   sets SLEEP_STOP to TRUE if load is too high */
void load_handler (
	int	sig)
{
	double  load_avg;

	load_avg = get_load_average ();
	if (load_avg < 0.0) return;
  
	if (SLEEP_STOP) {
		if (load_avg < LO_LOAD)
			SLEEP_STOP = FALSE;
	} else {
		if (load_avg > HI_LOAD)
			SLEEP_STOP = TRUE;
	}
}

/* init_load_check: initialises timer that calls load_handler
   every LOAD_CHECK_TIME seconds */
void init_load_check (void)
{
	struct itimerval timer, otimer;
	struct sigaction sigact;
	int	ret;

	timer.it_interval.tv_sec  =  LOAD_CHECK_TIME;
	timer.it_interval.tv_usec =  0;
	timer.it_value.tv_sec     =  LOAD_CHECK_TIME;
	timer.it_value.tv_usec    =  0;

	ret = setitimer (ITIMER_REAL, &timer, &otimer);
	if (ret < 0) return;
  
	sigact.sa_handler = &load_handler;
	sigemptyset(&sigact.sa_mask);
	sigact.sa_flags =  SA_RESTART;
	ret = sigaction(SIGALRM, &sigact, NULL);
	if (ret < 0) { /* clean up after ourselves */
		setitimer (ITIMER_REAL, &otimer, NULL);
	}
}

/* test_sleep: tests if SLEEP_STOP is set and sleeps until load is normal
   again or THREAD_STOP is set
*/
void test_sleep (void) 
{
	sigset_t newmask;

	while (SLEEP_STOP && !THREAD_STOP) {
		sigemptyset (&newmask);
		sigsuspend (&newmask);
	}
}
#endif

/* Main entry point! */

int main (
	int	argc,
	char	*argv[])
{
	char	buf[256];
	int	named_ini_files = -1;
	int	background = 0;
	int	contact_server = 0;
	int	torture_test = 0;
	int	i;
	char	*p;

/* catch termination signals */

	(void)signal(SIGTERM, sigterm_handler);
	(void)signal(SIGINT, sigterm_handler);

/* No buffering of output */

	setvbuf (stdout, NULL, _IONBF, 0);

/* Change to the executable's directory*/

	strcpy (buf, argv[0]);
	p = strrchr (buf, '/');
	if (p != NULL) {
		*p = 0;
		_chdir (buf);
	}

/* Initialize gwnum call back routines.  Using callback routines lets the */
/* gwnum library have a nice clean interface for users that do not need */
/* additional functionality that only prime95 uses. */

	StopCheckRoutine = stopCheck;
	OutputBothRoutine = OutputBoth;

/* Process command line switches */

	for (i = 1; i < argc; i++) {
		p = argv[i];

		if (*p++ != '-') break;
		switch (*p++) {

/* Accept a -A switch indicating an alternate set of INI files */
/* are to be used. */

		case 'A':
		case 'a':
			named_ini_files = 0;
			while (isspace (*p)) p++;
			while (isdigit (*p)) {
				named_ini_files = named_ini_files * 10 + (*p - '0');
				p++;
			}
			break;

/* -B - put in the background.  Accepts an optional CPU count. */

		case 'B':
		case 'b':
			while (isspace (*p)) p++;
			if (isdigit (*p)) {
				background = 0;
				while (isdigit (*p)) {
					background = background * 10 + (*p - '0');
					p++;
				}
			} else
				background = 1;
			break;

/* -C - contact the server now, then exit */

		case 'C':
		case 'c':
			contact_server = 1;
			VERBOSE = TRUE;
			NO_GUI = FALSE;
			break;
			
/* -D - debug */

		case 'D':
		case 'd':
			VERBOSE = TRUE;
			NO_GUI = FALSE;
			break;

/* -H - help */

		case 'H':
		case 'h':
		case '?':
			goto usage;

/* -M - Menu */

		case 'M':
		case 'm':
			MENUING = 1;
			NO_GUI = FALSE;
			break;

/* -S - status */

		case 'S':
		case 's':
			MENUING = 2;
			NO_GUI = FALSE;
			break;
		  

/* -T - Torture test */

		case 'T':
		case 't':
			torture_test = TRUE;
			break;

/* -V - version number */

		case 'V':
		case 'v':
			printf ("Mersenne Prime Test Program, Version %s.%d\n", VERSION, PORT);
			return (0); 

/* -W - use a different working directory */

		case 'W':
		case 'w':
			_chdir (p);
			break; 

/* Otherwise unknown switch */

		default:
			printf ("Invalid switch\n");
			goto usage;
		}
	}

/* Ignore background request on OS/2 */

#ifndef __IBMC__

/* Run in background if requested.  Code courtesy of Francois Gouget. */

	if (background) {
		int	i;

/* To enter daemon mode, close all the filedescs and detach from the tty */

		for (i = 0; i < background; i++) {
			int	fd;
			fd = fork ();
			if (fd == -1) {
				perror ("Could not fork to the background");
				exit (1);
			}
			if (fd != 0) {
				if (i == background - 1) exit (0);
				else continue;
			}
			for (fd = 0; fd < OPEN_MAX; fd++) {
				close (fd);
			}
			open ("/dev/null", O_APPEND);
			dup2 (0,1);
			dup2 (0,2);
			setsid ();

			if (i) named_ini_files = i;
			break;
		}
	}
#endif

/* Determine the names of the INI files */
/* Read the INI files */

	nameIniFiles (named_ini_files);
	readIniFiles ();

/* Read load averaging settings from INI files */

#ifdef MPRIME_LOADAVG
	IniGetString (INI_FILE, "MaxLoad", buf, sizeof (buf), "0");
	HI_LOAD = atof (buf);
	IniGetString (INI_FILE, "MinLoad", buf, sizeof (buf), "0");
	LO_LOAD = atof (buf);
	IniGetString (INI_FILE, "PauseTime", buf, sizeof (buf), "0");
	LOAD_CHECK_TIME = atol (buf);

/* Initialise load checking */

	if (HI_LOAD > 0.0 && LOAD_CHECK_TIME > 0)
		init_load_check ();
#endif

/* If a broadcast message from the server has been received but */
/* never viewed by the user, then try to display it now. */

	BroadcastMessage (NULL);

/* If this is a stress tester, then turn on menuing.  A stress tester */
/* ceases to be a stress tester if he ever turns on primenet or has work */
/* in his worktodo.ini file */

	if (IniGetInt (INI_FILE, "StressTester", 0)) {
		if (USE_PRIMENET || IniGetNumLines (WORKTODO_FILE)) {
			IniWriteInt (INI_FILE, "StressTester", 0);
		} else {
			MENUING = 1;
			VERBOSE = TRUE;
			NO_GUI = FALSE;
		}
	}

/* If running the torture test, do so now. */

	if (torture_test) {
		VERBOSE = TRUE;
		NO_GUI = FALSE;
		selfTest (1);
	}

/* On first run, get user name and email address before contacting server */
/* for a work assignment.  To make first time user more comfortable, we will */
/* display data to the screen, rather than running silently. */

	else if (USE_PRIMENET &&
		 USERID[0] == 0 &&
		 !IniGetInt (INI_FILE, "StressTester", 0)) {
		VERBOSE = TRUE;
		NO_GUI = FALSE;
		STARTUP_IN_PROGRESS = 1;
		test_welcome ();
	}

/* If we are to contact the server, do so now.  This option lets the */
/* user create a batch file that contacts the server at regular intervals */
/* or when the ISP is contacted, etc. */

	else if (contact_server) {
		MANUAL_COMM = 3;
		CHECK_WORK_QUEUE = 1;
		communicateWithServer ();
	}

/* Bring up the main menu */

	else if (MENUING == 1)
		main_menu ();
	else if (MENUING == 2)
		rangeStatus();
	
/* Continue testing */

	else
		linuxContinue ("Another mprime is already running!\n");

/* All done */

	return (0);

/* Invalid args message */

usage:	printf ("Usage: mprime [-cdhmstv] [-aN] [-b[N]] [-wDIR]\n");
	printf ("-c\tContact the PrimeNet server, then exit.\n");
	printf ("-d\tPrint detailed information to stdout.\n");
	printf ("-h\tPrint this.\n");
	printf ("-m\tMenu to configure mprime.\n");
	printf ("-s\tDisplay status.\n");
	printf ("-t\tRun the torture test.\n");
	printf ("-v\tPrint the version number.\n");
	printf ("-aN\tUse an alternate set of INI and output files.\n");
	printf ("-bN\tRun in the background.  N is number of CPUs.\n");
	printf ("-wDIR\tRun from a different working directory.\n");
	printf ("\n");
	return (1);
}

void title (char *msg)
{
}

void flashWindowAndBeep (void)
{
	printf ("\007");
}

/* Return TRUE if we should stop calculating */

int escapeCheck (void)
{
	if (THREAD_STOP) {
		THREAD_STOP = 0;
		return (TRUE);
	}
	return (FALSE);
}

void doMiscTasks (void)
{
#ifdef MPRIME_LOADAVG
	test_sleep ();
#endif
}

void OutputStr (char *buf)
{
	if (VERBOSE || MENUING) printf ("%s", buf);
}

unsigned long physical_memory (void)
{
	FILE	*fd;
	char	mem[80];

	fd = fopen ("/proc/meminfo", "r");
	if (fd == NULL) return (1024);
	for ( ; ; ) {
		if (fscanf (fd, "%s", mem) == EOF) {
			fclose (fd);
			return (1024);
		}
		if (isdigit (mem[0])) break;
	}
	fclose (fd);
	return (atoi (mem) >> 20);
}

/* Return a better guess for amount of memory to use in a torture test. */
/* Caller passes in its guess for amount of memory to use, but this routine */
/* can reduce that guess based on OS-specific code that looks at amount */
/* of available physical memory. */
/* This code was written by an anonymous GIMPS user. */

unsigned long GetSuggestedMemory (unsigned long nDesiredMemory)
{
	return (nDesiredMemory);
}

int getDefaultTimeFormat (void)
{
	return (2);
}

void Sleep (
	long	ms) 
{
#ifdef __IBMC__
	DosSleep(ms);
#else
	usleep (ms * 1000);
#endif
}

/* Set priority.  Map one (prime95's lowest priority) to 20 */
/* (linux's lowest priority).  Map eight (prime95's normal priority) to */
/* 0 (linux's normal priority). */

void SetPriority (void)
{
#ifdef __IBMC__
	DosSetPriority(PRTYS_PROCESS,
		(PRIORITY < 6) ? PRTYC_IDLETIME : PRTYC_REGULAR,
		(PRIORITY == 1 || PRIORITY == 6) ? PRTYD_MINIMUM :
		(PRIORITY == 2 || PRIORITY == 7) ? -10 :
		(PRIORITY == 3 || PRIORITY == 8) ? 0 :
		(PRIORITY == 4 || PRIORITY == 9) ? 10 :
		PRTYD_MAXIMUM,
		0);
#else
	int	p;
	p = (8 - (int) PRIORITY) * 20 / 7;
	setpriority (PRIO_PROCESS, getpid (), p);
#endif
}

void BlinkIcon (int x)
{
}

void ChangeIcon (int x)
{
}

void BroadcastMessage (
	char	*message)
{
	char	filename[33];
	int	fd, len;

/* Generate broadcast message file name */

        strcpy (filename, "bcastmsg");
        strcat (filename, EXTENSION);

/* If this is a call to check if a broadcast message exists, then do so */

	if (message == NULL) {
		if (! fileExists (filename)) return;
		fd = _open (filename, _O_TEXT | _O_RDONLY, 0);
		if (fd < 0) return;
		message = (char *) malloc (1024);
		len = _read (fd, message, 1024);
		_close (fd);
		if (len < 0) return;
		message[len] = 0;
		printf ("Important Message from PrimeNet Server:\n");
		printf ("%s", message);
		free (message);
		return;
	}

/* Otherwise, this is a new message - write it to the file */

	fd = _open (filename, _O_TEXT | _O_RDWR | _O_CREAT | _O_APPEND, 0666);
	if (fd < 0) return;
	_write (fd, message, strlen (message));
	_close (fd);
}


/* This routine calls primeContinue unless there is another copy of mprime */
/* already running.  In that case, it outputs an optional error message. */

void linuxContinue (
	char	*error_message)
{
#ifdef __linux__
#define PROCNAME	"/proc/%d/exe"
#endif
#ifdef __FreeBSD__
#define PROCNAME	"/proc/%d/file"
#endif
	pid_t	my_pid, running_pid;
	char	filename[30];
	int	fd;
	struct stat filedata;
	ino_t	inode1, inode2;

/* Compare this process' ID and the pid from the INI file */

	my_pid = getpid ();
	openIniFile (LOCALINI_FILE, 1);
	running_pid = IniGetInt (LOCALINI_FILE, "Pid", 0);
	if (running_pid == 0 || my_pid == running_pid) goto ok;

#ifdef __OS2__

        {
            USHORT handle1 = 0, handle2 = 0;
            unsigned char buf[0x2000];
            if( !DosQuerySysState(0x01, 0, 0, 0, (PCHAR)buf, 0x2000) ) {
                PQPROCESS p = ((PQTOPLEVEL)buf)->procdata;
                while(p && p->rectype == 1) {
                    if( p->pid == running_pid ) handle1 = p->hndmod;
                    if( p->pid == my_pid ) handle2 = p->hndmod;
                    p = (PQPROCESS)(p->threads + p->threadcnt);
                }
                if( handle1 != handle2 ) goto ok;
            }
        }

#else

/* See if the two pids are running the same executable */

	sprintf (filename, PROCNAME, my_pid);
	fd = _open (filename, _O_RDONLY);
	if (fd < 0) goto ok;
	fstat (fd, &filedata);
	inode1 = filedata.st_ino;
	_close (fd);
	sprintf (filename, PROCNAME, running_pid);
	fd = _open (filename, _O_RDONLY);
	if (fd < 0) goto ok;
	fstat (fd, &filedata);
	inode2 = filedata.st_ino;
	_close (fd);
	if (inode1 != inode2) goto ok;
#endif

/* The two pids are running the same executable, raise an error and return */

	if (error_message != NULL) printf ("%s", error_message);
	return;

/* All is OK, set affinity before running mprime. */

ok:
#ifdef __linux__
	{
		int	affinity_mask;
		affinity_mask = (CPU_AFFINITY == 99) ? -1 : (1 << CPU_AFFINITY);
		sched_setaffinity (0, sizeof (affinity_mask), &affinity_mask);
	}
#endif

/* Save our pid, run, then delete our pid */

	IniWriteInt (LOCALINI_FILE, "Pid", my_pid);
	primeContinue ();
	IniWriteInt (LOCALINI_FILE, "Pid", 0);
}

/* Load the PrimeNet DLL, make sure an internet connection is active */

int LoadPrimeNet (void)
{
	/* Init stuff */
	/* Set PRIMENET procedure pointer */
	/* return false if not connected to internet */

	int lines = 0;
#ifndef AOUT
	FILE* fd;
	char buffer[4096];
#ifdef __EMX__
	char command[128];
	char szProxyHost[120], *con_host;
	char *colon;
	
	IniGetString(INIFILENAME, "ProxyHost", szProxyHost, 120, NULL);
	if (*szProxyHost) {
		if ((colon = strchr(szProxyHost, ':'))) {
			*colon = 0;
		}
		con_host = szProxyHost;
	} else {
		con_host = szSITE;
	}

	sprintf(command,"host %s",con_host);
#ifdef __DEBUG
	fprintf(stderr,"Command = %s\n",command);
#endif
	fd = popen(command,"r");
	if (fd != NULL) {
	  fgets(buffer, 199, fd);
#ifdef __DEBUG
	  fprintf(stderr,"Response = %s\n",buffer);
#endif
	  if (strncmp(buffer,"host:",5) != 0) {
	    fclose(fd);
	    return TRUE;
	  }
	  fclose(fd);
	}
#else
#ifdef __linux__
	/* Open file that will hopefully tell us if we are connected to */
	/* the Internet.  There are four possible settings for RouteRequired */
	/* 0:	Always return TRUE */
	/* 1:   Use old version 19 code */
	/* 2:   Use new code supplied by Matthew Ashton. */
	/* 99:	Default.  Use case 2 above but if cannot open /proc/net/route*/
	/*	then assume you are connected (we probably do not have read */
	/*	permission or this is a funny Linux setup). */
	{
	  int RtReq = IniGetInt (INIFILENAME, "RouteRequired", 99);
	  if (RtReq == 0) return (TRUE);
	  fd = fopen("/proc/net/route","r");
	  if (fd == NULL) return (RtReq == 99);
	/* We have a readable /proc/net/route file.  Use the new check */
	/* for an Internet connection written by Matthew Ashton. However, */
	/* we still support the old style check (just in case) by setting */
	/* RouteRequired to 1. */
	  if (RtReq >= 2) {
	    while (fgets(buffer, sizeof(buffer), fd)) {
	      int dest;
	      if(sscanf(buffer, "%*s %x", &dest) == 1 && dest == 0) {
		fclose (fd);
		return (TRUE);
	      }
	    }
	  }
	/* The old code for testing an Internet connection is below */
	  else {
	    fgets(buffer, 199, fd);
	    fgets(buffer, 199, fd);
	    while (!feof(fd)) {
	      if (strncmp(buffer, "lo", 2)) {
	        fclose(fd);
	        return TRUE;
	      }
	      fgets(buffer, 199, fd);
	    }
	  }
	  fclose(fd);
	}
#endif
#if defined(__FreeBSD__) || defined(__WATCOMC__)
	/* The /proc/net/route test is not really meaningful under FreeBSD */
	/* There doesn't seem to be any meaningful test to see whether the */
	/* computer is connected to the Internet at the time using a non- */
	/* invasive test (which wouldn't, say, activate diald or ppp or */
	/* something else */
	return TRUE;
#endif                /* __FreeBSD__ */
#endif
#endif
	OutputStr ("You are not connected to the Internet.\n");
	return FALSE;
}

/* Unload the PrimeNet DLL */

void UnloadPrimeNet (void)
{
}

/* Check if a program is currently running - not implemented for OS/2 */

int checkPauseList ()
{
#ifndef __OS2__
	FILE	*fd;
	char	buf[80];

	fd = popen ("ps -eo comm", "r");
	if (fd != NULL) {
		while (fgets (buf, sizeof (buf), fd) != NULL) {
			int	len = strlen (buf);
			while (len && isspace (buf[len-1])) buf[--len] = 0;
			if (isInPauseList (buf)) {
				fclose (fd);
				return (TRUE);
			}
		}
		fclose (fd);
	}
#endif
	return (FALSE);
}
