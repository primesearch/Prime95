/* Copyright 1995-2000 Just For Fun Software, Inc. */
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
#include <dirent.h>
#include <fcntl.h>
#include <math.h>
#include <memory.h>
#include <signal.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <sys/stat.h>
#include <unistd.h>
#if defined (__linux__) || defined (__FreeBSD__)
#include <sys/time.h>
#include <sys/timeb.h>
#endif
#include <sys/resource.h>

/* Globals */

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

#ifdef __linux__
#define PORT	2
#endif
#ifdef __FreeBSD__
#define PORT	6
#endif
#ifdef __EMX__
#define PORT	7
#endif

#include "giants.c"
#include "gwnum.c"
#include "commona.c"
#include "commonb.c"
#include "commonc.c"
#include "ecm.c"

/* Signal handlers */

void sigterm_handler(int signo)
{
	THREAD_STOP = TRUE;
	if (signo != SIGINT) THREAD_KILL = TRUE;
	(void)signal(signo, sigterm_handler);
}

#ifdef MPRIME_LOADAVG

/* Routine to get the current load average */
double get_load_average ()
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
void init_load_check ()
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
			MENUING = TRUE;
			NO_GUI = FALSE;
			break;

/* -V - version number */

		case 'V':
		case 'v':
			printf ("Mersenne Prime Test Program, Version %s.2\n", VERSION);
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
				exit (0);
			}
			for (fd = 0; fd < OPEN_MAX; fd++) {
				close (fd);
			}
			open ("/dev/null", O_APPEND);
			dup2 (0,1);
			dup2 (0,2);
			setsid ();

			if (i) named_ini_files = i;
		}
	}

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

/* If we are to contact the server, do so now.  This option lets the */
/* user create a batch file that contacts the server at regular intervals */
/* or when the ISP is contacted, etc. */

	if (contact_server) {
		MANUAL_COMM = 3;
		CHECK_WORK_QUEUE = 1;
		communicateWithServer ();
	}

/* Bring up the main menu */

	else if (MENUING)
		main_menu ();

/* On first run, get user name and email address */
/* before contacting server for a work assignment */

	else if (USE_PRIMENET && USERID[0] == 0) {
		STARTUP_IN_PROGRESS = 1;
		test_user ();
	}

/* Continue testing the range */

	else
		linuxContinue ("Another mprime is already running!\n");

/* All done */

	return (0);

/* Invalid args message */

usage:	printf ("Usage: mprime [-aN] [-bcdhmv] [-wDIR]\n");
	printf ("-aN\tUse an alternate set of INI and output files.\n");
	printf ("-bN\tRun in the background.\n");
	printf ("-c\tContact the PrimeNet server, then exit.\n");
	printf ("-d\tPrint detailed information to stdout.\n");
	printf ("-h\tPrint this.\n");
	printf ("-m\tMenu to configure mprime.\n");
	printf ("-v\tPrint the version number.\n");
	printf ("-wDIR\tRun from a different working directory.\n");
	printf ("\n");
	return (1);
}

void title (char *msg)
{
}

void flashWindowAndBeep ()
{
	printf ("\007");
}

/* Return TRUE if we should stop calculating */

int escapeCheck ()
{
	if (THREAD_STOP) {
		THREAD_STOP = 0;
		return (TRUE);
	}
	return (FALSE);
}

void doMiscTasks ()
{
#ifdef MPRIME_LOADAVG
	test_sleep ();
#endif
}

void OutputStr (char *buf)
{
	if (VERBOSE || MENUING) printf ("%s", buf);
}

void guessCpuType ()
{
	FILE	*fd;
	char	buf[80];

	CPU_TYPE = isPentiumPro () ? 6 : isPentium () ? 5 : 4;
	CPU_SPEED = 100;
	fd = fopen ("/proc/cpuinfo", "r");
	if (fd == NULL) return;
	for ( ; ; ) {
		double	speed;
		if (fscanf (fd, "%s", buf) == EOF) break;
		if (strcmp (buf, "MHz") == 0) {
			fscanf (fd, " : %lf", &speed);
			if (speed > 25.0 && speed < 10000.0)
				CPU_SPEED = (unsigned long) (speed + 0.5);
			break;
		}
	}
	fclose (fd);
}

unsigned long physical_memory ()
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

int getDefaultTimeFormat ()
{
	return (2);
}

void Sleep (
	long	ms) 
{
	sleep (ms/1000);
}

/* Set priority.  Map one (prime95's lowest priority) to 20 */
/* (linux's lowest priority).  Map eight (prime95's normal priority) to */
/* 0 (linux's normal priority). */

void SetPriority ()
{
	int	p;
	p = (8 - (int) PRIORITY) * 20 / 7;
	setpriority (PRIO_PROCESS, getpid (), p);
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

/* The two pids are running the same executable, raise an error and return */

	if (error_message != NULL) printf ("%s", error_message);
	return;

/* All is OK.  Save our pid, run, then delete our pid */

ok:	IniWriteInt (LOCALINI_FILE, "Pid", my_pid);
	primeContinue ();
	IniWriteInt (LOCALINI_FILE, "Pid", 0);
}
