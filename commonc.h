/* Constants */

#define VERSION		"23.8"
#define VERSION_BIT	13      	/* Bit number in broadcast map */
/* The list of assigned version bits follows: */
/* Version 23.0    uses bit #13 */
/* Version 22.0    uses bit #12 */
/* Version 21.0    uses bit #11 */
/* Version 20.0    uses bit #10 */
/* Version 19.1    uses bit #9 */
/* Version 19.0    uses bit #8 */
/* The list of assigned OS ports follows: */
/* Win9x (prime95) uses bit #1 */
/* Linux (mprime)  uses bit #2 */
/* Solaris	   uses bit #3 */
/* WinNT (ntprime) uses bit #5 */
/* FreeBSD(mprime) uses bit #6 */
/* OS/2		   uses bit #7 */

#define MIN_PRIME	5L		/* Smallest testable prime */
#define MAX_FACTOR	MAX_PRIME	/* Largest factorable number */
#define ERROR_RATE	0.018		/* Estimated error rate on clean run */

/* Factoring limits based on complex formulas given the speed of the */
/* factoring code vs. the speed of the Lucas-Lehmer code */

#define FAC72	71000000L
#define FAC71	57020000L
#define FAC70	44150000L
#define FAC69	35200000L
#define FAC68	28130000L
#define FAC67	21590000L
#define FAC66	17890000L
#define FAC65	13380000L
#define FAC64	8250000L
#define FAC63	6515000L
#define FAC62	5160000L
#define FAC61	3960000L
#define FAC60	2950000L
#define FAC59	2360000L
#define FAC58	1930000L
#define FAC57	1480000L
#define FAC56	1000000L

/* Global variables */

extern char INI_FILE[80];		/* Name of the prime INI file */
extern char LOCALINI_FILE[80];		/* Name of the local INI file */
extern char WORKTODO_FILE[80];		/* Name of the work-to-do INI file */
extern char RESFILE[80];		/* Name of the results file */
extern char SPOOL_FILE[80];		/* Name of the spool file */
extern char LOGFILE[80];		/* Name of the server log file */
extern char EXTENSION[8];		/* Extension for several filenames */

extern char USERID[20];			/* User's ID */
extern char USER_PWD[20];		/* User's Password */
extern char OLD_USERID[20];		/* Previous User's ID */
extern char OLD_USER_PWD[20];		/* Previous User's Password */
extern char COMPID[20];			/* Computer ID */
extern char USER_NAME[80];		/* User's real name */
extern char USER_ADDR[80];		/* User's email address */
extern int NEWSLETTERS;			/* Send email about newsletters */
extern int USE_PRIMENET;		/* TRUE if we're using PrimeNet */
extern int DIAL_UP;			/* TRUE if we're dialing into */
					/* PrimeNet server */
extern short WORK_PREFERENCE;		/* Type of work (factoring, testing, */
					/* etc.) to get from the server. */
extern unsigned int DAYS_OF_WORK;	/* How much work to retrieve from */
					/* the primenet server */
extern time_t VACATION_END;		/* Date vacation is expected to end */
extern int ON_DURING_VACATION;		/* TRUE if computer is on while */
					/* vacationing */
extern int ADVANCED_ENABLED;		/* 1 if advanced menu is enabled */
extern int volatile ERRCHK;		/* 1 to turn on error checking */
extern unsigned int PRIORITY;		/* Desired priority level */
extern unsigned int CPU_AFFINITY;	/* NT Processor affinity */
extern int MANUAL_COMM;			/* Set on if user explicitly starts */
					/* all communication with the server */
EXTERNC unsigned int volatile CPU_TYPE;	/* 3=Cyrix, 4=486, 5=Pentium, */
					/* 6=Pro, 7=K6, 8=Celeron, 9=P-II */
					/* 10=P-III, 11=K7, 12=P4 */
extern unsigned int volatile CPU_HOURS;	/* Hours per day program will run */
extern unsigned int volatile DAY_MEMORY;/* Mem available in megabytes */
extern unsigned int volatile NIGHT_MEMORY;/* Mem available in megabytes */
extern unsigned int volatile DAY_START_TIME;/* When mem is first avail */
extern unsigned int volatile DAY_END_TIME;/* When mem is no longer avail */
extern unsigned long volatile ITER_OUTPUT;/* Iterations between outputs */
extern unsigned long volatile ITER_OUTPUT_RES;/* Iterations between results */
					/* file outputs */
extern unsigned long volatile DISK_WRITE_TIME;
					/* Number of minutes between writing */
					/* intermediate results to disk */
extern unsigned int MODEM_RETRY_TIME;	/* How often to try sending msgs */
					/* to primenet server whem modem off */
extern unsigned int NETWORK_RETRY_TIME;	/* How often to try sending msgs */
					/* to primenet server */
extern unsigned int DAYS_BETWEEN_CHECKINS; /* Days between sending updated */
					/* completion dates to the server */
extern int TWO_BACKUP_FILES;		/* TRUE for 2 backup files(qXXXXXXX) */
extern int SILENT_VICTORY;		/* Quiet find of new Mersenne prime */
extern int RUN_ON_BATTERY;		/* Run program even on battery power */
extern int TRAY_ICON;			/* Display tiny tray icon */
extern int HIDE_ICON;			/* Display no icon */
extern unsigned int ROLLING_AVERAGE;	/* Ratio of this computer's speed */
					/* compared to the expected speed */
					/* for this CPU */
extern unsigned int PRECISION;		/* Number of decimal places to output*/
					/* in percent complete lines */
extern time_t END_TIME;			/* Time at which we should reread */
					/* the INI files for new settings */
extern time_t SLEEP_TIME;		/* Time at which we should wake up */
					/* from sleep and read INI files */
extern int RDTSC_TIMING;		/* True if RDTSC is used to time */
extern int TIMESTAMPING;		/* True is timestamps to be output */
extern int CUMULATIVE_TIMING;		/* True if outputting cumulative time*/
extern int WELL_BEHAVED_WORK;		/* TRUE if undocumented feature */
					/* "well behaved worktodo file" */
					/* is on.  This reduces the number */
					/* of times worktodo.ini is read */
					/* and written. */
extern char **PAUSE_WHILE_RUNNING;	/* An array of program names that, */
					/* if running, prime95 should pause. */

extern unsigned long EXP_BEING_WORKED_ON; /* Exponent being tested */
extern int EXP_BEING_FACTORED;		/* TRUE is exp is being factored */
extern double EXP_PERCENT_COMPLETE;	/* Percent complete of factoring */
					/* and/or LL test */

extern int SPOOL_FILE_CHANGED;		/* Flag indicating the data in the */
					/* spool file has been updated. */
extern int CHECK_WORK_QUEUE;		/* Flag indicating we need to check */
					/* if our work queue is full. */
extern int GIMPS_QUIT;			/* TRUE if we just successfully */
					/* quit the GIMPS project */
extern time_t next_comm_time;		/* Next time to contact server */

extern char READFILEERR[];

/* Common routines */

void getCpuInfo (void);
void getCpuDescription (char *, int);

int isPrime (unsigned long p);
unsigned int max_mem (void);
unsigned int strToMinutes (char	*);
void minutesToStr (unsigned int, char *);

void nameIniFiles (int named_ini_files);
void readIniFiles (void);

void IniGetString (char *, char *, char *, unsigned int, char *);
long IniGetInt (char *, char *, long);
void IniWriteString (char *, char *, char *);
void IniWriteInt (char *, char *, long);

void IniFileOpen (char *, int);
void processTimedIniFile (char *);
void IniFileClose (char *);
unsigned int IniGetNumLines (char *);
void IniGetLineAsString (char *, unsigned int, char *, unsigned int,
			 char *, unsigned int);
void IniGetLineAsInt (char *, unsigned int, char *, unsigned int, long *);
void IniReplaceLineAsString (char *, unsigned int, char *, char *);
void IniReplaceLineAsInt (char *, unsigned int, char *, long);
void IniInsertLineAsString (char *, unsigned int, char *, char *);
void IniInsertLineAsInt (char *, unsigned int, char *, long);
void IniAppendLineAsString (char *, char *, char *);
void IniAppendLineAsInt (char *, char *, long);
void IniDeleteLine (char *, unsigned int);
void IniDeleteAllLines (char *);

void UpdateEndDates (void);
void ConditionallyUpdateEndDates (void);
void spoolMessage (short, void *);
void readMessage (int, long *, short *, void *);
int sendMessage (short, void *);
void spoolExistingResultsFile (void);
void OutputBoth (char *);
void OutputSomewhere (char *);
void LogMsg (char *);

struct work_unit {
	int	work_type;	/* Type of work to do */
	unsigned long p;	/* Exponent to work on */
	unsigned int bits;	/* How far factored */
	int	pminus1ed;	/* TRUE if has been P-1 factored */
	unsigned long B1;	/* ECM and P-1 - Stage 1 bound */
	unsigned long B2_start;	/* ECM and P-1 - Stage #2 start */
	unsigned long B2_end;	/* ECM and P-1 - Stage #2 end */
	unsigned int curves_to_do; /* ECM - curves to try */
	unsigned int curves_completed; /* ECM - curves done */
	double	curve;		/* ECM - Specific curve to test */
	int	plus1;		/* Flag for factoring 2^p+1 */
};
short default_work_type (void);
#define WORK_FACTOR		0
#define WORK_TEST		1
#define WORK_ADVANCEDTEST	2
#define WORK_DBLCHK		3
#define WORK_ECM		4
#define WORK_PMINUS1		5
#define WORK_PFACTOR		6
#define WORK_ADVANCEDFACTOR	7
int parseWorkToDoLine (unsigned int, struct work_unit *);
void addWorkToDoLine (struct work_unit *);
void checkResultsFile (unsigned long, int *, int *);
void getWorkFromDatabase (unsigned long, unsigned long, int, int);

unsigned long secondsUntilVacationEnds (void);
double pct_complete (int, unsigned long, unsigned long *);
unsigned long fftlen_from_ini_file (unsigned long);
unsigned long advanced_map_exponent_to_fftlen (unsigned long);
double raw_work_estimate (struct work_unit *);
double work_estimate (struct work_unit *);
unsigned int factorLimit (unsigned long, int);
void guess_pminus1_bounds (unsigned long, unsigned int, int, unsigned long *,
			   unsigned long *, unsigned long *, double *);

void strupper (char *);
void tempFileName (char	*, unsigned long);
int fileExists (char *);
int readFileHeader (char *, int *, short *, unsigned long *);
int writeResults (char	*);

int communicateWithServer (void);
void unreserve (unsigned long);

/* Routines called by common routines */

int LoadPrimeNet (void);
void UnloadPrimeNet (void);
int PRIMENET (short, void *);
int isHighResTimerAvailable (void);
double getHighResTimer (void);
double getHighResTimerFrequency (void);
void OutputStr (char *);
unsigned long physical_memory (void);
unsigned long num_cpus (void);
int getDefaultTimeFormat (void);
void doMiscTasks (void);
int escapeCheck (void);
void BroadcastMessage (char *);
#define	WORKING_ICON	0
#define	IDLE_ICON	1
void ChangeIcon (int);
void BlinkIcon (int);
