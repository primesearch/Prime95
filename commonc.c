/*----------------------------------------------------------------------
| This file contains routines and global variables that are common for
| all operating systems the program has been ported to.  It is included
| in one of the source code files of each port.  See common.h for the
| common #defines and common routine definitions.
|
| Commona contains information used only during setup
| Commonb contains information used only during execution
| Commonc contains information used during setup and execution
+---------------------------------------------------------------------*/

char JUNK[]="Copyright 1996-2004 Just For Fun Software, All rights reserved";

char	INI_FILE[80] = {0};
char	LOCALINI_FILE[80] = {0};
char	WORKTODO_FILE[80] = {0};
char	RESFILE[80] = {0};
char	SPOOL_FILE[80] = {0};
char	LOGFILE[80] = {0};
char	EXTENSION[8] = {0};

char	USERID[20] = {0};
char	USER_PWD[20] = {0};
char	OLD_USERID[20] = {0};
char	OLD_USER_PWD[20] = {0};
char	COMPID[20] = {0};
char	USER_NAME[80] = {0};
char	USER_ADDR[80] = {0};
int	NEWSLETTERS = 1;
int	USE_PRIMENET = 1;
int	DIAL_UP = 0;
short	WORK_PREFERENCE = 0;
unsigned int DAYS_OF_WORK = 20;
time_t	VACATION_END = 0;
int	ON_DURING_VACATION = 1;
int	ADVANCED_ENABLED = 0;
int volatile ERRCHK = 0;
unsigned int PRIORITY = 1;
unsigned int CPU_AFFINITY = 99;
int	MANUAL_COMM = 0;
EXTERNC unsigned int volatile CPU_TYPE = 0;
unsigned int volatile CPU_HOURS = 0;
unsigned int volatile DAY_MEMORY = 0;
unsigned int volatile NIGHT_MEMORY = 0;
unsigned int volatile DAY_START_TIME = 0;
unsigned int volatile DAY_END_TIME = 0;
unsigned long volatile ITER_OUTPUT = 0;
unsigned long volatile ITER_OUTPUT_RES = 999999999;
unsigned long volatile DISK_WRITE_TIME = 30;
unsigned int MODEM_RETRY_TIME = 2;
unsigned int NETWORK_RETRY_TIME = 60;
unsigned int DAYS_BETWEEN_CHECKINS = 28;
int	TWO_BACKUP_FILES = 1;
int	SILENT_VICTORY = 0;
int	RUN_ON_BATTERY = 1;
int	TRAY_ICON = TRUE;
int	HIDE_ICON = FALSE;
unsigned int ROLLING_AVERAGE = 0;
unsigned int PRECISION = 2;
time_t	END_TIME = 0;
time_t	SLEEP_TIME = 0;
int	RDTSC_TIMING = 1;
int	TIMESTAMPING = 1;
int	CUMULATIVE_TIMING = 0;
int	WELL_BEHAVED_WORK = 0;
char	**PAUSE_WHILE_RUNNING = NULL;

unsigned long EXP_BEING_WORKED_ON = 0;
int	EXP_BEING_FACTORED = 0;
double	EXP_PERCENT_COMPLETE = 0.0;

int	SPOOL_FILE_CHANGED = 0;
int	CHECK_WORK_QUEUE = 0;
int	GIMPS_QUIT = 0;
time_t	next_comm_time = 0;	/* Next time we should contact server */

char READFILEERR[] = "Error reading intermediate file: %s\n";
char BLACKMSG1[] = "The PrimeNet Server has temporarily requested your computer not contact the server.\n";
char BLACKMSG2[] = "The PrimeNet Server has permanently requested your computer not contact the server.\nTry downloading a new version from http://www.mersenne.org/freesoft.htm or contact woltman@alum.mit.edu\n";

/* Determine the CPU speed either empirically or by user overrides. */
/* getCpuType must be called prior to calling this routine. */

void getCpuSpeed (void)
{
	int	temp, old_cpu_speed;

/* Guess the CPU speed using the RDTSC instruction */

	guessCpuSpeed ();

/* Now let the user override the cpu speed from the local.ini file */

	if (IniGetInt (LOCALINI_FILE, "CpuOverride", 0)) {
		temp = IniGetInt (LOCALINI_FILE, "CpuSpeed", 99);
		if (temp != 99) CPU_SPEED = temp;
	}

/* Make sure the cpu speed is reasonable */

	if (CPU_SPEED > 50000) CPU_SPEED = 50000;
	if (CPU_SPEED < 25) CPU_SPEED = 25;

/* If CPU speed has changed greatly since the last time, then */
/* tell the server, recalculate new completion dates, and reset the */
/* rolling average.  Don't do this on the first run (before the Welcome */
/* dialog has been displayed). */

	old_cpu_speed = IniGetInt (LOCALINI_FILE, "OldCpuSpeed", 0);
	if (CPU_SPEED < (double) old_cpu_speed - 5.0 ||
	    CPU_SPEED > (double) old_cpu_speed + 5.0) {
		IniWriteInt (LOCALINI_FILE, "OldCpuSpeed", (int) (CPU_SPEED + 0.5));
		if (old_cpu_speed) {
			ROLLING_AVERAGE = 1000;
			IniWriteInt (LOCALINI_FILE, "RollingAverage", 1000);
			IniWriteInt (LOCALINI_FILE, "RollingStartTime", 0);
			spoolMessage (PRIMENET_SET_COMPUTER_INFO, NULL);
			UpdateEndDates ();
		}
	}
}

/* Set the CPU flags based on the CPUID instruction.  Also, the advanced */
/* user can override our guesses. */

void getCpuInfo (void)
{
	int	temp, old_cpu_type;

/* Get the CPU info using CPUID instruction */

	guessCpuType ();

/* Deduce the old-style CPU_TYPE variable given the CPU brand string */
/* returned by CPUID.  Our code should no longer use CPU_TYPE, but rather */
/* the individual CPU capability flags in CPU_FLAGS.  However, the old */
/* CPU_TYPE is still sent to the server for reporting reasons. */

	if (strstr (CPU_BRAND, "Intel")) {
		if (strstr (CPU_BRAND, "486")) CPU_TYPE = 4;
		else if (strstr (CPU_BRAND, "Pro")) CPU_TYPE = 6;
		else if (strstr (CPU_BRAND, "Celeron")) CPU_TYPE = 8;
		else if (strstr (CPU_BRAND, "III")) CPU_TYPE = 10;
		else if (strstr (CPU_BRAND, "II")) CPU_TYPE = 9;
		else if (CPU_FLAGS & CPU_SSE2) CPU_TYPE = 12;
		else if (strstr (CPU_BRAND, "Pentium")) CPU_TYPE = 5;
		else CPU_TYPE = 4;
	} else if (strstr (CPU_BRAND, "AMD")) {
		if (strstr (CPU_BRAND, "486")) CPU_TYPE = 3;
		else if (strstr (CPU_BRAND, "Unknown")) CPU_TYPE = 3;
		else if (strstr (CPU_BRAND, "K5")) CPU_TYPE = 3;
		else if (strstr (CPU_BRAND, "K6")) CPU_TYPE = 7;
		else CPU_TYPE = 11;
	} else
		CPU_TYPE = 3;

/* Now let the user override the cpu type and speed from the local.ini file */

	if (IniGetInt (LOCALINI_FILE, "CpuOverride", 0)) {
		CPU_TYPE = IniGetInt (LOCALINI_FILE, "CpuType", CPU_TYPE);
	}

/* Let the user override the cpu flags from the local.ini file */

	temp = IniGetInt (LOCALINI_FILE, "CpuSupportsRDTSC", 99);
	if (temp == 0) CPU_FLAGS &= ~CPU_RDTSC;
	if (temp == 1) CPU_FLAGS |= CPU_RDTSC;
	temp = IniGetInt (LOCALINI_FILE, "CpuSupportsCMOV", 99);
	if (temp == 0) CPU_FLAGS &= ~CPU_CMOV;
	if (temp == 1) CPU_FLAGS |= CPU_CMOV;
	temp = IniGetInt (LOCALINI_FILE, "CpuSupportsPrefetch", 99);
	if (temp == 0) CPU_FLAGS &= ~CPU_PREFETCH;
	if (temp == 1) CPU_FLAGS |= CPU_PREFETCH;
	temp = IniGetInt (LOCALINI_FILE, "CpuSupportsSSE", 99);
	if (temp == 0) CPU_FLAGS &= ~CPU_SSE;
	if (temp == 1) CPU_FLAGS |= CPU_SSE;
	temp = IniGetInt (LOCALINI_FILE, "CpuSupportsSSE2", 99);
	if (temp == 0) CPU_FLAGS &= ~CPU_SSE2;
	if (temp == 1) CPU_FLAGS |= CPU_SSE2;

/* Let the user override the L2 cache size in local.ini file */

	CPU_L2_CACHE_SIZE =
		IniGetInt (LOCALINI_FILE, "CpuL2CacheSize", CPU_L2_CACHE_SIZE);
	CPU_L2_CACHE_LINE_SIZE =
		IniGetInt (LOCALINI_FILE, "CpuL2CacheLineSize", CPU_L2_CACHE_LINE_SIZE);

/* If CPU has changed greatly since the last time prime95 ran, then */
/* tell the server, recalculate new completion dates, and reset the */
/* rolling average.  However, don't do this on the first call where */
/* the Welcome screen has yet to be displayed. */

	old_cpu_type = IniGetInt (LOCALINI_FILE, "OldCpuType", 0);
	if ((int) CPU_TYPE != old_cpu_type) {
		IniWriteInt (LOCALINI_FILE, "OldCpuType", CPU_TYPE);
		if (old_cpu_type) {
			ROLLING_AVERAGE = 1000;
			IniWriteInt (LOCALINI_FILE, "RollingAverage", 1000);
			IniWriteInt (LOCALINI_FILE, "RollingStartTime", 0);
			spoolMessage (PRIMENET_SET_COMPUTER_INFO, NULL);
			UpdateEndDates ();
		}
	}

/* Now get the CPU speed */

	getCpuSpeed ();
}

/* Format a long or very long textual cpu description */

void getCpuDescription (
	char	*buf,			/* A 512 byte buffer */
	int	long_desc)		/* True for a very long description */
{

/* Recalculate the CPU speed in case speed step has changed the original */
/* settings. */

	getCpuSpeed ();

/* Now format a pretty CPU description */

	sprintf (buf, "%s\nCPU speed: %.2f MHz\n", CPU_BRAND, CPU_SPEED);
	if (CPU_FLAGS) {
		strcat (buf, "CPU features: ");
		if (CPU_FLAGS & CPU_RDTSC) strcat (buf, "RDTSC, ");
		if (CPU_FLAGS & CPU_CMOV) strcat (buf, "CMOV, ");
		if (CPU_FLAGS & CPU_PREFETCH) strcat (buf, "PREFETCH, ");
		if (CPU_FLAGS & CPU_MMX) strcat (buf, "MMX, ");
		if (CPU_FLAGS & CPU_SSE) strcat (buf, "SSE, ");
		if (CPU_FLAGS & CPU_SSE2) strcat (buf, "SSE2, ");
		strcpy (buf + strlen (buf) - 2, "\n");
	}
	strcat (buf, "L1 cache size: ");
	if (CPU_L1_CACHE_SIZE < 0) strcat (buf, "unknown\n");
	else sprintf (buf + strlen (buf), "%d KB\n", CPU_L1_CACHE_SIZE);
	strcat (buf, "L2 cache size: ");
	if (CPU_L2_CACHE_SIZE < 0) strcat (buf, "unknown\n");
	else sprintf (buf + strlen (buf), "%d KB\n", CPU_L2_CACHE_SIZE);
	if (! long_desc) return;
	strcat (buf, "L1 cache line size: ");
	if (CPU_L1_CACHE_LINE_SIZE < 0) strcat (buf, "unknown\n");
	else sprintf (buf+strlen(buf), "%d bytes\n", CPU_L1_CACHE_LINE_SIZE);
	strcat (buf, "L2 cache line size: ");
	if (CPU_L2_CACHE_LINE_SIZE < 0) strcat (buf, "unknown\n");
	else sprintf (buf+strlen(buf), "%d bytes\n", CPU_L2_CACHE_LINE_SIZE);
	if (CPU_L1_DATA_TLBS > 0)
		sprintf (buf + strlen (buf), "L1 TLBS: %d\n", CPU_L1_DATA_TLBS);
	if (CPU_L2_DATA_TLBS > 0)
		sprintf (buf + strlen (buf), "%sTLBS: %d\n",
			 CPU_L1_DATA_TLBS > 0 ? "L2 " : "",
			 CPU_L2_DATA_TLBS);
}

/* Determine if a number is prime */

int isPrime (
	unsigned long p)
{
	unsigned long i;
	for (i = 2; i * i <= p; i = (i + 1) | 1)
		if (p % i == 0) return (FALSE);
	return (TRUE);
}

/* Return maximum available memory */

unsigned int max_mem (void)
{
	if (DAY_MEMORY > NIGHT_MEMORY) return (DAY_MEMORY);
	return (NIGHT_MEMORY);
}

/* Upper case a string */

void strupper (
	char	*p)
{
	for ( ; *p; p++) if (*p >= 'a' && *p <= 'z') *p = *p - 'a' + 'A';
}

/* Convert a string (e.g "11:30 AM") to minutes since midnight */

unsigned int strToMinutes (
	char	*buf)
{
	unsigned int hours, minutes, pm;

	pm = (strchr (buf, 'P') != NULL || strchr (buf, 'p') != NULL);
	hours = atoi (buf);
	while (isdigit (*buf)) buf++;
	while (*buf && ! isdigit (*buf)) buf++;
	minutes = atoi (buf);
	if (hours == 12) hours -= 12;
	if (pm) hours += 12;
	minutes = hours * 60 + minutes;
	if (minutes > 1440) minutes = 1440;
	return (minutes);
}

/* Convert minutes since midnight to a string (e.g "11:30 AM") */

void minutesToStr (
	unsigned int minutes,
	char	*buf)
{
	unsigned int fmt_type, hours, pm;

	hours = minutes / 60;
	fmt_type = IniGetInt (INI_FILE, "AMPM", 0);
	if (fmt_type == 0) fmt_type = getDefaultTimeFormat ();
	if (fmt_type != 1) {
		sprintf (buf, "%d:%02d", hours, minutes % 60);
	} else {
		if (hours >= 12) hours -= 12, pm = 1;
		else pm = 0;
		if (hours == 0) hours = 12;
		sprintf (buf, "%d:%02d %s", hours, minutes % 60, pm ? "PM" : "AM");
	}
}

/* Determine the names of the INI files */

void nameIniFiles (
	int	named_ini_files)
{
	char	buf[120];

	if (named_ini_files < 0) {
		strcpy (INI_FILE, "prime.ini");
		strcpy (LOCALINI_FILE, "local.ini");
		strcpy (SPOOL_FILE, "prime.spl");
		strcpy (WORKTODO_FILE, "worktodo.ini");
		strcpy (RESFILE, "results.txt");
		strcpy (LOGFILE, "prime.log");
		strcpy (EXTENSION, "");
	} else {
		sprintf (INI_FILE, "prim%04d.ini", named_ini_files);
		sprintf (LOCALINI_FILE, "loca%04d.ini", named_ini_files);
		sprintf (SPOOL_FILE, "prim%04d.spl", named_ini_files);
		sprintf (WORKTODO_FILE, "work%04d.ini", named_ini_files);
		sprintf (RESFILE, "resu%04d.txt", named_ini_files);
		sprintf (LOGFILE, "prim%04d.log", named_ini_files);
		sprintf (EXTENSION, ".%03d", named_ini_files);
	}

/* Let the user rename these files and pick a different working directory */

	IniGetString (INI_FILE, "WorkingDir", buf, sizeof(buf), NULL);
	IniGetString (INI_FILE, "local.ini", LOCALINI_FILE, 80, LOCALINI_FILE);
	IniGetString (INI_FILE, "prime.spl", SPOOL_FILE, 80, SPOOL_FILE);
	IniGetString (INI_FILE, "worktodo.ini", WORKTODO_FILE, 80, WORKTODO_FILE);
	IniGetString (INI_FILE, "results.txt", RESFILE, 80, RESFILE);
	IniGetString (INI_FILE, "prime.log", LOGFILE, 80, LOGFILE);
	IniGetString (INI_FILE, "prime.ini", INI_FILE, 80, INI_FILE);
	if (buf[0]) {
		_chdir (buf);
		IniFileOpen (INI_FILE, 0);
	}
}

/* Read the INI files */

void readIniFiles (void)
{
	int	temp;
	char	buf[512];

	processTimedIniFile (INI_FILE);
	IniGetString (INI_FILE, "UserID", USERID, sizeof (USERID), NULL);
	USERID[PRIMENET_USER_ID_LENGTH] = 0;
	IniGetString (INI_FILE, "UserPWD", USER_PWD, sizeof (USER_PWD), NULL);
	USER_PWD[PRIMENET_USER_PW_LENGTH] = 0;
	IniGetString (INI_FILE, "OldUserID", OLD_USERID, sizeof (OLD_USERID), NULL);
	OLD_USERID[PRIMENET_USER_ID_LENGTH] = 0;
	IniGetString (INI_FILE, "OldUserPWD", OLD_USER_PWD, sizeof (OLD_USER_PWD), NULL);
	OLD_USER_PWD[PRIMENET_USER_PW_LENGTH] = 0;
	IniGetString (INI_FILE, "UserName", USER_NAME,
		      sizeof (USER_NAME), NULL);
	IniGetString (INI_FILE, "UserEmailAddr", USER_ADDR,
		      sizeof (USER_ADDR), NULL);
	IniGetString (LOCALINI_FILE, "ComputerID", COMPID,
		      sizeof (COMPID), NULL);
	COMPID[PRIMENET_COMPUTER_ID_LENGTH] = 0;
	NEWSLETTERS = (int) IniGetInt (INI_FILE, "Newsletters", 0);
	USE_PRIMENET = (int) IniGetInt (INI_FILE, "UsePrimenet", 1);
	DIAL_UP = (int) IniGetInt (INI_FILE, "DialUp", 0);
	DAYS_OF_WORK = (unsigned int) IniGetInt (INI_FILE, "DaysOfWork", 5);
	if (DAYS_OF_WORK > 180) DAYS_OF_WORK = 180;
	WORK_PREFERENCE = (short) IniGetInt (INI_FILE, "WorkPreference", 0);

        VACATION_END = IniGetInt (LOCALINI_FILE, "VacationEnd", 0);
        ON_DURING_VACATION = (int) IniGetInt (LOCALINI_FILE, "VacationOn", 1);

	CPU_TYPE = (unsigned int) IniGetInt (LOCALINI_FILE, "CPUType", 0);
	CPU_SPEED = IniGetInt (LOCALINI_FILE, "CPUSpeed", 0);
	CPU_HOURS = (unsigned int) IniGetInt (LOCALINI_FILE, "CPUHours", 24);
	if (CPU_HOURS < 1) CPU_HOURS = 1;
	if (CPU_HOURS > 24) CPU_HOURS = 24;
	temp = IniGetInt (LOCALINI_FILE, "DayMemory", 8);
	DAY_MEMORY = (temp < 0) ? physical_memory () + temp : temp;
	if (DAY_MEMORY < 8) DAY_MEMORY = 8;
	if (DAY_MEMORY > physical_memory () - 8)
		DAY_MEMORY = physical_memory () - 8;
	temp = IniGetInt (LOCALINI_FILE, "NightMemory", 8);
	NIGHT_MEMORY = (temp < 0) ? physical_memory () + temp : temp;
	if (NIGHT_MEMORY < 8) NIGHT_MEMORY = 8;
	if (NIGHT_MEMORY > physical_memory () - 8)
		NIGHT_MEMORY = physical_memory () - 8;
	DAY_START_TIME = IniGetInt (LOCALINI_FILE, "DayStartTime", 450);
	DAY_END_TIME = IniGetInt (LOCALINI_FILE, "DayEndTime", 1410);

	ROLLING_AVERAGE = (unsigned int) IniGetInt (LOCALINI_FILE, "RollingAverage", 1000);
	if (ROLLING_AVERAGE < 10) ROLLING_AVERAGE = 10;
	if (ROLLING_AVERAGE > 4000) ROLLING_AVERAGE = 4000;

	PRECISION = (unsigned int) IniGetInt (INI_FILE, "PercentPrecision", 2);
	if (PRECISION > 6) PRECISION = 6;

	ITER_OUTPUT = IniGetInt (INI_FILE, "OutputIterations", 10000);
	if (ITER_OUTPUT > 999999999) ITER_OUTPUT = 999999999;
	if (ITER_OUTPUT <= 0) ITER_OUTPUT = 1;
	ITER_OUTPUT_RES = IniGetInt (INI_FILE, "ResultsFileIterations", 999999999);
	if (ITER_OUTPUT_RES > 999999999) ITER_OUTPUT_RES = 999999999;
	if (ITER_OUTPUT_RES < 1000) ITER_OUTPUT_RES = 1000;
	DISK_WRITE_TIME = IniGetInt (INI_FILE, "DiskWriteTime", 30);
	MODEM_RETRY_TIME = (unsigned int) IniGetInt (INI_FILE, "NetworkRetryTime", 2);
	if (MODEM_RETRY_TIME < 1) MODEM_RETRY_TIME = 1;
	if (MODEM_RETRY_TIME > 300) MODEM_RETRY_TIME = 300;
	NETWORK_RETRY_TIME = (unsigned int)
		IniGetInt (INI_FILE, "NetworkRetryTime2",
			   MODEM_RETRY_TIME > 60 ? MODEM_RETRY_TIME : 60);
	if (NETWORK_RETRY_TIME < 1) NETWORK_RETRY_TIME = 1;
	if (NETWORK_RETRY_TIME > 300) NETWORK_RETRY_TIME = 300;
	DAYS_BETWEEN_CHECKINS = (unsigned int) IniGetInt (INI_FILE, "DaysBetweenCheckins",28);
	if (DAYS_BETWEEN_CHECKINS > 60) DAYS_BETWEEN_CHECKINS = 60;
	if (DAYS_BETWEEN_CHECKINS < 1) DAYS_BETWEEN_CHECKINS = 1;
	TWO_BACKUP_FILES = (int) IniGetInt (INI_FILE, "TwoBackupFiles", 1);
	SILENT_VICTORY = (int) IniGetInt (INI_FILE, "SilentVictory", 0);
	RUN_ON_BATTERY = (int) IniGetInt (LOCALINI_FILE, "RunOnBattery", 1);

	ADVANCED_ENABLED = (int) IniGetInt (INI_FILE, "Advanced", 0);
	temp = (int) IniGetInt (INI_FILE, "ErrorCheck", 0);
	ERRCHK = (temp != 0);
	PRIORITY = (unsigned int) IniGetInt (INI_FILE, "Priority", 1);
	CPU_AFFINITY = (unsigned int) IniGetInt (INI_FILE, "Affinity", 99);
	if (CPU_AFFINITY != 99) {
		IniWriteString (INI_FILE, "Affinity", NULL);
		IniWriteInt (LOCALINI_FILE, "Affinity", CPU_AFFINITY);
	} else {
		CPU_AFFINITY = (unsigned int)
			IniGetInt (LOCALINI_FILE, "Affinity", 99);
	}
	MANUAL_COMM = (int) IniGetInt (INI_FILE, "ManualComm", 0);
	HIDE_ICON = (int) IniGetInt (INI_FILE, "HideIcon", 0);
	TRAY_ICON = (int) IniGetInt (INI_FILE, "TrayIcon", 1);

/* Get the CPU type, speed, and capabilities. */

	getCpuInfo ();

/* Get the option controlling which timer to use.  If the high resolution */
/* performance counter is not available on this machine, then add 10 to */
/* the RDTSC_TIMING value. */

	RDTSC_TIMING = IniGetInt (INI_FILE, "RdtscTiming", 1);
	if (RDTSC_TIMING < 10 && ! isHighResTimerAvailable ())
		RDTSC_TIMING += 10;

/* Other oddball options */

	TIMESTAMPING = IniGetInt (INI_FILE, "TimeStamp", 1);
	CUMULATIVE_TIMING = IniGetInt (INI_FILE, "CumulativeTiming", 0);
	WELL_BEHAVED_WORK = IniGetInt (INI_FILE, "WellBehavedWork", 0);
	IniGetString (INI_FILE, "PauseWhileRunning", buf, sizeof (buf), NULL);
	if (buf[0]) {
		char	*p;
		int	i, cnt;
		for (cnt = 1, p = buf; p = strchr (p, ','); cnt++) *p++ = 0;
		PAUSE_WHILE_RUNNING = (char **)
			malloc ((cnt + 1) * sizeof (char *));
		for (i = 0, p = buf; i < cnt; i++) {
			PAUSE_WHILE_RUNNING[i] = (char *)
				malloc (strlen (p) + 1);
			strupper (p);
			strcpy (PAUSE_WHILE_RUNNING[i], p);
			p += strlen (p) + 1;
		}
		PAUSE_WHILE_RUNNING[i] = NULL;
	} else
		PAUSE_WHILE_RUNNING = NULL;
}

/*----------------------------------------------------------------------
| Portable routines to read and write ini files!  NOTE:  These only
| work if you open no more than 5 ini files.  Also you must not
| change the working directory at any time during program execution.
+---------------------------------------------------------------------*/

struct IniLine {
	char	*keyword;
	char	*value;
	int	active;
};
struct IniCache {
	char	*filename;
	int	immediate_writes;
	int	dirty;
	unsigned int num_lines;
	unsigned int array_size;
	struct IniLine **lines;
};

void growIniLineArray (
	struct IniCache *p)
{
	struct IniLine **newlines;

	if (p->num_lines != p->array_size) return;

	newlines = (struct IniLine **)
		malloc ((p->num_lines + 100) * sizeof (struct IniLine **));
	if (p->num_lines) {
		memcpy (newlines, p->lines, p->num_lines * sizeof (struct IniLine *));
		free (p->lines);
	}
	p->lines = newlines;
	p->array_size = p->num_lines + 100;
}

struct IniCache *openIniFile (
	char	*filename,
	int	forced_read)
{
static	struct IniCache *cache[10] = {0};
	struct IniCache *p;
	FILE	*fd;
	unsigned int i;
	char	line[80];
	char	*val;

/* See if file is cached */

	for (i = 0; i < 10; i++) {
		p = cache[i];
		if (p == NULL) {
			p = (struct IniCache *) malloc (sizeof (struct IniCache));
			p->filename = (char *) malloc (strlen (filename) + 1);
			strcpy (p->filename, filename);
			p->immediate_writes = 1;
			p->dirty = 0;
			p->num_lines = 0;
			p->array_size = 0;
			p->lines = NULL;
			forced_read = 1;
			cache[i] = p;
			break;
		}
		if (strcmp (filename, p->filename) == 0)
			break;
	}

/* Skip reading the ini file if appropriate */

	if (!forced_read) return (p);
	if (p->dirty) return (p);

/* Free the data if we've already read some in */

	for (i = 0; i < p->num_lines; i++) {
		free (p->lines[i]->keyword);
		free (p->lines[i]->value);
		free (p->lines[i]);
	}
	p->num_lines = 0;

/* Read the IniFile */
	
	fd = fopen (filename, "r");
	if (fd == NULL) return (p);

	while (fgets (line, 80, fd)) {
		if (line[strlen(line)-1] == '\n') line[strlen(line)-1] = 0;
		if (line[0] == 0) continue;
		if (line[strlen(line)-1] == '\r') line[strlen(line)-1] = 0;
		if (line[0] == 0) continue;

		if (line[0] == ';') continue;
		if (line[0] == '#') continue;
		if (line[0] == '[') continue;

		val = strchr (line, '=');
		if (val == NULL) {
			char	buf[130];
			sprintf (buf, "Illegal line in INI file: %s\n", line);
			OutputSomewhere (buf);
			continue;
		}
		*val++ = 0;

		growIniLineArray (p);
		
/* Allocate and fill in a new line structure */

		i = p->num_lines++;
		p->lines[i] = (struct IniLine *) malloc (sizeof (struct IniLine));
		p->lines[i]->keyword = (char *) malloc (strlen (line) + 1);
		p->lines[i]->value = (char *) malloc (strlen (val) + 1);
		p->lines[i]->active = TRUE;
		strcpy (p->lines[i]->keyword, line);
		strcpy (p->lines[i]->value, val);
	}
	fclose (fd);

	return (p);
}

void writeIniFile (
	struct IniCache *p)
{
	int	fd;
	unsigned int j;
	char	buf[100];

/* Delay writing the file unless this INI file is written */
/* to immediately */

	if (!p->immediate_writes) {
		p->dirty = 1;
		return;
	}

/* Create and write out the INI file */

	fd = _open (p->filename, _O_CREAT | _O_TRUNC | _O_WRONLY | _O_TEXT, 0666);
	if (fd < 0) return;
	for (j = 0; j < p->num_lines; j++) {
		strcpy (buf, p->lines[j]->keyword);
		strcat (buf, "=");
		strcat (buf, p->lines[j]->value);
		strcat (buf, "\n");
		_write (fd, buf, strlen (buf));
	}
	p->dirty = 0;
	_close (fd);
}

/* Routines to help analyze a Time= line in an INI file */

void parseTimeLine (
	char	**line,
	int	*start_day,
	int	*end_day,
	int	*start_time,
	int	*end_time)
{
	char	*p;

/* Get the days of the week, e.g. 1-5 */

	p = *line;
	*start_day = atoi (p); while (isdigit (*p)) p++;
	if (*p == '-') {
		p++;
		*end_day = atoi (p); while (isdigit (*p)) p++;
	} else
		*end_day = *start_day;

/* Now do time portion.  If none present, then assume the numbers we */
/* parsed above were times, not days of the week. */

	if (*p == '/')
		p++;
	else {
		p = *line;
		*start_day = 1;
		*end_day = 7;
	} 
	*start_time = atoi (p) * 60; while (isdigit (*p)) p++;
	if (*p == ':') {
		p++;
		*start_time += atoi (p); while (isdigit (*p)) p++;
	}
	if (*p == '-') p++;			/* Skip '-' */
	*end_time = atoi (p) * 60; while (isdigit (*p)) p++;
	if (*p == ':') {
		p++;
		*end_time += atoi (p); while (isdigit (*p)) p++;
	}

/* Return ptr to next time interval on the line */

	if (*p == ',') p++;
	*line = p;
}

int analyzeTimeLine (
	char	*line,
	time_t	current_t)
{
	char	*p;
	struct tm *x;
	int	current_time;
	int	day, start_day, end_day, start_time, end_time;

/* Break current time into a more easily maniupulated form */

	x = localtime (&current_t);
	current_time = (x->tm_wday ? x->tm_wday : 7) * 24 * 60;
	current_time += x->tm_hour * 60 + x->tm_min;

/* Process each interval on the line */

	p = line;
	while (*p) {
		parseTimeLine (&p, &start_day, &end_day, &start_time, &end_time);

/* Treat each day in the range as a separate time interval to process */

		for (day = start_day; day <= end_day; day++) {
			int	temp;
			time_t	wakeup_t;

/* Is the current time in this interval? */

			if (current_time >= day * 24 * 60 + start_time &&
			    current_time < day  * 24 * 60 + end_time)
				goto winner;

/* No, see if this start time should be our new wakeup time. */

			temp = day * 24 * 60 + start_time;
			if (temp >= current_time)
				wakeup_t = current_t + (temp - current_time) * 60;
			else
				wakeup_t = current_t + (temp + 7 * 24 * 60 - current_time) * 60;
			if (SLEEP_TIME == 0 || SLEEP_TIME > wakeup_t)
				SLEEP_TIME = wakeup_t;
		}
	}

/* Current time was not in any of the intervals */

	return (FALSE);

/* Current time is in this interval, then set the END_TIME. */

winner:	END_TIME = current_t + (day * 24 * 60 + end_time - current_time) * 60;

/* Also, look for a start time that matches the end time and replace */
/* the end time.  For example, if current time is 18:00 and the */
/* Time= entry is 0:00-8:00,17:00-24:00, then the */
/* end time of 24:00 should be replaced with 8:00 of the next day. */

	p = line;
	while (*p) {
		parseTimeLine (&p, &start_day, &end_day, &start_time, &end_time);

/* Treat each day in the range as a separate time interval to process */

		for (day = start_day; day <= end_day; day++) {
			int	temp;
			time_t	start_t;

/* If this start time is the same as the END_TIME, then set the new */
/* END_TIME to be the end of this interval */

			temp = day * 24 * 60 + start_time;
			if (temp >= current_time)
				start_t = current_t + (temp - current_time) * 60;
			else
				start_t = current_t + (temp + 7 * 24 * 60 - current_time) * 60;
			if (END_TIME == start_t) {
				END_TIME += (end_time - start_time) * 60;
				p = line;
				break;
			}
		}
	}

/* Return indicator that current time was covered by one of the intervals */

	return (TRUE);
}

/* INI files can contain conditional sections bracketed with Time= lines */
/* that tell when that section should be active.  For example, this INI */
/* file has different properties during the work week, and sleeps for an */
/* hour each workday. */
/*	UserID=foo						*/
/*	Time=1-5/8:30-17:30					*/
/*	Priority=1						*/
/*	Time=1-5/0:00-8:30,1-5/17:30-23:00,6-7/0:00-24:00	*/
/*	Priority=5						*/

void processTimedIniFile (
	char	*filename)
{
	struct IniCache *p;
	time_t	current_time;
	unsigned int i;
	int	active;

/* Open ini file */

	p = openIniFile (filename, 0);

/* Get the current time - to compare to any Time= lines */
	
	time (&current_time);

/* Assume this will be an untimed run */

	END_TIME = 0;
	SLEEP_TIME = 0;

/* Process lines until we run into next Time= line. */
/* Mark any lines prior to the first Time= line active. */

	active = TRUE;
	for (i = 0; i < p->num_lines; i++) {
		p->lines[i]->active = active;

/* If this is a Time= line, see if the current time is part of the */
/* interval described in the Time= line.  If the current time is part */
/* of the interval then set following lines active, otherwise set the */
/* following lines inactive. */
		
		if (stricmp (p->lines[i]->keyword, "Time")) continue;
		active = analyzeTimeLine (p->lines[i]->value, current_time);
	}

/* If there were Time= lines where the current time was included in a */
/* time interval, then we won't be sleeping. */

	if (END_TIME) SLEEP_TIME = 0;
}

void truncated_strcpy (
	char	*buf,
	unsigned int bufsize,
	char	*val)
{
	if (strlen (val) >= bufsize) {
		memcpy (buf, val, bufsize-1);
		buf[bufsize-1] = 0;
	} else {
		strcpy (buf, val);
	}
}

void IniGetString (
	char	*filename,
	char	*keyword,
	char	*val,
	unsigned int val_bufsize,
	char	*default_val)
{
	struct IniCache *p;
	unsigned int i;

/* Open ini file */

	p = openIniFile (filename, 0);

/* Look for the keyword */

	for (i = 0; ; i++) {
		if (i == p->num_lines) {
			if (default_val == NULL) {
				val[0] = 0;
			} else {
				truncated_strcpy (val, val_bufsize, default_val);
			}
			return;
		}
		if (p->lines[i]->active &&
		    stricmp (keyword, p->lines[i]->keyword) == 0) break;
	}

/* Copy info from the line structure to the user buffers */

	truncated_strcpy (val, val_bufsize, p->lines[i]->value);
}

long IniGetInt (
	char	*filename,
	char	*keyword,
	long	default_val)
{
	char	buf[20], defval[20];
	sprintf (defval, "%ld", default_val);
	IniGetString (filename, keyword, buf, 20, defval);
	return (atol (buf));
}

void IniWriteString (
	char	*filename,
	char	*keyword,
	char	*val)
{
	struct IniCache *p;
	unsigned int i, j;

/* Open ini file */

	p = openIniFile (filename, 1);

/* Look for the keyword */

	for (i = 0; ; i++) {
		if (i == p->num_lines ||
		    stricmp (p->lines[i]->keyword, "Time") == 0) {

/* Ignore request if we are deleting line */

			if (val == NULL) return;

/* Make sure the line array has room for the new line */

			growIniLineArray (p);

/* Shuffle entries down to make room for this entry */

			for (j = p->num_lines; j > i; j--)
				p->lines[j] = p->lines[j-1];

/* Allocate and fill in a new line structure */

			p->lines[i] = (struct IniLine *) malloc (sizeof (struct IniLine));
			p->lines[i]->keyword = (char *) malloc (strlen (keyword) + 1);
			strcpy (p->lines[i]->keyword, keyword);
			p->lines[i]->value = NULL;
			p->num_lines++;
			break;
		}
		if (p->lines[i]->active &&
		    stricmp (keyword, p->lines[i]->keyword) == 0) {
			if (val != NULL && strcmp (val, p->lines[i]->value) == 0) return;
			break;
		}
	}

/* Delete the line if requested */

	if (val == NULL) {
		IniDeleteLine (filename, i+1);
		return;
	}

/* Replace the value associated with the keyword */

	free (p->lines[i]->value);
	p->lines[i]->value = (char *) malloc (strlen (val) + 1);
	strcpy (p->lines[i]->value, val);

/* Write the INI file back to disk */

	writeIniFile (p);
}

void IniWriteInt (
	char	*filename,
	char	*keyword,
	long	val)
{
	char	buf[20];
	sprintf (buf, "%ld", val);
	IniWriteString (filename, keyword, buf);
}

void IniFileOpen (
	char	*filename,
	int	immediate_writes)
{
	struct IniCache *p;
	p = openIniFile (filename, 1);
	p->immediate_writes = immediate_writes;
}

void IniFileClose (
	char	*filename)
{
	struct IniCache *p;
	p = openIniFile (filename, 0);
	if (p->dirty) {
		p->immediate_writes = 1;
		writeIniFile (p);
		p->immediate_writes = 0;
	}
}

int IniFileWritable (
	char	*filename)
{
	struct IniCache *p;
	int	fd;
	unsigned int j;
	char	buf[100];

/* Create and write out the INI file */

	p = openIniFile (filename, 0);
	fd = _open (p->filename, _O_CREAT | _O_TRUNC | _O_WRONLY | _O_TEXT, 0666);
	if (fd < 0) return (FALSE);
	for (j = 0; j < p->num_lines; j++) {
		strcpy (buf, p->lines[j]->keyword);
		strcat (buf, "=");
		strcat (buf, p->lines[j]->value);
		strcat (buf, "\n");
		if (_write (fd, buf, strlen (buf)) != (int) strlen (buf)) {
			_close (fd);
			return (FALSE);
		}
	}
	if (p->num_lines == 0) {
		if (_write (fd, "DummyLine=XXX\n", 14) != 14) {
			_close (fd);
			return (FALSE);
		}
		p->dirty = 1;
	}
	_close (fd);
	return (TRUE);
}

unsigned int IniGetNumLines (
	char	*filename)
{
	struct IniCache *p;
	p = openIniFile (filename, 0);
	return (p->num_lines);
}

void IniGetLineAsString (
	char	*filename,
	unsigned int line,
	char	*keyword,
	unsigned int keyword_bufsize,
	char	*val,
	unsigned int val_bufsize)
{
	struct IniCache *p;

/* Open ini file */

	p = openIniFile (filename, 0);

/* Copy info from the line structure to the user buffers */

	truncated_strcpy (keyword, keyword_bufsize, p->lines[line-1]->keyword);
	truncated_strcpy (val, val_bufsize, p->lines[line-1]->value);
}

void IniGetLineAsInt (
	char	*filename,
	unsigned int line,
	char	*keyword,
	unsigned int keyword_bufsize,
	long	*val)
{
	char	buf[20];
	IniGetLineAsString (filename, line, keyword, keyword_bufsize, buf, 20);
	*val = atol (buf);
}

void IniReplaceLineAsString (
	char	*filename,
	unsigned int line,
	char	*keyword,
	char	*val)
{
	IniDeleteLine (filename, line);
	IniInsertLineAsString (filename, line, keyword, val);
}

void IniReplaceLineAsInt (
	char	*filename,
	unsigned int line,
	char	*keyword,
	long	val)
{
	char	buf[20];
	sprintf (buf, "%ld", val);
	IniReplaceLineAsString (filename, line, keyword, buf);
}

void IniInsertLineAsString (
	char	*filename,
	unsigned int line,
	char	*keyword,
	char	*val)
{
	struct IniCache *p;
	unsigned int i;

/* Open ini file, do not reread it as that could change the line numbers! */

	p = openIniFile (filename, 0);

/* Adjust line number if it doesn't make sense */

	if (line == 0) line = 1;
	if (line > p->num_lines+1) line = p->num_lines+1;

/* Make sure the line array has room for the new line */

	growIniLineArray (p);

/* Shuffle lines down in the array to make room for the new line */

	for (i = p->num_lines; i >= line; i--) p->lines[i] = p->lines[i-1];
	p->num_lines++;

/* Allocate and fill in a new line structure */

	p->lines[line-1] = (struct IniLine *) malloc (sizeof (struct IniLine));
	p->lines[line-1]->keyword = (char *) malloc (strlen (keyword) + 1);
	p->lines[line-1]->value = (char *) malloc (strlen (val) + 1);
	p->lines[line-1]->active = TRUE;
	strcpy (p->lines[line-1]->keyword, keyword);
	strcpy (p->lines[line-1]->value, val);

/* Write the INI file back to disk */

	writeIniFile (p);
}

void IniInsertLineAsInt (
	char	*filename,
	unsigned int line,
	char	*keyword,
	long	val)
{
	char	buf[20];
	sprintf (buf, "%ld", val);
	IniInsertLineAsString (filename, line, keyword, buf);
}

void IniAppendLineAsString (
	char	*filename,
	char	*keyword,
	char	*val)
{
	struct IniCache *p;
	p = openIniFile (filename, 0);
	IniInsertLineAsString (filename, p->num_lines+1, keyword, val);
}

void IniAppendLineAsInt (
	char	*filename,
	char	*keyword,
	long	val)
{
	char	buf[20];
	sprintf (buf, "%ld", val);
	IniAppendLineAsString (filename, keyword, buf);
}

void IniDeleteLine (
	char	*filename,
	unsigned int line)
{
	struct IniCache *p;
	unsigned int i;

/* Open ini file, do not reread it as that could change the line numbers! */

	p = openIniFile (filename, 0);
	if (line == 0 || line > p->num_lines) return;

/* Free the data associated with the given line */

	free (p->lines[line-1]->keyword);
	free (p->lines[line-1]->value);
	free (p->lines[line-1]);

/* Delete the line from the lines array */

	for (i = line; i < p->num_lines; i++) p->lines[i-1] = p->lines[i];
	p->num_lines--;

/* Write the INI file back to disk */

	writeIniFile (p);
}

void IniDeleteAllLines (
	char	*filename)
{
	struct IniCache *p;
	unsigned int i;

/* Open ini file! */

	p = openIniFile (filename, 0);

/* Free the data associated with the given line */

	for (i = 0; i < p->num_lines; i++) {
		free (p->lines[i]->keyword);
		free (p->lines[i]->value);
		free (p->lines[i]);
	}
	p->num_lines = 0;

/* Write the INI file back to disk */

	writeIniFile (p);
}

/*----------------------------------------------------------------------+
| Portable routines to read and write messages in the spool file.	|
+----------------------------------------------------------------------*/

#include "security.c"
#ifndef IPSHASH
void hash_packet (short operation, void *pkt)
{
	struct primenetUserInfo *z = (struct primenetUserInfo *) pkt;
	if (operation == PRIMENET_PING_SERVER_INFO) return;
	z->salt = 0;
	z->hash = 0;
}
#endif

/* Update completion dates on the server.  Set a flag in */
/* the spool file saying this is necessary. */

void UpdateEndDates (void)
{
	spoolMessage (PRIMENET_COMPLETION_DATE, NULL);
}

/* Update completion dates on the server if it has been a */
/* month since we last updated the server. */

void ConditionallyUpdateEndDates (void)
{
#ifndef SERVER_TESTING
	time_t	start_time, current_time;

/* Get the current time and when the completion dates were last sent */

	time (&current_time);
	start_time = IniGetInt (LOCALINI_FILE, "LastEndDatesSent", 0);

/* If it's been the correct number of days, then update the end dates */

	if (current_time < start_time ||
	    current_time > (time_t) (start_time + DAYS_BETWEEN_CHECKINS * 86400))
		UpdateEndDates ();
#endif
}

/* Write a message to the spool file */

void spoolMessage (
	short	msgType,
	void	*msg)
{
	int	fd;
	char	header_byte;

/* If we're not using primenet, ignore this call */

	if (!USE_PRIMENET) goto leave;

/* Open the spool file */

	fd = _open (SPOOL_FILE, _O_RDWR | _O_BINARY | _O_CREAT, 0666);
	if (fd < 0) {
		LogMsg ("ERROR: Unable to open spool file.\n");
		goto leave;
	}

/* Get the current header byte */

	header_byte = 0;
	_read (fd, &header_byte, 1);

/* If this is a maintain user info message, then set the header */
/* byte appropriately.  The 0x80 bit is the "Create a team" flag. */

	if ((msgType & 0x7F) == PRIMENET_MAINTAIN_USER_INFO) {
		header_byte |= 0x4;
		if (msgType & 0x80) header_byte |= 0x80;
		msgType = PRIMENET_MAINTAIN_USER_INFO;
	}

/* If this is a set computer info message, then set */
/* the header byte appropriately */

	else if (msgType == PRIMENET_SET_COMPUTER_INFO)
		header_byte |= 0x20;

/* If this is an update completion dates message, then set */
/* the header byte appropriately */

	else if (msgType == PRIMENET_COMPLETION_DATE)
		header_byte |= 0x8;

/* Ugly little hack when quitting GIMPS */

	else if (msgType == 999)
		header_byte |= 0x10;

/* See if this is an important message */

	else if (msgType == PRIMENET_ASSIGNMENT_RESULT)
		header_byte |= 0x2;
	else
		header_byte |= 0x1;

/* Write the new header byte */

	_lseek (fd, 0, SEEK_SET);
	_write (fd, &header_byte, 1);

/* Write out a full message */

	if (msgType == PRIMENET_RESULT_MESSAGE ||
	    msgType == PRIMENET_ASSIGNMENT_RESULT) {
		char	buf[1024];
		short	datalen;
		union {
			struct primenetResultMessage msg;
			struct primenetAssignmentResult res;
		} x;

/* Skip the remaining messages */

		while (_read (fd, buf, sizeof (buf)));

/* Append the latest message */

		if (msgType == PRIMENET_RESULT_MESSAGE) {
			memset (&x.msg, 0, sizeof (x.msg));
			if (strlen ((char *) msg) >= sizeof (x.msg.message)) {
				memcpy (x.msg.message, msg, sizeof (x.msg.message));
				x.msg.message[sizeof(x.msg.message)-1] = 0;
			} else
				strcpy (x.msg.message, (char *) msg);
			datalen = sizeof (struct primenetResultMessage);
		} else {
			datalen = sizeof (struct primenetAssignmentResult);
			memcpy (&x, msg, datalen);
		}
		_write (fd, &msgType, sizeof (short));
		_write (fd, &datalen, sizeof (short));
		_write (fd, &x, datalen);
	}

/* Set flag that the spool file has changed */

	SPOOL_FILE_CHANGED = 1;

/* Close the spool file */

	_close (fd);

/* If this is a maintain user info message, then also */
/* write a message to the results file */

leave:	if (msgType == PRIMENET_MAINTAIN_USER_INFO &&
	    (USER_NAME[0] != 0 || USER_ADDR[0] != 0)) {
		char	buf[200];
		if (USERID[0])
			sprintf (buf, "UID: %s, User: %s, %s\n",
				 USERID, USER_NAME, USER_ADDR);
		else
			sprintf (buf, "User: %s, %s\n", USER_NAME, USER_ADDR);
		writeResults (buf);
		spoolMessage (PRIMENET_RESULT_MESSAGE, buf);
		if (OLD_USERID[0]) {
			sprintf (buf, "Change UID from %s to %s\n",
				 OLD_USERID, USERID);
			writeResults (buf);
			spoolMessage (PRIMENET_RESULT_MESSAGE, buf);
		}
	}
}

/* Read a spooled message */

void readMessage (
	int	fd,
	long	*offset,	/* Offset of the message */
	short	*msgType,	/* 0 = no message */
	void	*msg)
{
	short	datalen;

/* Loop until a message that hasn't already been sent is found */

	for ( ; ; ) {
		*offset = _lseek (fd, 0, SEEK_CUR);
		if (_read (fd, msgType, sizeof (short)) != sizeof (short))
			goto err;
		if (_read (fd, &datalen, sizeof (short)) != sizeof (short))
			goto err;

/* Upgrade old style spooled messages */

		if (*msgType == PRIMENET_RESULT_MESSAGE &&
		    datalen == sizeof (struct primenet2ResultMessage)) {
			struct primenet2ResultMessage v2msg;
			struct primenetResultMessage *v3msg;
			v3msg = (struct primenetResultMessage *) msg;
			if (_read (fd, &v2msg, datalen) != datalen) goto err;
			memset (msg, 0, sizeof (struct primenetResultMessage));
			strcpy (v3msg->message, v2msg.message);
		} else if (*msgType == PRIMENET_ASSIGNMENT_RESULT &&
			   datalen==sizeof (struct primenet2AssignmentResult)){
			struct primenet2AssignmentResult v2msg;
			struct primenetAssignmentResult *v3msg;
			v3msg = (struct primenetAssignmentResult *) msg;
			if (_read (fd, &v2msg, datalen) != datalen) goto err;
			memset (msg, 0,
				sizeof (struct primenetAssignmentResult));
			v3msg->exponent = v2msg.exponent;
			v3msg->resultType = v2msg.resultType;
			memcpy (&v3msg->resultInfo, &v2msg.resultInfo,
				sizeof (v2msg.resultInfo));
		}

/* Read the body of the message */

		else
			if (_read (fd, msg, datalen) != datalen) goto err;

/* Loop if message has already been sent */

		if (*msgType != -1) break;
	}
	return;
err:	*msgType = 0;
}


/* Send a message that was read from the spool file */

int sendMessage (
	short	msgType,
	void	*msg)
{
	struct primenetResultMessage *pkt;
	struct primenetResultMessage local_pkt;
	short	structSize;
	int	return_code;
	char	buf[200];

/* Load the primenet library and do any special handling */
/* required prior to calling primenet */

	if (!LoadPrimeNet ()) return (PRIMENET_ERROR_MODEM_OFF);

/* Prepend all messages with the userid and computer id */

	pkt = (struct primenetResultMessage *) msg;
	if (msgType == PRIMENET_RESULT_MESSAGE) {
		memcpy (&local_pkt, pkt, sizeof(struct primenetResultMessage));
		if (USERID[0] == 0) 
			sprintf (local_pkt.message, "%s", pkt->message);
		else if (COMPID[0] == 0)
			sprintf (local_pkt.message, "UID: %s, %s",
				 USERID, pkt->message);
		else
			sprintf (local_pkt.message, "UID: %s/%s, %s",
				 USERID, COMPID, pkt->message);
		pkt = &local_pkt;
	}

/* Print a message on the screen and in the log file */

	switch (msgType) {
	case PRIMENET_PING_SERVER_INFO:
		OutputStr ("Contacting PrimeNet Server.\n");
		structSize = sizeof (struct primenetPingServerInfo);
		break;
	case PRIMENET_MAINTAIN_USER_INFO:
		LogMsg ("Updating user information on the server\n");
		structSize = sizeof (struct primenetUserInfo);
		break;
	case PRIMENET_SET_COMPUTER_INFO:
		LogMsg ("Updating computer information on the server\n");
		structSize = sizeof (struct primenetComputerInfo);
		break;
	case PRIMENET_GET_ASSIGNMENT:
		LogMsg ("Getting exponents from server\n");
		structSize = sizeof (struct primenetGetAssignment);
		break;
	case PRIMENET_COMPLETION_DATE:
		{
			time_t	this_time;
			char	timebuf[30];
			time (&this_time);
			this_time += ((struct primenetCompletionDate *)pkt)->days * 86400;
			strcpy (timebuf, ctime (&this_time)+4);
			strcpy (timebuf+6, timebuf+15);
			sprintf (buf, "Sending expected completion date for M%ld: %s",
				((struct primenetCompletionDate *)pkt)->exponent, timebuf);
			LogMsg (buf);
		}
		structSize = sizeof (struct primenetCompletionDate);
		break;
	case PRIMENET_RESULT_MESSAGE:
		LogMsg ("Sending text message to server:\n");
		LogMsg (pkt->message);
		structSize = sizeof (struct primenetResultMessage) - 200 +
			strlen (pkt->message) + 1;
		break;
	case PRIMENET_ASSIGNMENT_RESULT:
		if (((struct primenetAssignmentResult *)pkt)->resultType == PRIMENET_RESULT_UNRESERVE)
			sprintf (buf, "Unreserving exponent %ld\n",
				 ((struct primenetAssignmentResult *)pkt)->exponent);
		else
			sprintf (buf, "Sending result to server for exponent %ld\n",
				 ((struct primenetAssignmentResult *)pkt)->exponent);
		LogMsg (buf);
		structSize = sizeof (struct primenetAssignmentResult);
		break;
	}

/* Fill in the common header fields */

	pkt->structSize = structSize;
	if (msgType != PRIMENET_PING_SERVER_INFO) {
		pkt->versionNumber = PRIMENET_VERSION;
		strcpy (pkt->userID, USERID);
		strcpy (pkt->userPW, USER_PWD);
		strcpy (pkt->computerID, COMPID);
		pkt->hashFunction = PRIMENET_HASHFUNC_WOLTMAN_1;
	}
	hash_packet (msgType, pkt);

/* Send the message and for safety's sake re-init the FPU */

	return_code = PRIMENET (msgType, pkt);
	fpu_init ();

/* Check for the special do-not-talk-to-server-anymore return codes */

	if (return_code == PRIMENET_ERROR_KILL_CLIENT ||
	    (return_code >= PRIMENET_ERROR_SLEEP_CLIENT_1 &&
	     return_code <= PRIMENET_ERROR_SLEEP_CLIENT_200)) {
		time_t	end_time;
		time (&end_time);
		if (return_code == PRIMENET_ERROR_KILL_CLIENT)
			end_time = 1999999999L;
		else
			end_time += 3600 *
				(return_code - PRIMENET_ERROR_SLEEP_CLIENT_1 + 1);
		IniWriteInt (LOCALINI_FILE, "BlackoutEnd", end_time);
		IniWriteString (LOCALINI_FILE, "BlackoutVersion", VERSION);
		sprintf (buf, "ERROR %d: PrimeNet Blackout!!\n", return_code);
		LogMsg (buf);
		return (PRIMENET_ERROR_BLACKOUT);
	}

/* Check for a broadcast message! */

	if (return_code == PRIMENET_ERROR_OK &&
	    msgType == PRIMENET_PING_SERVER_INFO &&
	    pkt->versionNumber >= PRIMENET_PING_INFO_MESSAGE) {
		struct primenetPingServerInfo *p;
		p = (struct primenetPingServerInfo *) pkt;
		if (p->u.messageInfo.targetMap & (1L << PORT) &&
		    p->u.messageInfo.targetMap & (1L << VERSION_BIT) &&
		    p->u.messageInfo.msgID > IniGetInt (LOCALINI_FILE, "BroadcastID", 0)) {
			IniWriteInt (LOCALINI_FILE, "BroadcastID", p->u.messageInfo.msgID);
			strcat (p->u.messageInfo.szMsgText, "\n");
			LogMsg (p->u.messageInfo.szMsgText);
			writeResults (p->u.messageInfo.szMsgText);
			BroadcastMessage (p->u.messageInfo.szMsgText);
		}
	}

/* Print message based on return code */

	switch (return_code) {
	case PRIMENET_ERROR_OK:
		break;
	case 12029:
		sprintf (buf, "ERROR %d: Cannot find server.  See http://www.mersenne.org/ips/faq.html\n", return_code);
		if (msgType == PRIMENET_PING_SERVER_INFO)
			OutputStr (buf);
		else
			LogMsg (buf);
		break;
	case PRIMENET_ERROR_CONNECT_FAILED:
		sprintf (buf, "ERROR %d: Server unavailable\n", return_code);
		if (msgType == PRIMENET_PING_SERVER_INFO)
			OutputStr (buf);
		else
			LogMsg (buf);
		break;
	case PRIMENET_ERROR_SERVER_BUSY:
		sprintf (buf, "ERROR %d: Server busy\n", return_code);
		if (msgType == PRIMENET_PING_SERVER_INFO)
			OutputStr (buf);
		else
			LogMsg (buf);
		return_code = PRIMENET_ERROR_CONNECT_FAILED;
		break;
	case PRIMENET_ERROR_NO_ASSIGNMENT:
		sprintf (buf, "ERROR %d: Server has run out of exponents to assign.\n", return_code);
		LogMsg (buf);
		break;
	case PRIMENET_ERROR_ACCESS_DENIED:
		sprintf (buf, "ERROR %d: UserID/Password error.\n", return_code);
		LogMsg (buf);
		if (OLD_USERID[0]) {
			LogMsg ("Old userID and password will be used.\n");
			strcpy (USERID, OLD_USERID);
			strcpy (USER_PWD, OLD_USER_PWD);
			OLD_USERID[0] = 0;
			OLD_USER_PWD[0] = 0;
		} else {
			LogMsg ("A new userID and password will be generated.\n");
			USERID[0] = 0;
			USER_PWD[0] = 0;
		}
		spoolMessage (PRIMENET_MAINTAIN_USER_INFO, NULL);
		break;
	case PRIMENET_ERROR_ALREADY_TESTED:
		sprintf (buf, "ERROR %d: Exponent already tested.\n", return_code);
		LogMsg (buf);
		if (msgType == PRIMENET_ASSIGNMENT_RESULT) return_code = 0;
		break;
	case PRIMENET_ERROR_UNASSIGNED:
	case PRIMENET_ERROR_NOT_PERMITTED:
		sprintf (buf, "ERROR %d: Exponent not assigned to this computer.\n", return_code);
		LogMsg (buf);
		if (msgType == PRIMENET_ASSIGNMENT_RESULT) return_code = 0;
		break;
	default:
		sprintf (buf, "ERROR: Primenet error %d\n", return_code);
		LogMsg (buf);
		break;
	}

/* Print out message saying where more help may be found. */

	if (return_code)
		OutputStr ("The FAQ at http://www.mersenne.org/ips/faq.html may have more information.\n");

/* Return the return code */

	return (return_code);
}


/* Copy an existing results file to the spool file */
/* This is used when converting from manual to automatic mode */
/* It provides an extra chance that the existing results file */
/* will get sent to us. */

void spoolExistingResultsFile (void)
{
	int	i;
	char	*filename;
	char	line[256];
	FILE	*fd;

	for (i = 0; i <= 2; i++) {
		if (i == 0) filename = "results";
		if (i == 1) filename = "results.txt";
		if (i == 2) {
			if (!strcmp (RESFILE, "results")) continue;
			if (!strcmp (RESFILE, "results.txt")) continue;
			filename = RESFILE;
		}
		fd = fopen (filename, "r");
		if (fd == NULL) continue;
		while (fgets (line, sizeof (line) - 1, fd)) {
			if (line[0] == '[') continue;
			if (strstr (line, "Res64") == NULL &&
			    strstr (line, "factor:") == NULL &&
			    strstr (line, "completed P-1") == NULL &&
			    strstr (line, "no factor") == NULL) continue;
			if (line[0] == 'U' && strchr (line, ',') != NULL)
				strcpy (line, strchr (line, ',') + 2);
			spoolMessage (PRIMENET_RESULT_MESSAGE, line);
		}
		fclose (fd);
	}
}

/* Output string to screen or results file */

void OutputSomewhere (
	char	*buf)
{
	if (NO_GUI) writeResults (buf);
	else OutputStr (buf);
}

/* Output string to both the screen and results file */

void OutputBoth (
	char	*buf)
{
	OutputStr (buf);
	writeResults (buf);
}

/* Output message to screen and prime.log file */

void LogMsg (
	char	*str)
{
	int	fd;
	unsigned long filelen;
static	time_t	last_time = 0;
	time_t	this_time;

/* Output it to the screen */

	OutputStr (str);

/* Open the log file and position to the end */

	fd = _open (LOGFILE, _O_TEXT | _O_RDWR | _O_CREAT, 0666);
	if (fd < 0) {
		OutputStr ("Unable to open log file.\n");
		return;
	}
	filelen = _lseek (fd, 0L, SEEK_END);

/* Output to the log file only if it hasn't grown too big */

	if (filelen < (unsigned long) IniGetInt (INI_FILE, "MaxLogFileSize", 250000)) {

/* If it has been at least 5 minutes since the last time stamp */
/* was output, then output a new timestamp */

		time (&this_time);
		if (this_time - last_time > 300) {
			char	buf[48];
			last_time = this_time;
			buf[0] = '[';
			strcpy (buf+1, ctime (&this_time));
			sprintf (buf+25, " - ver %s]\n", VERSION);
			_write (fd, buf, strlen (buf));
		}

/* Output the message */

		_write (fd, str, strlen (str));
	}

/* Display message about full log file */
	
	else {
		char	*fullmsg = "Prime.log file full.  Please delete it.\n";
		OutputStr (fullmsg);
		if (filelen < 251000)
			_write (fd, fullmsg, strlen (fullmsg));
	}
	_close (fd);
}


/* Determine the preferred work type (based on CPU speed) */
/* Slower machines will be relegated to just factoring.  Medium speed */
/* machines will be given double-checking assignments.  Faster machines */
/* get first time LL tests. */

short default_work_type (void)
{
	double	speed;

/* Get estimated iteration time for a 512K length FFT */

	speed = map_fftlen_to_timing (524288, GW_MERSENNE_MOD, CPU_TYPE, CPU_SPEED) *
			24.0 / CPU_HOURS * 1000.0 / ROLLING_AVERAGE;

/* PIII-900 and faster get first-time LL tests */
/* P-233 and faster get double-checking work */
/* All slower machines get factoring work */

	if (speed < 0.115) return (PRIMENET_ASSIGN_TEST);
	if (speed < 0.555) return (PRIMENET_ASSIGN_DBLCHK);
	return (PRIMENET_ASSIGN_FACTOR);
}

/* Parse a line from the work-to-do file */

int parseWorkToDoLine (
	unsigned int line,	/* Line number to parse */
	struct work_unit *w)	/* Returned type of work */
{
	char	keyword[20];
	char	value[80];
	char	*comma;
		
/* Get the line from work-to-do file */

loop:	if (line > IniGetNumLines (WORKTODO_FILE)) return (FALSE);
	IniGetLineAsString (WORKTODO_FILE, line, keyword, sizeof (keyword),
			    value, sizeof (value));

/* Parse the keyword */

	if (stricmp (keyword, "Factor") == 0)
		w->work_type = WORK_FACTOR;
	else if (stricmp (keyword, "PFactor") == 0)
		w->work_type = WORK_PFACTOR;
	else if (stricmp (keyword, "Test") == 0)
		w->work_type = WORK_TEST;
	else if (stricmp (keyword, "AdvancedTest") == 0)
		w->work_type = WORK_ADVANCEDTEST;
	else if (stricmp (keyword, "DoubleCheck") == 0)
		w->work_type = WORK_DBLCHK;
	else if (stricmp (keyword, "ECM") == 0) {
		int	i;
		unsigned long j;
		char	*q;
		w->work_type = WORK_ECM;
		sscanf (value, "%ld,%lu,%lu,%ld,%ld",
			 &w->p, &w->B1, &w->B2_end, &w->curves_to_do,
			 &w->curves_completed);
		w->B2_start = w->B1;
		w->curve = 0;
		w->plus1 = 0;
		for (i = 1, q = value; i <= 5; i++, q++)
			if ((q = strchr (q, ',')) == NULL) return (TRUE);
		w->curve = atof (q);
		if ((q = strchr (q, ',')) == NULL) return (TRUE);
		w->plus1 = atoi (q+1);
		if ((q = strchr (q+1, ',')) == NULL) return (TRUE);
		j = atoi (q+1);
		if (j > w->B1) w->B2_start = j;
		return (TRUE);
	} else if (stricmp (keyword, "Pminus1") == 0) {
		w->work_type = WORK_PMINUS1;
		sscanf (value, "%ld,%lu,%lu,%ld,%ld",
			&w->p, &w->B1, &w->B2_end, &w->plus1, &w->B2_start);
		if (w->B2_start < w->B1) w->B2_start = w->B1;
		return (TRUE);
	} else if (stricmp (keyword, "AdvancedFactor") == 0) {
		w->work_type = WORK_ADVANCEDFACTOR;
		w->p = 0;
		sscanf (value, "%d,%d,%d,%d",
			&w->bits, &w->B1, &w->B2_start, &w->B2_end);
		return (TRUE);
	} else {
		char	buf[130];
		sprintf (buf, "Illegal line in INI file: %s=%s\n", keyword, value);
		OutputSomewhere (buf);
		IniDeleteLine (WORKTODO_FILE, line);
		goto loop;
	}

/* If work type involves factoring, get the number of bits factored */

	comma = strchr (value, ',');
	if (comma == NULL) w->bits = 0;
	else {
		*comma++ = 0;
		w->bits = (unsigned int) atol (comma);
		w->pminus1ed = 0;

/* And get the "has been P-1 factored" bit too. */

		comma = strchr (comma, ',');
		if (comma != NULL) {
			*comma++ = 0;
			w->pminus1ed = (int) atol (comma);
		}
	}

/* Get the exponent being tested */

	w->p = atol (value);
	return (TRUE);
}

/* Add a line of work to the work-to-do INI file.  */
/* This assumes the caller will open and close the ini file. */

void addWorkToDoLine (
	struct work_unit *w)	/* Type of work */
{
	char	buf[80];
	unsigned int i;

/* Add Advanced/Factor work to end of the file */

	if (w->work_type == WORK_ADVANCEDFACTOR) {
		sprintf (buf, "%d,%d,%d,%d", w->bits, w->B1, w->B2_start, w->B2_end);
		IniAppendLineAsString (WORKTODO_FILE, "AdvancedFactor", buf);
	}

/* Add new ECM and P-1 work after all other ECM and P-1 work, but before */
/* any other type of work */

	if (w->work_type == WORK_ECM || w->work_type == WORK_PMINUS1) {
		for (i = 1; ; i++) {
			struct work_unit temp;
			if (! parseWorkToDoLine (i, &temp)) break;
			if (temp.work_type != WORK_ECM &&
			    temp.work_type != WORK_PMINUS1) break;
		}

/* Format and insert line into ini file */

		if (w->work_type == WORK_ECM) {
			sprintf (buf, "%lu,%lu,%lu,%lu,%lu,%.0f,%lu,%lu",
				 w->p, w->B1, w->B2_end, w->curves_to_do,
				 w->curves_completed, w->curve, w->plus1,
				 w->B2_start);
			IniAppendLineAsString (WORKTODO_FILE, "ECM", buf);
		} else {
			sprintf (buf, "%lu,%lu,%lu,%lu,%lu",
				 w->p, w->B1, w->B2_end, w->plus1,
				 w->B2_start);
			IniAppendLineAsString (WORKTODO_FILE, "Pminus1", buf);
		}
		return;
	}

/* See if this exponent is already in the worktodo file. */

	for (i = 1; ; i++) {
		struct work_unit temp;

/* Read the line of the work file */

		if (! parseWorkToDoLine (i, &temp)) break;

/* Skip the line if p is not in the range */

		if (w->p != temp.p) continue;

/* Ignore this add request if it is a subset of the */
/* entry that is already in the work-to-do file */

		if (w->work_type == temp.work_type) return;
		if (w->work_type == WORK_FACTOR) return;

/* Delete the line */

		IniDeleteLine (WORKTODO_FILE, i);
		break;
	}

/* Format the number, how far factored, a P-1 factoring flag */
/* Actually, the pminus1ed flag doubles as a first vs. second LL flag */
/* in the WORK_PFACTOR case. */

	if (w->work_type == WORK_FACTOR)
		sprintf (buf, "%ld,%d", w->p, w->bits);
	else
		sprintf (buf, "%ld,%d,%d", w->p, w->bits, w->pminus1ed);

/* Append a line */

	if (w->work_type == WORK_FACTOR)
		IniAppendLineAsString (WORKTODO_FILE, "Factor", buf);
	if (w->work_type == WORK_PFACTOR)
		IniAppendLineAsString (WORKTODO_FILE, "PFactor", buf);
	if (w->work_type == WORK_TEST)
		IniAppendLineAsString (WORKTODO_FILE, "Test", buf);
	if (w->work_type == WORK_DBLCHK)
		IniAppendLineAsString (WORKTODO_FILE, "DoubleCheck", buf);
}

/* Open the results file and see if a given exponent has */
/* already been worked on. */ 

void checkResultsFile (
	unsigned long p,
	int	*tested,
	int	*factored)
{
	int	fd;
	char	buf[80];
	long	offset = 0;

/* Assume prime has not been lucas tested or factored. */

	*tested = FALSE;
	*factored = FALSE;

/* Loop until entire file has been processed */

	fd = _open (RESFILE, _O_BINARY | _O_RDONLY);
	if (fd < 0) return;
	for ( ; ; ) {
		unsigned long file_p = 0;
		int	i = 0;
		_lseek (fd, offset, SEEK_SET);
		buf[_read (fd, buf, sizeof (buf) - 1)] = 0;
		while (buf[i]) if (buf[i++] == 'M') break;
		if (buf[i] == 0) break;
		for (file_p = 0; isdigit(buf[i]); i++)
			file_p = file_p * 10 + buf[i] - '0';
		while (isspace(buf[i])) i++;
		if (file_p == p) {
			if (buf[i] == 'h') *tested = TRUE;
			if (buf[i] == 'n') *factored = TRUE;
			if (buf[i] == 'i') *tested = TRUE;
		}
		while (buf[i] && buf[i] != '\n') i++;
		offset += i;
	}
	_close (fd);
}

/* Return the number of seconds until vacation is over */

unsigned long secondsUntilVacationEnds (void)
{
	time_t	current_time;

	if (VACATION_END == 0) return (0);
	time (&current_time);
	if (current_time >= VACATION_END) {
		VACATION_END = 0;
		IniWriteInt (LOCALINI_FILE, "VacationEnd", 0);
		return (0);
	}
	return (VACATION_END - current_time);
}

/* Return a double estimating current exponent's percentage complete. */

double pct_complete (
	int	work_type,
	unsigned long p,
	unsigned long *iteration)
{
	unsigned long counter;
	short	type;
	int	fd;
	char	filename[20];

/* If we are working on this exponent now, then global variables are */
/* set every now and then which tell us our percent completion. */

	if (p == EXP_BEING_WORKED_ON) {
		if (EXP_PERCENT_COMPLETE == 999.0) {
			*iteration = p;
			return (0.999);
		}
		if (work_type == WORK_FACTOR && EXP_BEING_FACTORED) {
			*iteration = (unsigned long) (EXP_PERCENT_COMPLETE * 16 + 1);
			return (EXP_PERCENT_COMPLETE);
		}
		if (work_type != WORK_FACTOR && !EXP_BEING_FACTORED) {
			*iteration = (unsigned long) (EXP_PERCENT_COMPLETE * p);
			return (EXP_PERCENT_COMPLETE);
		}
	}

/* Otherwise, lets see if there is a intermediate file.  If there is */
/* read it to try and figure out our percent completion */

	tempFileName (filename, p);
	if (fileExists (filename) &&
	    readFileHeader (filename, &fd, &type, &counter)) {

/* A LL test save file.  The counter field will tell us how much */
/* is already done. */

		if (type >= 4) {
			_close (fd);
			*iteration = counter;
			return ((double) counter / (double) p);
		}

/* A factoring intermediate file */

		if (type == 2) {
			short	pass, result, old_style, bits, limit;
			double	pct;

/* If this is an LL assignment, then the factoring work is pretty */
/* much negligible. */

			if (work_type != WORK_FACTOR) {
				_close (fd);
				*iteration = 1;
				return (0.001);
			}

/* Read the rest of the intermediate file and ignore any old-style */
/* factoring files. */

			if (_read (fd, &pass, sizeof (short)) != sizeof (short))
				goto no_est;
			if (pass != 999)
				goto no_est;
			if (_read (fd, &result, sizeof (short)) != sizeof (short))
				goto no_est;
			if (_read (fd, &old_style, sizeof (short)) != sizeof (short))
				goto no_est;
			if (old_style)
				goto no_est;
			if (_read (fd, &bits, sizeof (short)) != sizeof (short))
				goto no_est;
			if (_read (fd, &pass, sizeof (short)) != sizeof (short))
				goto no_est;
			_close (fd);

/* Now return the percent complete.  We assume that the how_far_factored */
/* value in worktodo.ini is not close to the factorLimit value. */

			*iteration = pass;
			if (result) return ((double) pass / 16.0);
			limit = factorLimit (p, work_type);
			pct = pow (2.0, bits - limit);
			return (pct * (1.0 + (double) pass / 16.0));

/* An error or intermediate file that we do not recognize. */
/* Return a special code indicating an unknown percent complete. */

no_est:			_close (fd);
			*iteration = 1;
			return (999.0);
		}
	}

/* We have not started working on this item */

	*iteration = 0;
	return (0.0);
}

/* Given an exponent, return the fft length from the local.ini file. */
/* Used in implementing the "soft" FFT crossover points in version 22.8 */

unsigned long fftlen_from_ini_file (
	unsigned long p)
{
	char	buf[80];

	IniGetString (LOCALINI_FILE, "SoftCrossoverData",
		      buf, sizeof (buf), "0");
	if (p == (unsigned long) atol (buf)) {
		char	*comma;
		unsigned long fftlen;
		comma = strchr (buf, ',');
		if (comma != NULL) {
			*comma++ = 0;
			fftlen = atol (comma);
			comma = strchr (comma, ',');
			if (comma != NULL) {
				unsigned long sse2;
				*comma++ = 0;
				sse2 = atol (comma);
				if ((sse2 && (CPU_FLAGS & CPU_SSE2)) ||
				    (!sse2 && !(CPU_FLAGS & CPU_SSE2)))
					return (fftlen);
			}
		}
	}
	return (0);
}

/* Given an exponent, determine the fft length.  Handles the "soft" */
/* FFT crossover points implemented in version 22.8 */

unsigned long advanced_map_exponent_to_fftlen (
	unsigned long p)
{
	unsigned long fftlen = fftlen_from_ini_file (p);
	return (fftlen ? fftlen : map_exponent_to_fftlen (p, GW_MERSENNE_MOD));
}

/* Make a guess as to how long a chore should take. */

double raw_work_estimate (
	struct work_unit *w)
{
	double	timing;
	double	est;

/* Not very accurate estimates for ECM and P-1! */

	if (w->work_type == WORK_ECM ||
	    w->work_type == WORK_PMINUS1 ||
	    w->work_type == WORK_ADVANCEDFACTOR)
		est = 86400.0;

/* If factoring, guess how long that will take.  Timings are based on */
/* how long it takes my PII-400 to process the exponent 12,000,017. */
/* Below 2^60, prime95 runs through 0.004093*2^58 factors in 3.198 seconds. */
/* Below 2^62, prime95 runs through 0.004093*2^58 factors in 3.204 seconds. */
/* Below 2^64, prime95 runs through 0.004093*2^58 factors in 5.949 seconds. */
/* Above 2^64, prime95 runs through 0.004093*2^58 factors in 13.511 seconds. */
/* Compute timing * 2^limit / (0.004093 * 2^58) * (12,000,017 / p) */
/* Which simplifies to: timing * 2^(limit-44) * 178945.25 / p */

	if (w->work_type == WORK_FACTOR) {
		unsigned int limit;
		limit = factorLimit (w->p, w->work_type);
		timing = (limit > 64) ? 13.511 :
			 (limit > 62) ? 5.949 :
			 (limit > 60) ? 3.204 : 3.198;
		est = timing * 178945.25 * (1L << (limit - 44)) / w->p;
		est = est * 400.0 / CPU_SPEED;
		if (CPU_TYPE <= 4) est *= REL_486_SPEED;
		if (CPU_TYPE == 7) est *= REL_K6_SPEED;
		if (CPU_TYPE == 11) est *= REL_K7_SPEED;
		if (CPU_TYPE == 12) est *= REL_P4_SPEED;
		if (w->bits >= limit) est = 0.0;
		else if (w->bits == limit-1) est *= 0.5;
		else if (w->bits == limit-2) est *= 0.75;
		else if (w->bits == limit-3) est *= 0.875;
	}

/* If P-1 factoring, then estimate the time.  Remember the pminus1ed */
/* flag doubles as a first vs. second LL flag. */

	if (w->work_type == WORK_PFACTOR) {
		unsigned int bits;
		unsigned long bound1, bound2, squarings;
		double	prob;
		bits = factorLimit (w->p, w->work_type);
		if (w->bits > bits) bits = w->bits;
		guess_pminus1_bounds (w->p, bits, w->pminus1ed,
				      &bound1, &bound2, &squarings, &prob);
		est = squarings * map_fftlen_to_timing (
			advanced_map_exponent_to_fftlen (w->p),
			GW_MERSENNE_MOD, CPU_TYPE, CPU_SPEED);
	}

/* If testing add in the Lucas-Lehmer testing time */

	if (w->work_type == WORK_TEST ||
	    w->work_type == WORK_ADVANCEDTEST ||
	    w->work_type == WORK_DBLCHK) {
		est = w->p * map_fftlen_to_timing (
			advanced_map_exponent_to_fftlen (w->p),
			GW_MERSENNE_MOD, CPU_TYPE, CPU_SPEED);
	}

/* Return the total estimated time in seconds */

	return (est);
}


/* Make a guess as to how long a chore should take. */

double work_estimate (
	struct work_unit *w)
{
	return (raw_work_estimate (w) *
		24.0 / CPU_HOURS * 1000.0 / ROLLING_AVERAGE);
}


/* Determine how much we should factor (in bits) */
/* Don't let feeble 486 any Cyrix machines factor higher than 2^62 */
/* as this would execute some really slow floating point code */

unsigned int factorLimit (
	unsigned long p,
	int	work_type)
{
	unsigned int test, override;

	if (p > FAC72) test = 72;	/* Test all 72 bit factors */
	else if (p > FAC71) test = 71;	/* Test all 71 bit factors */
	else if (p > FAC70) test = 70;	/* Test all 70 bit factors */
	else if (p > FAC69) test = 69;	/* Test all 69 bit factors */
	else if (p > FAC68) test = 68;	/* Test all 68 bit factors */
	else if (p > FAC67) test = 67;	/* Test all 67 bit factors */
	else if (p > FAC66) test = 66;	/* Test all 66 bit factors */
	else if (p > FAC65) test = 65;	/* Test all 65 bit factors */
	else if (p > FAC64) test = 64;	/* Test all 64 bit factors */
	else if (p > FAC63) test = 63;	/* Test all 63 bit factors */
	else if (p > FAC62) test = 62;	/* Test all 62 bit factors */
	else if (p > FAC61) test = 61;	/* Test all 61 bit factors */
	else if (p > FAC60) test = 60;	/* Test all 60 bit factors */
	else if (p > FAC59) test = 59;	/* Test all 59 bit factors */
	else if (p > FAC58) test = 58;	/* Test all 58 bit factors */
	else if (p > FAC57) test = 57;	/* Test all 57 bit factors */
	else if (p > FAC56) test = 56;	/* Test all 56 bit factors */
	else test = 40;			/* Test all 40 bit factors */
	if (work_type == WORK_FACTOR && !USE_PRIMENET) {
		override = (unsigned int)
			IniGetInt (INI_FILE, "FactorOverride", 99);
		if (override != 99) test = override;
	}
	if (work_type == WORK_DBLCHK) test--;
	return (test);
}

/* Generate temporary file name */

void tempFileName (
	char	*buf,
	unsigned long p)
{
	char	c;
	if (p == 0) {
		IniGetString (INI_FILE, "AdvFacFileName", buf, 10, "p0000000");
		return;
	}
	sprintf (buf, "p%07li", p % 10000000);
	if (p >= 10000000) buf[1] = (char) ('A' + (p / 1000000) - 10);
	if (p >= 35000000) c=buf[1], buf[1]=buf[2], buf[2]=(char)(c - 25);
	if (p >= 70000000) c=buf[2], buf[2]=buf[3], buf[3]=(char)(c - 25);
}

/* See if the given file exists */

int fileExists (
	char	*filename)
{
	int	fd;
	fd = _open (filename, _O_RDONLY | _O_BINARY);
	if (fd < 0) return (0);
	_close (fd);
	return (1);
}

/* Read the header of an intermediate Lucas-Lehmer results file */

int readFileHeader (
	char	*filename,
	int	*fd,
	short	*type,
	unsigned long *counter)
{
	*fd = _open (filename, _O_BINARY | _O_RDONLY);
	if (*fd <= 0) return (FALSE);
	if (_read (*fd, type, sizeof (short)) != sizeof (short)) goto readerr;
	if (_read (*fd, counter, sizeof (long)) != sizeof (long)) goto readerr;
	return (TRUE);
readerr:_close (*fd);
	return (FALSE);
}

/* Open the results file and write a line to the end of it. */

int writeResults (
	char	*msg)
{
static	time_t	last_time = 0;
	time_t	this_time;
	int	fd;

/* Open file, position to end */

	fd = _open (RESFILE, _O_TEXT | _O_RDWR | _O_CREAT | _O_APPEND, 0666);
	if (fd < 0) {
		LogMsg ("Error opening results file to output this message:\n");
		LogMsg (msg);
		return (FALSE);
	}

/* If it has been at least 5 minutes since the last time stamp */
/* was output, then output a new timestamp */

	time (&this_time);
	if (this_time - last_time > 300) {
		char	buf[32];
		last_time = this_time;
		buf[0] = '[';
		strcpy (buf+1, ctime (&this_time));
		buf[25] = ']';
		buf[26] = '\n';
		_write (fd, buf, 27);
	}

/* Output a USERID prefix for result messages */

	if ((msg[0] == 'M' || msg[0] == 'P') &&
	    isdigit (msg[1]) && USERID[0]) {
		if (_write (fd, "UID: ", 5) < 0) goto fail;
		if (_write (fd, USERID, strlen (USERID)) < 0) goto fail;
		if (COMPID[0]) {
			if (_write (fd, "/", 1) < 0) goto fail;
			if (_write (fd, COMPID, strlen(COMPID)) < 0) goto fail;
		}
		if (_write (fd, ", ", 2) < 0) goto fail;
	}

/* Output the message */

	if (_write (fd, msg, strlen (msg)) < 0) goto fail;
	_close (fd);
	return (TRUE);

/* On a write error, close file and return error flag */

fail:	_close (fd);
	LogMsg ("Error writing message to results file:\n");
	LogMsg (msg);
	return (FALSE);
}


/* Unreserve an exponent */

void unreserve (
	unsigned long p)
{
	struct primenetAssignmentResult pkt;
	unsigned int i;

/* Build a packet and spool message */

	memset (&pkt, 0, sizeof (pkt));
	pkt.exponent = p;
	pkt.resultType = PRIMENET_RESULT_UNRESERVE;
	spoolMessage (PRIMENET_ASSIGNMENT_RESULT, &pkt);

/* Find exponent in worktodo.ini and delete it if present */

	IniFileOpen (WORKTODO_FILE, 0);
	for (i = 1; ; i++) {
		struct work_unit temp;

/* Read the line of the work file */

		if (! parseWorkToDoLine (i, &temp)) break;

/* Skip the line if exponent is not a match */

		if (p != temp.p) continue;

/* Delete the line */

		IniDeleteLine (WORKTODO_FILE, i);
		break;
	}
	IniFileClose (WORKTODO_FILE);
}

/* Send any queued up messages to the server.  See if we have enough */
/* work queued up.  If we have too much work, give some back */

int communicateWithServer (void)
{
static	int	msgs_to_send = 0;/* 0 = no messages, 0x1 = messages */
				/* 0x2 = critical messages */
	char	header_byte;	/* Cached flag byte from spool file */
				/* 0 = no spool file */
				/* 0x1 = informational msgs in file */
				/* 0x2 = important msgs in file */
				/* 0x4 = user info has changed */
				/* 0x8 = completion dates need updating */
				/* 0x10 = return all work */
				/* 0x20 = computer info has changed */
	int	fd;		/* Spool file handle */
	time_t	this_time, blackout_end_time;
	unsigned int i;
	int	have_enough_work, retries, est_is_accurate;
	double	est, work_to_get, unreserve_threshold;
	unsigned long vacation_time;
	int	rc;
	int	talked_to_server = 0;
	int	manual_comm = 0;
static	int	work_queue_counter = 0;

/* Use this opportunity to perform other miscellaneous tasks that may */
/* be required by this particular OS port */

	doMiscTasks ();

/* If we are operating manually, return */

	if (!USE_PRIMENET) return (TRUE);

/* Every now and then (remember this routine is called frequently), see */
/* if we need to update server completion dates or check the work queue. */

	if (++work_queue_counter > 1000) {
		ConditionallyUpdateEndDates ();
		CHECK_WORK_QUEUE = 1;
		work_queue_counter = 0;
	}

/* If all communication with the server is manually initiated, */
/* then handle that case here */

	if (MANUAL_COMM & 0x2) {
		MANUAL_COMM &= 0xFD;
		manual_comm = 1;
		next_comm_time = 0;
		talked_to_server = 1;
	} else if (MANUAL_COMM & 0x1) {
		return (TRUE);
	}

/* If a computer ID has not been given, then generate one here */

	if (COMPID[0] == 0) {
		unsigned long id, hi, lo;
		srand ((unsigned) time (NULL));
		id = ((unsigned long) rand () << 16) + rand ();
		if (CPU_FLAGS & CPU_RDTSC) { rdtsc (&hi,&lo); id += lo; }
		sprintf (COMPID, "C%08X", id);
		IniWriteString (LOCALINI_FILE, "ComputerID", COMPID);
		spoolMessage (PRIMENET_SET_COMPUTER_INFO, NULL);
	}

/* Open and read the spool file to set msgs_to_send correctly */

	if (SPOOL_FILE_CHANGED) {
		SPOOL_FILE_CHANGED = 0;
		fd = _open (SPOOL_FILE, _O_RDONLY | _O_BINARY);
		if (fd >= 0) {

/* Read the first byte.  This indicates if there are */
/* "important" messages queued up.  Important messages */
/* will cause us to try and contact the server every */
/* few minutes until we get through. */

			_read (fd, &header_byte, 1);
			if (header_byte > 1) msgs_to_send |= 0x1;
			if (header_byte & 0x10) msgs_to_send |= 0x2;
			_close (fd);
		}
	}

/* If no important messages queued up AND we're not */
/* checking to see if we have enough work to do, then return */

	if (!msgs_to_send && !CHECK_WORK_QUEUE) return (TRUE);

/* Make sure MODEM_RETRY_TIME or NETWORK_RETRY_TIME minutes have */
/* passed since our last attempt to communicate with the server. */

loop:	time (&this_time);
#ifndef SERVER_TESTING
	if (this_time < next_comm_time) goto exit_or_loop_nomsg;
#endif

/* Make sure we don't pummel the server with data.  Suppose user uses */
/* Advanced/Factor to find tiny factors again.  At least, make sure */
/* he only sends the data once every 5 minutes. */

	next_comm_time = this_time + 300;

/* If server has requested that we do not contact it for the next N hours */
/* then obey the server's request. */

	blackout_end_time = IniGetInt (LOCALINI_FILE, "BlackoutEnd", 0);
	if (blackout_end_time) {
		char	blackout_version[9];
		IniGetString (LOCALINI_FILE, "BlackoutVersion",
			      (char *) blackout_version, 8, "99.9");
		if (this_time < blackout_end_time &&
		    strcmp (VERSION, blackout_version) <= 0) {
			OutputStr (blackout_end_time == 1999999999L ? BLACKMSG2 : BLACKMSG1);
			goto exit_or_loop_nomsg;
		}
		IniWriteString (LOCALINI_FILE, "BlackoutEnd", NULL);
		IniWriteString (LOCALINI_FILE, "BlackoutVersion", NULL);
	}

/* If we have messages to send, see if the server is available */
/* Use Ping so that we don't write to the log file every few minutes */
/* when the server is unavailable */

	if (msgs_to_send) {
		struct primenetPingServerInfo pkt;
		memset (&pkt, 0, sizeof (pkt));
		pkt.u.serverInfo.versionNumber = PRIMENET_PING_INFO_MESSAGE + PRIMENET_VERSION;
		rc = sendMessage (PRIMENET_PING_SERVER_INFO, &pkt);
		if (rc) goto exit_or_loop;
		talked_to_server = 1;
	}

/* Send all the queued up messages */

	header_byte = 0;
	fd = _open (SPOOL_FILE, _O_RDWR | _O_BINARY);
	if (fd >= 0) {

/* Read the first byte.  This indicates if there are */
/* "important" messages queued up.  Important messages */
/* will cause us to try and contact the server every */
/* few minutes until we get through. */

		_read (fd, &header_byte, 1);

/* Send any changes to the user info first */

		if (header_byte & 0x4) {
			struct primenetUserInfo pkt;

			memset (&pkt, 0, sizeof (pkt));
			strcpy (pkt.userName, USER_NAME);
			strcpy (pkt.userEmailAddr, USER_ADDR);
			strcpy (pkt.oldUserID, OLD_USERID);
			strcpy (pkt.oldUserPW, OLD_USER_PWD);
			if (NEWSLETTERS)
				pkt.bUserOptions |= PRIMENET_USER_OPTION_SENDEMAIL;
			if (header_byte & 0x80)
				pkt.bUserOptions |= PRIMENET_USER_OPTION_TEAMACCT;
			rc = sendMessage (PRIMENET_MAINTAIN_USER_INFO, &pkt);
			if (rc) {
				_close (fd);
				goto exit_or_loop;
			}
			talked_to_server = 1;

/* Delete old User id info */

			OLD_USERID[0] = 0;
			OLD_USER_PWD[0] = 0;
			IniWriteString (INI_FILE, "OldUserID", OLD_USERID);
			IniWriteString (INI_FILE, "OldUserPWD", OLD_USER_PWD);

/* If a user ID has just been assigned to us, then write it the INI file */

			if (strcmp (USERID, pkt.userID)) {
				strcpy (USERID, pkt.userID);
				IniWriteString (INI_FILE, "UserID", USERID);
			}
			if (strcmp (USER_PWD, pkt.userPW)) {
				strcpy (USER_PWD, pkt.userPW);
				IniWriteString (INI_FILE, "UserPWD", USER_PWD);
			}
		}

/* Send any changes to the computer info next. */
/* Do this whenever the computer info changes, the user info changes */
/* or at Scott's request whenever new end dates are sent. */

		if (header_byte & 0x2C) {
			struct primenetComputerInfo pkt;

			memset (&pkt, 0, sizeof (pkt));
			pkt.cpu_type = CPU_TYPE;
			pkt.speed = (short) CPU_SPEED;
			pkt.hours_per_day = CPU_HOURS;
			rc = sendMessage (PRIMENET_SET_COMPUTER_INFO, &pkt);
			if (rc) {
				_close (fd);
				goto exit_or_loop;
			}
			talked_to_server = 1;

/* Flag the message as successfully sent */

			header_byte &= 0xDB;
			_lseek (fd, 0, SEEK_SET);
			_write (fd, &header_byte, 1);
		}

/* Send the messages */

		for ( ; ; ) {
			short	msgType;
			union {
				struct primenetResultMessage msg;
				struct primenetAssignmentResult res;
			} msg;
			long	old_offset, new_offset;

			readMessage (fd, &old_offset, &msgType, &msg);
			if (msgType == 0) break;
			rc = sendMessage (msgType, &msg);
			if (rc) {
				_close (fd);
				goto exit_or_loop;
			}
			talked_to_server = 1;

/* Flag the message as successfully sent */

			new_offset = _lseek (fd, 0, SEEK_CUR);
			_lseek (fd, old_offset, SEEK_SET);
			msgType = -1;
			_write (fd, &msgType, sizeof (short));
			_lseek (fd, new_offset, SEEK_SET);

/* Stop sending to server if so requested */

			if (escapeCheck ()) {
				_close (fd);
				return (FALSE);
			}
		}
		header_byte &= 0x7C;

/* Close the spool file */

		_close (fd);
	}

/* At this point, we can usually delete the spool file */

	if (header_byte == 0) _unlink (SPOOL_FILE);

/* Clear flag indicating work queue needs checking */
/* We will set msgs_to_send if we find we need to return */
/* or request exponents.  NOTE: We must be extra careful here */
/* to avoid opening the work-to-do file every few minutes when */
/* the server is unavailable. */

	CHECK_WORK_QUEUE = 0;

/* Get work to do until we've accumulated enough to keep */
/* us busy for a while (DAYS_OF_WORK or the length of a vacation). */
/* If we have too much work to do, lets give some back. */

	have_enough_work = 0;
	unreserve_threshold = IniGetInt (INI_FILE, "UnreserveDays", 30) * 86400.0;
	est = 0.0; est_is_accurate = TRUE;
	vacation_time = secondsUntilVacationEnds ();
	work_to_get = DAYS_OF_WORK * 86400.0;
	if (ON_DURING_VACATION && work_to_get < vacation_time)
		work_to_get = vacation_time;
	IniFileOpen (WORKTODO_FILE, 0);
	for (i = 1; ; i++) {
		struct work_unit w;
		unsigned long iteration;
		double	work, pct;

/* Read the line of the work file */

		if (! parseWorkToDoLine (i, &w)) break;

/* Adjust our time estimate */

		work = work_estimate (&w);
		pct = pct_complete (w.work_type, w.p, &iteration);
		if (pct == 999.0) est_is_accurate = FALSE;
		else work *= (1.0 - pct);
		est += work;

/* Don't trouble server with work types that it does not keep track of */

		if (w.work_type != WORK_FACTOR &&
		    w.work_type != WORK_PFACTOR &&
		    w.work_type != WORK_TEST &&
		    w.work_type != WORK_DBLCHK);

/* If we have too much work queued, give it back */

		else if (header_byte & 0x10 ||
			 (have_enough_work &&
			  est_is_accurate &&
			  est - work >= work_to_get + unreserve_threshold &&
			  pct == 0.0)) {
			struct primenetAssignmentResult pkt;
			memset (&pkt, 0, sizeof (pkt));
			pkt.exponent = w.p;
			pkt.resultType = PRIMENET_RESULT_UNRESERVE;
			msgs_to_send |= 0x1;
			rc = sendMessage (PRIMENET_ASSIGNMENT_RESULT, &pkt);
			if (rc) {
				IniFileClose (WORKTODO_FILE);
				goto exit_or_loop;
			}
			talked_to_server = 1;
			IniDeleteLine (WORKTODO_FILE, i);
			i--;
		}

/* Acknowledge the exponent by sending a projected completion date */
/* If an unexpected error occurs, ignore it.  We can get "not-assigned" */
/* errors when testing manually reserved exponents and we've just */
/* enabled primenet. */

		else if (header_byte & 0x8) {
			struct primenetCompletionDate pkt2;
			double	comp_date;
			memset (&pkt2, 0, sizeof (pkt2));
			pkt2.exponent = w.p;
			if (!ON_DURING_VACATION)
				comp_date = est + vacation_time;
			else if (est < vacation_time)
				comp_date = vacation_time;
			else
				comp_date = est;
			pkt2.days = (unsigned long) (comp_date / 86400.0) + 1;
			pkt2.requestType =
				(w.work_type == WORK_FACTOR) ? PRIMENET_ASSIGN_FACTOR :
				(w.work_type == WORK_PFACTOR) ? PRIMENET_ASSIGN_PFACTOR :
				(w.work_type == WORK_TEST) ? PRIMENET_ASSIGN_TEST :
					PRIMENET_ASSIGN_DBLCHK;
			pkt2.programType = PRIMENET_PROGRAM_WOLTMAN;
			pkt2.iteration = iteration;
			if (ON_DURING_VACATION &&
			    vacation_time > (unsigned long)
						DAYS_BETWEEN_CHECKINS * 86400)
				pkt2.nextMsg = vacation_time / 86400;
			else
				pkt2.nextMsg = DAYS_BETWEEN_CHECKINS;
			rc = sendMessage (PRIMENET_COMPLETION_DATE, &pkt2);
			if (rc == PRIMENET_ERROR_CONNECT_FAILED ||
			    rc == PRIMENET_ERROR_ACCESS_DENIED) {
				IniFileClose (WORKTODO_FILE);
				goto exit_or_loop;
			}
			if (i > 1 &&
			    (rc == PRIMENET_ERROR_ALREADY_TESTED ||
			     rc == PRIMENET_ERROR_NOT_PERMITTED)) {
				est -= work;
				IniDeleteLine (WORKTODO_FILE, i);
				i--;
			}
			talked_to_server = 1;
		}

/* Set flag if we have enough work */

		if (est >= work_to_get)
			have_enough_work = 1;
	}

/* Delete the spool file. */

	if (header_byte) _unlink (SPOOL_FILE);

/* If we just returned all work to the server, enter manual mode */
/* and exit - the user is quitting GIMPS */

	if (header_byte & 0x10 ||
	    (IniGetInt (INI_FILE, "NoMoreWork", 0) &&
	     IniGetNumLines (WORKTODO_FILE) == 0)) {
		USE_PRIMENET = 0;
		IniWriteInt (INI_FILE, "UsePrimenet", 0);
		IniWriteInt (INI_FILE, "NoMoreWork", 0);
		GIMPS_QUIT = 1;
		OutputSomewhere ("Successfully quit GIMPS.\n");
		goto exit;
	}

/* After sending new completion dates remember the current time */
/* so that we can send new completion dates in a month. */

	if (header_byte & 0x8) {
		time_t current_time;
		time (&current_time);
		IniWriteInt (LOCALINI_FILE, "LastEndDatesSent", current_time);
	}

/* Before we get work from the server, make absolutely sure that we */
/* have write access to the worktodo.ini file. */

	if (!have_enough_work && !IniFileWritable (WORKTODO_FILE)) {
		char	buf[80];
		sprintf (buf, "Cannot write to file %s.\n", WORKTODO_FILE);
		OutputBoth (buf);
		IniFileClose (WORKTODO_FILE);
		goto exit_or_loop;
	}

/* If we don't have enough work to do, get more work from the server. */

	retries = 0;
	while (!have_enough_work &&
	       est_is_accurate &&
	       ! IniGetInt (INI_FILE, "NoMoreWork", 0) &&
	       IniGetNumLines (WORKTODO_FILE) < (unsigned int)
			IniGetInt (INI_FILE, "MaxExponents", 20)) {
		struct primenetGetAssignment pkt1;
		struct primenetCompletionDate pkt2;
		struct work_unit w;
		double	new_est, comp_date;

/* Check for ESC */

		if (escapeCheck ()) {
			IniFileClose (WORKTODO_FILE);
			return (FALSE);
		}

/* Get an exponent to work on */

		memset (&pkt1, 0, sizeof (pkt1));
		if (WORK_PREFERENCE)
			pkt1.requestType = WORK_PREFERENCE;
		else
			pkt1.requestType = default_work_type ();
		if (pkt1.requestType & PRIMENET_ASSIGN_BIGONES)
			pkt1.requestType =
				PRIMENET_ASSIGN_BIGONES | PRIMENET_ASSIGN_TEST;
		pkt1.programType = PRIMENET_PROGRAM_WOLTMAN;
		pkt1.how_far_factored = (CPU_TYPE < 5) ? 486.0 : 64.0;
		msgs_to_send |= 0x1;
		rc = sendMessage (PRIMENET_GET_ASSIGNMENT, &pkt1);
		if (rc) {
			IniFileClose (WORKTODO_FILE);
			goto exit_or_loop;
		}
		talked_to_server = 1;

/* For some strange reason the server has been known to send us */
/* bogus exponents (zero) to test */

		if (pkt1.exponent < 1000000) {
			char	buf[80];
			sprintf (buf, "Server sent bad exponent: %ld.\n",
				 pkt1.exponent);
			LogMsg (buf);
			if (++retries < 3) continue;
			IniFileClose (WORKTODO_FILE);
			goto exit_or_loop;
		}

/* Add the exponent to our time estimate */

		w.work_type = (pkt1.requestType == PRIMENET_ASSIGN_FACTOR) ?
						WORK_FACTOR :
			      (pkt1.requestType == PRIMENET_ASSIGN_PFACTOR) ?
						WORK_PFACTOR :
			      (pkt1.requestType == PRIMENET_ASSIGN_TEST) ?
						WORK_TEST : WORK_DBLCHK;
		w.p = pkt1.exponent;
		w.bits = (unsigned int) pkt1.how_far_factored;
		/* Kludge that uses 0.5 to indicate P-1 factoring has been done */
		w.pminus1ed = (pkt1.how_far_factored - (double) w.bits) == 0.5;
		new_est = est + work_estimate (&w);

/* Acknowledge the exponent by sending a projected completion date */

		memset (&pkt2, 0, sizeof (pkt2));
		pkt2.exponent = pkt1.exponent;
		if (!ON_DURING_VACATION)
			comp_date = new_est + vacation_time;
		else if (new_est < vacation_time)
			comp_date = vacation_time;
		else
			comp_date = new_est;
		pkt2.days = (unsigned long) (comp_date / 86400.0) + 1;
		pkt2.requestType =
			(w.work_type == WORK_FACTOR) ? PRIMENET_ASSIGN_FACTOR :
			(w.work_type == WORK_PFACTOR) ? PRIMENET_ASSIGN_PFACTOR :
			(w.work_type == WORK_TEST) ? PRIMENET_ASSIGN_TEST :
				PRIMENET_ASSIGN_DBLCHK;
		pkt2.programType = PRIMENET_PROGRAM_WOLTMAN;
		pkt2.iteration = 0;
		if (ON_DURING_VACATION &&
		    vacation_time > (unsigned long)
					DAYS_BETWEEN_CHECKINS * 86400)
			pkt2.nextMsg = vacation_time / 86400;
		else
			pkt2.nextMsg = DAYS_BETWEEN_CHECKINS;
		rc = sendMessage (PRIMENET_COMPLETION_DATE, &pkt2);
		if (rc) {
			if (rc != PRIMENET_ERROR_BLACKOUT && ++retries < 3)
				continue;
			IniFileClose (WORKTODO_FILE);
			goto exit_or_loop;
		}
		talked_to_server = 1;

/* Once the server accepts our acknowledgement, write the exponent */
/* to our worktodo file */

		addWorkToDoLine (&w);

/* Set flag if we have enough work */

		est = new_est;
		if (est >= work_to_get)
			have_enough_work = 1;
	}

/* Close the work-to-do file */

exit:	IniFileClose (WORKTODO_FILE);

/* Tell user we're done communicating */

	if (talked_to_server)
		OutputStr ("Done communicating with server.\n");

/* Clear flag indicating msgs need to be sent to primenet server */

	msgs_to_send = 0;
	return (TRUE);

/* We could not contact server, either exit or loop */
/* We loop when we are quitting GIMPS */

exit_or_loop:
	if (!manual_comm) {
		char	buf[80];
		unsigned int retry_time;
		retry_time = (rc == PRIMENET_ERROR_MODEM_OFF) ?
			MODEM_RETRY_TIME : NETWORK_RETRY_TIME;
		sprintf (buf, "Will try contacting server again in %d %s.\n",
			 retry_time, retry_time == 1 ? "minute" : "minutes");
		OutputStr (buf);
		next_comm_time = this_time + 60 * retry_time;
	}
exit_or_loop_nomsg:
	if (msgs_to_send & 0x2) {
		Sleep (5000);
		if (escapeCheck ()) return (FALSE);
		next_comm_time = 0;
		goto loop;
	}
	return (TRUE);
}

/**************************************************************/
/*             The All New (v20) P-1 factoring stage!         */
/**************************************************************/

/* This table gives the values of Dickman's function given an input */
/* between 0.000 and 0.500.  These values came from a different program */
/* that did a numerical integration. */

static double savedF[501] = {
	0, 0, 0, 0, 0, 0, 3.3513e-215, 5.63754e-208, 4.00865e-201,
	1.65407e-194, 4.53598e-188, 8.93587e-182, 1.33115e-175,
	1.55557e-169, 1.46609e-163, 1.13896e-157, 7.42296e-152,
	3.80812e-146, 1.56963e-140, 5.32886e-135, 1.51923e-129,
	3.69424e-124, 7.76066e-119, 1.42371e-113, 2.30187e-108,
	3.30619e-103, 4.24793e-098, 4.80671e-093, 4.78516e-088,
	4.22768e-083, 3.33979e-078, 2.37455e-073, 1.52822e-068,
	8.94846e-064, 4.78909e-059, 4.65696e-057, 4.49802e-055, 4.31695e-053,
	4.07311e-051, 3.81596e-049, 3.61043e-047, 1.73046e-045, 8.26375e-044,
	3.9325e-042, 1.86471e-040, 8.8102e-039, 4.14402e-037, 1.99497e-035,
	1.83001e-034, 1.59023e-033, 1.45505e-032, 1.24603e-031, 1.15674e-030,
	9.70832e-030, 9.23876e-029, 4.20763e-028, 4.24611e-027, 1.61371e-026,
	6.59556e-026, 3.17069e-025, 1.12205e-024, 4.65874e-024, 2.01267e-023,
	6.2941e-023, 3.02604e-022, 7.84622e-022, 2.3526e-021, 6.7049e-021,
	1.88634e-020, 4.59378e-020, 1.37233e-019, 4.00682e-019, 8.34209e-019,
	2.21612e-018, 4.84252e-018, 1.02457e-017, 2.03289e-017, 4.07704e-017,
	1.33778e-016, 2.4263e-016, 4.14981e-016, 7.0383e-016, 1.20511e-015,
	3.85644e-015, 6.52861e-015, 1.06563e-014, 1.67897e-014, 2.79916e-014,
	4.54319e-014, 9.83296e-014, 1.66278e-013, 2.61858e-013, 4.03872e-013,
	5.98967e-013, 1.09674e-012, 1.70553e-012, 2.56573e-012, 3.72723e-012,
	6.14029e-012, 9.33636e-012, 1.36469e-011, 1.89881e-011, 2.68391e-011,
	4.12016e-011, 5.94394e-011, 8.43746e-011, 1.12903e-010, 1.66987e-010,
	2.36959e-010, 3.11726e-010, 4.28713e-010, 5.90781e-010, 7.79892e-010,
	1.05264e-009, 1.4016e-009, 1.87506e-009, 2.42521e-009, 3.14508e-009,
	4.38605e-009, 5.43307e-009, 6.96737e-009, 8.84136e-009, 1.16286e-008,
	1.42343e-008, 1.79697e-008, 2.30867e-008, 2.88832e-008, 3.52583e-008,
	4.31032e-008, 5.46444e-008, 6.66625e-008, 8.06132e-008, 1.00085e-007,
	1.20952e-007, 1.4816e-007, 1.80608e-007, 2.13125e-007, 2.5324e-007,
	3.094e-007, 3.64545e-007, 4.31692e-007, 5.19078e-007, 6.03409e-007,
	7.21811e-007, 8.53856e-007, 9.71749e-007, 1.13949e-006, 1.37042e-006,
	1.53831e-006, 1.79066e-006, 2.15143e-006, 2.40216e-006, 2.76872e-006,
	3.20825e-006, 3.61263e-006, 4.21315e-006, 4.76404e-006, 5.43261e-006,
	6.2041e-006, 6.96243e-006, 7.94979e-006, 8.89079e-006, 1.01387e-005,
	1.13376e-005, 1.2901e-005, 1.44183e-005, 1.59912e-005, 1.79752e-005,
	1.99171e-005, 2.22665e-005, 2.47802e-005, 2.7678e-005, 3.0492e-005,
	3.34189e-005, 3.71902e-005, 4.12605e-005, 4.54706e-005, 4.98411e-005,
	5.48979e-005, 6.06015e-005, 6.61278e-005, 7.22258e-005, 7.97193e-005,
	8.66574e-005, 9.48075e-005, 0.00010321, 0.000112479, 0.000121776,
	0.000133344, 0.000144023, 0.000156667, 0.000168318, 0.000183192,
	0.000196527, 0.00021395, 0.000228389, 0.000249223, 0.000264372,
	0.000289384, 0.000305707, 0.000333992, 0.000353287, 0.000379868,
	0.000408274, 0.00043638, 0.000465319, 0.000496504, 0.000530376,
	0.000566008, 0.000602621, 0.000642286, 0.000684543, 0.000723853,
	0.000772655, 0.000819418, 0.000868533, 0.000920399, 0.000975529,
	0.00103188, 0.00109478, 0.00115777, 0.00122087, 0.00128857,
	0.00136288, 0.00143557, 0.00151714, 0.00159747, 0.00167572,
	0.00176556, 0.00186199, 0.00195063, 0.00205239, 0.00216102,
	0.00225698, 0.00236962, 0.00249145, 0.00259636, 0.00272455,
	0.00287006, 0.00297545, 0.00312346, 0.0032634, 0.00340298,
	0.00355827, 0.00371195, 0.00387288, 0.00404725, 0.00420016,
	0.00439746, 0.00456332, 0.00475936, 0.00495702, 0.00514683,
	0.00535284, 0.00557904, 0.00578084, 0.00601028, 0.00623082,
	0.00647765, 0.00673499, 0.00696553, 0.00722529, 0.00748878,
	0.00775537, 0.00803271, 0.00832199, 0.00861612, 0.00889863,
	0.00919876, 0.00953343, 0.00985465, 0.0101993, 0.0105042, 0.0108325,
	0.0112019, 0.0115901, 0.0119295, 0.0123009, 0.0127191, 0.0130652,
	0.0134855, 0.0139187, 0.0142929, 0.0147541, 0.0151354, 0.0156087,
	0.0160572, 0.0165382, 0.0169669, 0.0174693, 0.017946, 0.0184202,
	0.0189555, 0.0194336, 0.0200107, 0.0204863, 0.0210242, 0.0216053,
	0.0221361, 0.0226858, 0.0232693, 0.0239027, 0.0244779, 0.025081,
	0.0257169, 0.0263059, 0.0269213, 0.0275533, 0.0282065, 0.0289028,
	0.029567, 0.0302268, 0.0309193, 0.0316619, 0.0323147, 0.0330398,
	0.0338124, 0.0345267, 0.0353038, 0.0360947, 0.0368288, 0.0376202,
	0.0383784, 0.0391894, 0.0399684, 0.0408148, 0.0416403, 0.042545,
	0.0433662, 0.0442498, 0.0451003, 0.046035, 0.0468801, 0.0478059,
	0.0487442, 0.0496647, 0.0505752, 0.0515123, 0.0524792, 0.0534474,
	0.0544682, 0.0554579, 0.0565024, 0.0574619, 0.0584757, 0.0595123,
	0.0605988, 0.0615874, 0.062719, 0.0637876, 0.064883, 0.0659551,
	0.0670567, 0.0681256, 0.0692764, 0.0704584, 0.0715399, 0.0727237,
	0.0738803, 0.0750377, 0.0762275, 0.0773855, 0.0785934, 0.0797802,
	0.0810061, 0.0822205, 0.0834827, 0.084714, 0.0858734, 0.0871999,
	0.0884137, 0.0896948, 0.090982, 0.0922797, 0.093635, 0.0948243,
	0.0961283, 0.0974718, 0.0988291, 0.100097, 0.101433, 0.102847,
	0.104222, 0.105492, 0.106885, 0.10833, 0.109672, 0.111048, 0.112438,
	0.113857, 0.115311, 0.11673, 0.118133, 0.119519, 0.12099, 0.122452,
	0.123905, 0.125445, 0.126852, 0.128326, 0.129793, 0.131277, 0.132817,
	0.134305, 0.135772, 0.137284, 0.138882, 0.140372, 0.14192, 0.143445,
	0.14494, 0.146515, 0.148145, 0.149653, 0.151199, 0.152879, 0.154368,
	0.155958, 0.157674, 0.159211, 0.160787, 0.16241, 0.164043, 0.165693,
	0.167281, 0.168956, 0.170589, 0.172252, 0.173884, 0.175575, 0.177208,
	0.178873, 0.180599, 0.18224, 0.183975, 0.185654, 0.187363, 0.189106,
	0.190729, 0.19252, 0.194158, 0.195879, 0.197697, 0.199391, 0.201164,
	0.202879, 0.204602, 0.206413, 0.20818, 0.209911, 0.211753, 0.213484,
	0.215263, 0.21705, 0.218869, 0.220677, 0.222384, 0.224253, 0.226071,
	0.227886, 0.229726, 0.231529, 0.233373, 0.235234, 0.237081, 0.238853,
	0.240735, 0.242606, 0.244465, 0.246371, 0.248218, 0.250135, 0.251944,
	0.253836, 0.255708, 0.257578, 0.259568, 0.261424, 0.263308, 0.265313,
	0.26716, 0.269073, 0.271046, 0.272921, 0.274841, 0.276819, 0.278735,
	0.280616, 0.282653, 0.284613, 0.286558, 0.288478, 0.290472, 0.292474,
	0.294459, 0.296379, 0.298382, 0.300357, 0.302378, 0.30434, 0.306853
};

/* This evaluates Dickman's function for any value.  See Knuth vol. 2 */
/* for a description of this function and its use. */

double F (double x)
{
	int	i;

	if (x >= 1.0) return (1.0);
	if (x >= 0.5) return (1.0 + log (x));
	i = (int) (x * 1000.0);
	return (savedF[i] + (x * 1000.0 - i) * (savedF[i+1] - savedF[i]));
}

/* Analyze how well P-1 factoring will perform */

void guess_pminus1_bounds (
	unsigned long p,
	unsigned int how_far_factored,
	int	double_check,		/* True if double-checking */
	unsigned long *bound1,
	unsigned long *bound2,
	unsigned long *squarings,
	double	*success_rate)
{
	unsigned long B1, B2, fftlen, vals;
	unsigned int h;
	double	pass1_squarings, pass2_squarings;
	double	logB1, logB2, k, logk, temp, logtemp, log2;
	double	size, prob, gcd_cost, ll_tests, numprimes;
	struct {
		unsigned long B1;
		unsigned long B2;
		double	prob;
		double	pass1_squarings;
		double	pass2_squarings;
	} best[2];

/* Balance P-1 against 1 or 2 LL tests (actually more since we get a */
/* corrupt result reported some of the time). */

	ll_tests = 2.0 + 2 * ERROR_RATE;
	if (double_check) ll_tests = 1.0 + 2 * ERROR_RATE;

/* Precompute the cost of a GCD.  We used Excel to come up with the */
/* formula GCD is equivalent to 861 * Ln (p) - 7775 transforms. */
/* Since one squaring equals two transforms we get the formula below. */
/* NOTE: In version 22, the GCD speed has approximately doubled.  I've */
/* adjusted the formula accordingly. */

	gcd_cost = (430.5 * log (p) - 3887.5) / 2.0;
	if (gcd_cost < 50.0) gcd_cost = 50.0;

/* Compute how many temporaries we can use given our memory constraints. */

	fftlen = map_exponent_to_fftlen (p, 0);
	temp = (double) max_mem () * 1000000.0 -
		(double) map_fftlen_to_memused (fftlen, 0);
	size = (double) gwnum_size (fftlen);
	if (temp <= size) vals = 1;
	else vals = (unsigned long) (temp / size);

/* Find the best B1 */

	log2 = log (2);
	for (B1 = 10000; ; B1 += 5000) {

/* Constants */

	logB1 = log (B1);

/* Compute how many squarings will be required in pass 1 */

	pass1_squarings = ceil (1.44 * B1);

/* Try a lot of B2 values */

	for (B2 = B1; ; B2 += B1 >> 2) {

/* Compute how many squarings will be required in pass 2.  In the */
/* low-memory cases, assume choose_pminus1_plan will pick D = 210, E = 1 */
/* If more memory is available assume choose_pminus1_plan will pick */
/* D = 2310, E = 2.  This will provide an accurate enough cost for our */
/* purposes even if different D and E values are picked.  See */
/* choose_pminus1_plan for a description of the costs of P-1 stage 2. */

	logB2 = log (B2);
	numprimes = (unsigned long) (B2 / (logB2 - 1.0) - B1 / (logB1 - 1.0));
	if (B2 <= B1) {
		pass2_squarings = 0.0;
	} else if (vals <= 8) {		/* D = 210, E = 1, passes = 48/temps */
		unsigned long num_passes;
		num_passes = (unsigned long) ceil (48.0 / (vals - 3));
		pass2_squarings = ceil ((B2 - B1) / 210.0) * num_passes;
		pass2_squarings += numprimes * 1.1;
	} else {
		unsigned long num_passes;
		double	numpairings;
		num_passes = (unsigned long) ceil (480.0 / (vals - 5));
		numpairings = (unsigned long)
			(numprimes / 2.0 * numprimes / ((B2-B1) * 480.0/2310.0));
		pass2_squarings = 2400.0 + num_passes * 90.0; /* setup costs */
		pass2_squarings += ceil ((B2-B1) / 4620.0) * 2.0 * num_passes;
		pass2_squarings += numprimes - numpairings;
	}

/* Pass 2 FFT multiplications seem to be at least 20% slower than */
/* the squarings in pass 1.  This is probably due to several factors. */
/* These include: better L2 cache usage and no calls to the faster */
/* gwsquare routine. */

	pass2_squarings *= 1.2;

/* What is the average value of k in 2kp+1 s.t. you get a number */
/* between 2^how_far_factored and 2^(how_far_factored+1)? */

	k = 1.5 * pow (2.0, how_far_factored) / 2.0 / p;
	logk = log (k);

/* Set temp to the number that will need B1 smooth if k has an */
/* average-sized factor found in stage 2 */

	temp = k / ((B1 + B2) / 2);
	logtemp = log (temp);

/* Loop over increasing bit lengths for the factor */

	prob = 0.0;
	for (h = how_far_factored; ; ) {
		double	prob1, prob2;

/* See how many smooth k's we should find using B1 */
/* Using Dickman's function (see Knuth pg 382-383) we want k^a <= B1 */

		prob1 = F (logB1 / logk);

/* See how many smooth k's we should find using B2 */
/* Adjust this slightly to eliminate k's that have two primes > B1 and < B2 */
/* Do this by assuming the largest factor is the average of B1 and B2 */
/* and the remaining cofactor is B1 smooth */

		prob2 = prob1 + (F (logB2 / logk) - prob1) *
				(F (logB1 / logtemp) / F (logB2 / logtemp));
		if (prob2 < 0.001) break;

/* Add this data in to the total chance of finding a factor */

		h++;
		logk += log2;
		logtemp += log2;
		prob += prob2 / h;
	}

/* See if this is a new best case scenario */

	if (B2 == B1 ||
	    prob * ll_tests * p - pass2_squarings >
			best[0].prob * ll_tests * p - best[0].pass2_squarings){
		best[0].B2 = B2;
		best[0].prob = prob;
		best[0].pass2_squarings = pass2_squarings;
		if (vals < 4) break;
		continue;
	}

	if (prob * ll_tests * p - pass2_squarings <
		0.9 * (best[0].prob * ll_tests * p - best[0].pass2_squarings))
		break;
	continue;
	}

/* Is this the best B1 thusfar? */

	if (B1 == 10000 ||
	    best[0].prob * ll_tests * p -
			(pass1_squarings + best[0].pass2_squarings) >
		best[1].prob * ll_tests * p -
			(best[1].pass1_squarings + best[1].pass2_squarings)) {
		best[1].B1 = B1;
		best[1].B2 = best[0].B2;
		best[1].prob = best[0].prob;
		best[1].pass1_squarings = pass1_squarings;
		best[1].pass2_squarings = best[0].pass2_squarings;
		continue;
	}
	if (best[0].prob * ll_tests * p -
			(pass1_squarings + best[0].pass2_squarings) <
	    0.9 * (best[1].prob * ll_tests * p -
			(best[1].pass1_squarings + best[1].pass2_squarings)))
		break;
	continue;
	}

/* Return the final best choice */

	if (best[1].prob * ll_tests * p >
		best[1].pass1_squarings + best[1].pass2_squarings + gcd_cost) {
		*bound1 = best[1].B1;
		*bound2 = best[1].B2;
		*squarings = (unsigned long)
			(best[1].pass1_squarings +
			 best[1].pass2_squarings + gcd_cost);
		*success_rate = best[1].prob;
	} else {
		*bound1 = 0;
		*bound2 = 0;
		*squarings = 0;
		*success_rate = 0.0;
	}
}
