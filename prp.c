/*----------------------------------------------------------------------
| This file contains routines and global variables that are common for
| all operating systems the program has been ported to.  It is included
| in one of the source code files of each port.  See common.h for the
| common #defines and common routine definitions.
+---------------------------------------------------------------------*/

char JUNK[]="Copyright 1996-2002 Just For Fun Software, All rights reserved";
char	INI_FILE[80] = {0};
char	RESFILE[80] = {0};
char	LOGFILE[80] = {0};
char	EXTENSION[8] = {0};

int volatile ERRCHK = 0;
unsigned int PRIORITY = 1;
unsigned int CPU_AFFINITY = 99;
unsigned long volatile ITER_OUTPUT = 0;
unsigned long volatile ITER_OUTPUT_RES = 99999999;
unsigned long volatile DISK_WRITE_TIME = 30;
int	TWO_BACKUP_FILES = 1;
int	RUN_ON_BATTERY = 1;
int	TRAY_ICON = TRUE;
int	HIDE_ICON = FALSE;
unsigned int PRECISION = 2;
int	CUMULATIVE_TIMING = 0;
int	HIGH_RES_TIMER = 0;

int isPRP (char	*, unsigned long, int *);

/* PRP global variables */

giant	N = NULL;		/* Number being factored */


/* Utility output routines */

void LineFeed ()
{
	OutputStr ("\n");
}

void OutputNum (
	unsigned long num)
{
	char	buf[20];
	sprintf (buf, "%lu", num);
	OutputStr (buf);
}

/* Sleep five minutes before restarting */

char ERROK[] = "Disregard last error.  Result is reproducible and thus not a hardware problem.\n";
char ERRMSG0[] = "Bit: %ld/%ld, %s";
char ERRMSG1A[] = "ERROR: ILLEGAL SUMOUT\n";
char ERRMSG1B[] = "ERROR: SUM(INPUTS) != SUM(OUTPUTS), %.16g != %.16g\n";
char ERRMSG1C[] = "ERROR: ROUND OFF (%.10g) > 0.40\n";
char ERRMSG2[] = "Possible hardware failure, consult the readme file.\n";
char ERRMSG3[] = "Continuing from last save file.\n";
char ERRMSG4[] = "Waiting five minutes before restarting.\n";
char WRITEFILEERR[] = "Error writing intermediate file: %s\n";

int SleepFive ()
{
	int	i;

	OutputStr (ERRMSG4);
	BlinkIcon (10);			/* Blink icon for 10 seconds */
	Sleep (10000);
	ChangeIcon (IDLE_ICON);		/* Idle icon for rest of 5 minutes */
	for (i = 0; i < 290; i++) {
		Sleep (1000);
		if (escapeCheck ()) return (FALSE);
	}
	ChangeIcon (WORKING_ICON);	/* And back to the working icon */
	return (TRUE);
}

/* Truncate a percentage to the requested number of digits. */
/* Truncating prevents 99.5% from showing up as 100% complete. */

double trunc_percent (
	double	percent)
{
	if (percent > 100.0) percent = 100.0;
	percent -= 0.5 * pow (10.0, - (double) PRECISION);
	if (percent < 0.0) return (0.0);
	return (percent);
}

/* Routines used to time code chunks */

double timers[10];			/* Up to ten separate timers */

void clear_timers () {
	int	i;
	for (i = 0; i < 10; i++) timers[i] = 0.0;
}

void clear_timer (
	int	i)
{
	timers[i] = 0.0;
}

void start_timer (
	int	i)
{
	if (HIGH_RES_TIMER) {
		timers[i] -= getHighResTimer ();
	} else {
		struct _timeb timeval;
		_ftime (&timeval);
		timers[i] -= (double) timeval.time * 1000.0 + timeval.millitm;
	}
}

void end_timer (
	int	i)
{
	if (HIGH_RES_TIMER) {
		timers[i] += getHighResTimer ();
	} else {
		struct _timeb timeval;
		_ftime (&timeval);
		timers[i] += (double) timeval.time * 1000.0 + timeval.millitm;
	}
}

void divide_timer (
	int	i,
	int	j)
{
	timers[i] = timers[i] / j;
}

double timer_value (
	int	i)
{
	if (HIGH_RES_TIMER)
		return (timers[i] / getHighResTimerFrequency ());
	else
		return (timers[i] / 1000.0);
}

#define TIMER_NL	0x1
#define TIMER_CLR	0x2
#define TIMER_OPT_CLR	0x4

void print_timer (
	int	i,
	int	flags)
{
	char	buf[40];
	double	t;

	t = timer_value (i);
	if (t >= 1.0)
		sprintf (buf, "%.3f sec.", t);
	else
		sprintf (buf, "%.3f ms.", t * 1000.0);
	OutputStr (buf);
	if (flags & TIMER_NL) LineFeed ();
	if (flags & TIMER_CLR) timers[i] = 0.0;
	if ((flags & TIMER_OPT_CLR) && !CUMULATIVE_TIMING) timers[i] = 0.0;
}

/* Determine the CPU speed either empirically or by user overrides. */
/* getCpuType must be called prior to calling this routine. */

void getCpuSpeed (void)
{
	int	temp;

/* Guess the CPU speed using the RDTSC instruction */

	guessCpuSpeed ();

/* Now let the user override the cpu speed from the local.ini file */

	if (IniGetInt (INI_FILE, "CpuOverride", 0)) {
		temp = IniGetInt (INI_FILE, "CpuSpeed", 99);
		if (temp != 99) CPU_SPEED = temp;
	}

/* Make sure the cpu speed is reasonable */

	if (CPU_SPEED > 50000) CPU_SPEED = 50000;
	if (CPU_SPEED < 25) CPU_SPEED = 25;
}

/* Set the CPU flags based on the CPUID data.  Also, the */
/* advanced user can override our guesses. */

void getCpuInfo (void)
{
	int	temp;

/* Get the CPU info using CPUID instruction */

	guessCpuType ();

/* Let the user override the cpu flags from the local.ini file */

	temp = IniGetInt (INI_FILE, "CpuSupportsRDTSC", 99);
	if (temp == 0) CPU_FLAGS &= ~CPU_RDTSC;
	if (temp == 1) CPU_FLAGS |= CPU_RDTSC;
	temp = IniGetInt (INI_FILE, "CpuSupportsCMOV", 99);
	if (temp == 0) CPU_FLAGS &= ~CPU_CMOV;
	if (temp == 1) CPU_FLAGS |= CPU_CMOV;
	temp = IniGetInt (INI_FILE, "CpuSupportsPrefetch", 99);
	if (temp == 0) CPU_FLAGS &= ~CPU_PREFETCH;
	if (temp == 1) CPU_FLAGS |= CPU_PREFETCH;
	temp = IniGetInt (INI_FILE, "CpuSupportsSSE", 99);
	if (temp == 0) CPU_FLAGS &= ~CPU_SSE;
	if (temp == 1) CPU_FLAGS |= CPU_SSE;
	temp = IniGetInt (INI_FILE, "CpuSupportsSSE2", 99);
	if (temp == 0) CPU_FLAGS &= ~CPU_SSE2;
	if (temp == 1) CPU_FLAGS |= CPU_SSE2;

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

/* Determine the names of the INI files */

void nameIniFiles (
	int	named_ini_files)
{
	char	buf[120];

	if (named_ini_files < 0) {
		strcpy (INI_FILE, "prp.ini");
		strcpy (RESFILE, "results.txt");
		strcpy (LOGFILE, "prime.log");
		strcpy (EXTENSION, "");
	} else {
		sprintf (INI_FILE, "prp%04d.ini", named_ini_files);
		sprintf (RESFILE, "resu%04d.txt", named_ini_files);
		sprintf (LOGFILE, "prim%04d.log", named_ini_files);
		sprintf (EXTENSION, ".%03d", named_ini_files);
	}

/* Let the user rename these files and pick a different working directory */

	IniGetString (INI_FILE, "WorkingDir", buf, sizeof(buf), NULL);
	IniGetString (INI_FILE, "results.txt", RESFILE, 80, RESFILE);
	IniGetString (INI_FILE, "prime.log", LOGFILE, 80, LOGFILE);
	IniGetString (INI_FILE, "prime.ini", INI_FILE, 80, INI_FILE);
	if (buf[0]) {
		_chdir (buf);
		IniFileOpen (INI_FILE, 0);
	}
}

/* Read the INI files */

void readIniFiles ()
{
	int	temp;

	getCpuInfo ();

	PRECISION = (unsigned int) IniGetInt (INI_FILE, "PercentPrecision", 2);
	if (PRECISION > 6) PRECISION = 6;

	ITER_OUTPUT = IniGetInt (INI_FILE, "OutputIterations", 10000);
	if (ITER_OUTPUT <= 0) ITER_OUTPUT = 1;
	ITER_OUTPUT_RES = IniGetInt (INI_FILE, "ResultsFileIterations",
				     99999999);
	if (ITER_OUTPUT_RES < 1000) ITER_OUTPUT_RES = 1000;
	DISK_WRITE_TIME = IniGetInt (INI_FILE, "DiskWriteTime", 30);
	TWO_BACKUP_FILES = (int) IniGetInt (INI_FILE, "TwoBackupFiles", 1);
	RUN_ON_BATTERY = (int) IniGetInt (INI_FILE, "RunOnBattery", 1);

	temp = (int) IniGetInt (INI_FILE, "ErrorCheck", 0);
	ERRCHK = (temp != 0);
	PRIORITY = (unsigned int) IniGetInt (INI_FILE, "Priority", 1);
	CPU_AFFINITY = (unsigned int) IniGetInt (INI_FILE, "Affinity", 99);
	HIDE_ICON = (int) IniGetInt (INI_FILE, "HideIcon", 0);
	TRAY_ICON = (int) IniGetInt (INI_FILE, "TrayIcon", 1);

/* Guess the CPU type if it isn't known.  Otherwise, validate it. */

	getCpuInfo ();

/* Other oddball options */

	CUMULATIVE_TIMING = IniGetInt (INI_FILE, "CumulativeTiming", 0);
	HIGH_RES_TIMER = isHighResTimerAvailable ();
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

	if (filelen < 250000) {

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


/* Generate temporary file name */

void tempFileName (
	char	*buf)
{
	giant	tmp, tmp2;

	tmp = popg ((PARG >> 5) + 1);
	tmp2 = popg (2);
	gtog (N, tmp);
	itog (19999981, tmp2);
	modg (tmp2, tmp);
	sprintf (buf, "z%07li", tmp->n[0] % 10000000);
	pushg (2);
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
		LogMsg ("Error opening the results file.\n");
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

/* Output the message */

	if (_write (fd, msg, strlen (msg)) < 0) goto fail;
	_close (fd);
	return (TRUE);

/* On a write error, close file and return error flag */

fail:	_close (fd);
	return (FALSE);
}





/* Convert a value from FFT-ready format to giant */

void gwtobinary (
	gwnum	gg,
	giant	v)
{
	long	val;
	int	j, bits, bitsout, carry;
	unsigned long i, *outptr;

/* Collect bits until we have all of them */

	carry = 0;
	bitsout = 0;
	outptr = v->n;
	*outptr = 0;
	for (i = 0; i < FFTLEN; i++) {
		get_fft_value (gg, i, &val);
		bits = BITS_PER_WORD;
		if (is_big_word (i)) bits++;
		val += carry;
		for (j = 0; j < bits; j++) {
			*outptr >>= 1;
			if (val & 1) *outptr += 0x80000000;
			val >>= 1;
			bitsout++;
			if (bitsout == 32) {
				outptr++;
				*outptr = 0;
				bitsout = 0;
			}
		}
		carry = val;
	}
	*outptr >>= (32 - bitsout);
	outptr++;
	v->sign = (outptr - v->n);

/* Wrap carry from most significant word to least significant word */

	if (PLUS1) carry = -carry;
	outptr = v->n;
	while (carry) {
		RES = *outptr;
		CARRYL = 0;
		if (carry > 0) addhlp (carry);
		else subhlp (-carry);
		*outptr++ = RES;
		carry = CARRYL;
	}

/* Set the length */

	while (v->sign && v->n[v->sign-1] == 0) v->sign--;
}

/* Convert a giant back to the FFT format */

void binarytogw (
	giant	a,
	gwnum	g)
{
	unsigned long i, mask1, mask2, e1len;
	int	neg, bits1, bits2, bits_in_next_binval;
	unsigned long *e1, binval, carry;

	e1len = abs (a->sign);
	e1 = a->n;
	neg = (a->sign < 0);

	if (e1len == 0) {
		dbltogw (0.0, g);
		return;
	}

	bits1 = BITS_PER_WORD;
	bits2 = bits1 + 1;
	mask1 = (1L << bits1) - 1;
	mask2 = (1L << bits2) - 1;
	binval = *e1++; e1len--; bits_in_next_binval = 32;
	carry = (PLUS1 && neg) ? 2 : 0;
	for (i = 0; i < FFTLEN; i++) {
		int	big_word, bits;
		long	value, mask;
		big_word = is_big_word (i);
		bits = big_word ? bits2 : bits1;
		mask = big_word ? mask2 : mask1;
		if (i == FFTLEN - 1) value = binval;
		else value = binval & mask;
		if (neg) value = value ^ mask;
		value = value + carry;
		if (value > (mask >> 1) && bits > 1 && i != FFTLEN - 1) {
			value = value - (mask + 1);
			carry = 1;
		} else {
			carry = 0;
		}
		set_fft_value (g, i, value);

		binval >>= bits;
		if (e1len == 0) continue;
		if (bits_in_next_binval < bits) {
			if (bits_in_next_binval)
				binval |= (*e1 >> (32 - bits_in_next_binval)) << (32 - bits);
			bits -= bits_in_next_binval;
			e1++; e1len--; bits_in_next_binval = 32;
			if (e1len == 0) continue;
		}
		if (bits) {
			binval |= (*e1 >> (32 - bits_in_next_binval)) << (32 - bits);
			bits_in_next_binval -= bits;
		}
	}
	((long *) g)[-1] = 0;	/* Clear needs-normalize counter */
}

/* Special gwtobinary code.  The standard code does not handle negative */
/* gwnums.  This is because it returns a result modulo a Mersenne number - */
/* thanks to the historical roots of gwnums for use in GIMPS. */
/* This code adds a multiple of N to the gwnum and then does a modg call */

void specialgwtobinary (
	gwnum	gw,
	giant	g)
{
	gwnum	tmp;

	tmp = gwalloc ();
	gtog (N, g);
	gshiftleft (BITS_PER_WORD + 1, g);
	binarytogw (g, tmp);
	gwaddquick (gw, tmp);
	gwtobinary (tmp, g);
	modg (N, g);
	gwfree (tmp);
}

/* Read and write intermediate results to a file */

int read_gwnum (
	int	fd,
	gwnum	g,
	long	*sum)
{
	giant	tmp;
	long	i, len, bytes;

	tmp = popg ((PARG >> 5) + 1);
	if (_read (fd, &len, sizeof (long)) != sizeof (long)) return (FALSE);
	bytes = len * sizeof (long);
	if (_read (fd, tmp->n, bytes) != bytes) return (FALSE);
	tmp->sign = len;
	*sum += len;
	for (i = 0; i < len; i++) *sum += tmp->n[i];
	binarytogw (tmp, g);
	pushg (1);
	return (TRUE);
}

int write_gwnum (
	int	fd,
	gwnum	g,
	long	*sum)
{
	giant	tmp;
	long	i, len, bytes;

	tmp = popg ((PARG >> 5) + 1);
	specialgwtobinary (g, tmp);
	len = tmp->sign;
	if (_write (fd, &len, sizeof (long)) != sizeof (long)) return (FALSE);
	bytes = len * sizeof (long);
	if (_write (fd, tmp->n, bytes) != bytes) return (FALSE);
	*sum += len;
	for (i = 0; i < len; i++) *sum += tmp->n[i];
	pushg (1);
	return (TRUE);
}

int read_long (
	int	fd,
	unsigned long *val,
	long	*sum)
{
	if (_read (fd, val, sizeof (long)) != sizeof (long)) return (FALSE);
	*sum += *val;
	return (TRUE);
}

int write_long (
	int	fd,
	unsigned long val,
	long	*sum)
{
	if (_write (fd, &val, sizeof (long)) != sizeof (long)) return (FALSE);
	*sum += val;
	return (TRUE);
}

int writeToFile (
	char	*filename,
	unsigned long j,
	gwnum	x)
{
	char	newfilename[16];
	int	fd;
	unsigned long magicnum, version;
	long	sum = 0;

/* If we are allowed to create multiple intermediate files, then */
/* write to a file called yNNNNNNN. */

	strcpy (newfilename, filename);
	if (TWO_BACKUP_FILES) newfilename[0] = 'y';

/* Create the intermediate file */

	fd = _open (newfilename, _O_BINARY|_O_WRONLY|_O_TRUNC|_O_CREAT, 0666);
	if (fd < 0) return (FALSE);

/* Write the file header. */

	magicnum = 0x9f2b3cd4;
	if (_write (fd, &magicnum, sizeof (long)) != sizeof (long))
		goto writeerr;
	version = 1;
	if (_write (fd, &version, sizeof (long)) != sizeof (long))
		goto writeerr;

/* Write the file data */

	if (! write_long (fd, j, &sum)) goto writeerr;

/* Write the data values */

	if (! write_gwnum (fd, x, &sum)) goto writeerr;

/* Write the checksum */

	if (_write (fd, &sum, sizeof (long)) != sizeof (long)) goto writeerr;
	_commit (fd);
	_close (fd);

/* Now rename the intermediate files */

	if (TWO_BACKUP_FILES) {
		_unlink (filename);
		rename (newfilename, filename);
	}
	return (TRUE);

/* An error occured.  Close and delete the current file. */

writeerr:
	_close (fd);
	_unlink (newfilename);
	return (FALSE);
}

int readFromFile (
	char	*filename,
	unsigned long *j,
	gwnum	x)
{
	int	fd;
	unsigned long magicnum, version;
	long	sum = 0, i;

/* Open the intermediate file */

	fd = _open (filename, _O_BINARY | _O_RDONLY);
	if (fd < 0) goto error;

/* Read the file header */

	if (_read (fd, &magicnum, sizeof (long)) != sizeof (long))
		goto readerr;
	if (magicnum != 0x9f2b3cd4) goto readerr;

	if (_read (fd, &version, sizeof (long)) != sizeof (long)) goto readerr;
	if (version != 1) goto readerr;

/* Read the file data */

	if (! read_long (fd, j, &sum)) goto readerr;

/* Read the values */

	if (! read_gwnum (fd, x, &sum)) goto readerr;

/* Read and compare the checksum */

	if (_read (fd, &i, sizeof (long)) != sizeof (long)) goto readerr;
	if (i != sum) goto readerr;
	_close (fd);
	return (TRUE);

/* An error occured.  Delete the current intermediate file. */
/* Set stage to -1 to indicate an error. */

readerr:
	OutputStr ("Error reading PRP save file.\n");
	_close (fd);
error:
	_unlink (filename);
	return (FALSE);
}

/* Generate the 64-bit residue of a Lucas-Lehmer test */

/*#define GTEST*/
#ifdef GTEST
void gen64 (char *str, gwnum x)
{
	giant	tmp;
	char	buf[200];
	tmp = popg ((PARG >> 5) + 3);
	specialgwtobinary (x, tmp); /* Compensate for possible negative x */
	sprintf (buf, "%s res64: %08lX%08lX\n", str, tmp->n[1], tmp->n[0]);
	writeResults (buf);
	pushg (1);
}
#endif



/* Test if N is a probable proth prime.  The number N can be of ANY form. */
/* Used to test code for Chris Nash's PFORM program. */

int isProthPRP (
	unsigned long a,
	unsigned long k,
	unsigned long n,
	int	inc,			/* Plus or minus one */
	int	*res)
{
	unsigned long bit, len, iters;
	unsigned long p, fftlen, ptmp, pmax;
	unsigned long bits_per_word;
	gwnum	x;
	giant	tmp;
	char	filename[20], buf[100], str[40], res64[17];
	long	write_time = DISK_WRITE_TIME * 60;
	int	echk, saving, stopping;
	time_t	start_time, current_time;
	double	reallyminerr = 1.0;
	double	reallymaxerr = 0.0;

/* Get the current time */

restart:
	*res = TRUE;		/* Assume it is a probable prime */
	time (&start_time);

/* Compute the number we are testing */

	sprintf (str, "%lu*2^%lu%c1", k, n, inc < 0 ? '-' : '+');
	N = newgiant ((n + 64) >> 4);
	ultog (k, N);
	gshiftleft (n, N);
	iaddg (inc, N);
	len = bitlen (N);

/* Assume intermediate results of twice the length of N, plus a few */
/* spare bits.  Then round up so that we are using a rational (integral */
/* number of bits per FFT value).  Finally, setup the assembly code. */

	clear_timers ();
	start_timer (0);
	start_timer (1);
	p = len + len + 64;
	for (ptmp = p / 2; ; ptmp = pmax + 1) {
		fftlen = map_exponent_to_fftlen (ptmp, GW_MERSENNE_MOD);
		pmax = map_fftlen_to_max_exponent (fftlen, GW_MERSENNE_MOD);
		if (p <= (pmax + fftlen / 2) / fftlen * fftlen) break;
	}
	bits_per_word = (p + fftlen - 1) / fftlen;
	p = bits_per_word * fftlen;
	gwsetup (p, fftlen, GW_MERSENNE_MOD);
	x = gwalloc ();

/* Setup the proth mod code */

	gwprothsetup (k, n, inc);

/* Init filename */

	tempFileName (filename);

/* Optionally resume from save file and output a message */
/* indicating we are resuming a test */

	if (fileExists (filename) && readFromFile (filename, &bit, x)) {
		char	fmt_mask[80];
		double	pct;
		pct = trunc_percent (bit * 100.0 / len);
		sprintf (fmt_mask,
			 "Resuming probable prime test of %%s at bit %%ld [%%.%df%%%%]\n",
			 PRECISION);
		sprintf (buf, fmt_mask, str, bit, pct);
		OutputStr (buf);
	}

/* Otherwise, output a message indicating we are starting test */

	else {
		sprintf (buf, "Starting probable prime test of %s\n", str);
		OutputStr (buf);
		bit = 1;
		dbltogw ((double) a, x);
	}

	ReplaceableLine (1);	/* Remember where replaceable line is */

/* Init the title */

	title (str);

/* Do the PRP test */

	gwsetmulbyconst (a);
	iters = 0;
	while (bit < len) {

/* Error check the last 50 iterations, before writing an */
/* intermediate file (either user-requested stop or a */
/* 30 minute interval expired), and every 128th iteration. */

		stopping = stopCheck ();
		echk = stopping || ERRCHK || (bit >= len - 50);
		if ((bit & 127) == 0) {
			echk = 1;
			time (&current_time);
			saving = (current_time - start_time > write_time);
		} else
			saving = 0;

/* Process this bit */

		if (bitval (N, len-bit-1)) {
			gwsetnormroutine (0, echk, 1);
			gwsquare (x);
		} else {
//{giant t1, t2;
//t1 = popg ((PARG >> 4) + 3);
//t2 = popg ((PARG >> 4) + 3);
//specialgwtobinary (x, t1);
//binarytogw (t1, x);
			gwsetnormroutine (0, echk, 0);
			gwsquare (x);
//squareg (t1);
//gwtobinary (x, t2);
//if (gcompg (t1, t2) != 0)
//bit++;
//pushg(2);}
		}
//{giant t1, t2;
//t1 = popg ((PARG >> 4) + 3);
//t2 = popg ((PARG >> 4) + 3);
//gwtobinary (x, t1);
//			end_timer (0);
		gwprothmod (x);
//			start_timer (0);
//specialgwtobinary (x, t2);
//modg (N, t1);
//modg (N, t2);
//if (gcompg (t1, t2) != 0)
//bit++;
//pushg(2);}
		bit++;
		iters++;
#ifdef GTEST
		if (bit == 100) {
			gen64 (str, x);
			term_giants ();
			gwdone ();
			return (isPRP (str, 3, res));
		}
#endif

/* If the sum of the output values is an error (such as infinity) */
/* then raise an error. */

		if (gw_test_illegal_sumout ()) {
			sprintf (buf, ERRMSG0, bit, len, ERRMSG1A);
			OutputBoth (buf);
			goto error;
		}

/* Check that the sum of the input numbers squared is approximately */
/* equal to the sum of unfft results.  Since this check may not */
/* be perfect, check for identical results after a restart. */

		if (gw_test_mismatched_sums ()) {
			static unsigned long last_bit = 0;
			static double last_suminp = 0.0;
			static double last_sumout = 0.0;
			double suminp, sumout;
			suminp = gwsuminp (x);
			sumout = gwsumout (x);
			if (bit == last_bit &&
			    suminp == last_suminp &&
			    sumout == last_sumout) {
				writeResults (ERROK);
				saving = 1;
			} else {
				char	msg[80];
				sprintf (msg, ERRMSG1B, suminp, sumout);
				sprintf (buf, ERRMSG0, bit, len, msg);
				OutputBoth (buf);
				last_bit = bit;
				last_suminp = suminp;
				last_sumout = sumout;
				goto error;
			}
		}

/* Check for excessive roundoff error  */

		if (echk && MAXERR > 0.40) {
			static unsigned long last_bit = 0;
			static double last_maxerr = 0.0;
			if (bit == last_bit &&
			    MAXERR == last_maxerr) {
				writeResults (ERROK);
				saving = 1;
			} else {
				char	msg[80];
				sprintf (msg, ERRMSG1C, MAXERR);
				sprintf (buf, ERRMSG0, bit, len, msg);
				OutputBoth (buf);
				last_bit = bit;
				last_maxerr = MAXERR;
				goto error;
			}
		}

		if (ERRCHK) {
			if (MAXERR < reallyminerr && bit > 30)
				reallyminerr = MAXERR;
			if (MAXERR > reallymaxerr)
				reallymaxerr = MAXERR;
		}

/* Print a message every so often */

		if (bit % ITER_OUTPUT == 0) {
			char	fmt_mask[80];
			double	pct;
			pct = trunc_percent (bit * 100.0 / len);
			sprintf (fmt_mask, "%%.%df%%%% of %%s", PRECISION);
			sprintf (buf, fmt_mask, pct, str);
			title (buf);
			ReplaceableLine (2);	/* Replace line */
			sprintf (fmt_mask,
				 "%%s, bit: %%ld / %%ld [%%.%df%%%%]",
				 PRECISION);
			sprintf (buf, fmt_mask, str, bit, len, pct);
			OutputStr (buf);
			if (ERRCHK && bit > 30) {
				OutputStr (".  Round off: ");
				sprintf (buf, "%10.10f", reallyminerr);
				OutputStr (buf);
				sprintf (buf, " to %10.10f", reallymaxerr);
				OutputStr (buf);
			}
			end_timer (0);
			if (CUMULATIVE_TIMING) {
				OutputStr (".  Time thusfar: ");
			} else {
				OutputStr (".  Time per bit: ");
				divide_timer (0, iters);
				iters = 0;
			}
			print_timer (0, TIMER_NL | TIMER_OPT_CLR);
			start_timer (0);
		}

/* Print a results file message every so often */

		if (bit % ITER_OUTPUT_RES == 0 || (NO_GUI && stopping)) {
			sprintf (buf, "Bit %ld / %ld\n", bit, len);
			writeResults (buf);
		}

/* Write results to a file every DISK_WRITE_TIME minutes */
/* On error, retry in 10 minutes (it could be a temporary */
/* disk-full situation) */

		if (saving || stopping) {
			write_time = DISK_WRITE_TIME * 60;
			if (! writeToFile (filename, bit, x)) {
				sprintf (buf, WRITEFILEERR, filename);
				OutputBoth (buf);
				if (write_time > 600) write_time = 600;
			}
			time (&start_time);

/* If an escape key was hit, write out the results and return */

			if (stopping) return (FALSE);
		}
	}
	tmp = popg ((PARG >> 5) + 3);
	specialgwtobinary (x, tmp); /* Compensate for possible negative x */
	ulsubg (a, tmp);
	if (!isZero (tmp)) {
		*res = FALSE;	/* Not a prime */
		sprintf (res64, "%08lX%08lX", tmp->n[1], tmp->n[0]);
	}
	pushg (1);
	gwfree (x);

/* Cleanup */

	end_timer (1);
	if (*res)
		sprintf (buf, "%s is a probable prime.\n", str);
	else
		sprintf (buf, "%s is not prime.  Res64: %s\n", str, res64);
	if ((*res && IniGetInt (INI_FILE, "OutputPrimes", 1)) ||
	    (!*res && IniGetInt (INI_FILE, "OutputComposites", 1)))
		writeResults (buf);
	sprintf (buf+strlen(buf)-1, "  Time: ");
	ReplaceableLine (2);	/* Replace line */
	OutputStr (buf);
	print_timer (1, TIMER_CLR | TIMER_NL);
	free (N);
	term_giants ();
	gwdone ();
	_unlink (filename);
	return (TRUE);

/* An error occured, sleep, then try restarting at last save point. */

error:	free (N);
	term_giants ();
	gwdone ();

/* Output a message saying we are restarting */

	OutputBoth (ERRMSG2);
	OutputBoth (ERRMSG3);

/* Sleep five minutes before restarting */

	if (! SleepFive ()) return (FALSE);

/* Restart */

	goto restart;
}


/* Test if N is a probable prime.  The number N can be of ANY form. */
/* Used to test code for Chris Nash's PFORM program. */

int isPRP (
	char	*str,		/* string representation of N */
	unsigned long a,
	int	*res)
{
#define EB	10		/* Extra bits of precision */
	unsigned long bit, len, iters;
	unsigned long p, fftlen, ptmp, pmax, zerowordslow, zerowordshigh;
	unsigned long bits_per_word;
	gwnum	x, y, recip, n;
	giant	tmp;
	char	filename[20], buf[100], res64[17];
	long	write_time = DISK_WRITE_TIME * 60;
	int	echk, saving, stopping;
	time_t	start_time, current_time;
	double	reallyminerr = 1.0;
	double	reallymaxerr = 0.0;

/* Get the current time */

restart:
	*res = TRUE;		/* Assume it is a probable prime */
	time (&start_time);

/* Assume intermediate results of twice the length of N, plus a few */
/* spare bits.  Then round up so that we are using a rational (integral */
/* number of bits per FFT value).  Finally, setup the assembly code. */

	clear_timers ();
	start_timer (0);
	start_timer (1);
	len = bitlen (N);
	p = len + len + 2*EB + 64;
	for (ptmp = p / 2; ; ptmp = pmax + 1) {
		fftlen = map_exponent_to_fftlen (ptmp, GW_MERSENNE_MOD);
		pmax = map_fftlen_to_max_exponent (fftlen, GW_MERSENNE_MOD);
		if (p <= (pmax + fftlen / 2) / fftlen * fftlen) break;
	}
	bits_per_word = (p + fftlen - 1) / fftlen;
	p = bits_per_word * fftlen;
	gwsetup (p, fftlen, GW_MERSENNE_MOD);

	n = gwalloc ();
	binarytogw (N, n);
	gwfft (n, n);
	x = gwalloc ();

/* Precompute the reciprocal */

	tmp = newgiant ((p >> 4) + 1);
	itog (1, tmp);
	gshiftleft (len + len + EB, tmp);
	divg (N, tmp);			/* computes len+EB+1 bits of reciprocal */
	gshiftleft (p - len - len - EB, tmp);/* shift so gwmul routines wrap */
					/* quotient to lower end of fft */
	recip = gwalloc ();
	binarytogw (tmp, recip);
	gwfft (recip, recip);
	free (tmp);

/* Init filename */

	tempFileName (filename);

/* Optionally resume from save file and output a message */
/* indicating we are resuming a test */

	if (fileExists (filename) && readFromFile (filename, &bit, x)) {
		char	fmt_mask[80];
		double	pct;
		pct = trunc_percent (bit * 100.0 / len);
		sprintf (fmt_mask,
			 "Resuming probable prime test of %%s at bit %%ld [%%.%df%%%%]\n",
			 PRECISION);
		sprintf (buf, fmt_mask, str, bit, pct);
		OutputStr (buf);
	}

/* Otherwise, output a message indicating we are starting test */

	else {
		sprintf (buf, "Starting probable prime test of %s\n", str);
		OutputStr (buf);
		bit = 1;
		dbltogw ((double) a, x);
	}

	ReplaceableLine (1);	/* Remember where replaceable line is */

/* Init the title */

	title (str);

/* Do the PRP test */

	zerowordslow = (len - EB) / bits_per_word;
	zerowordshigh = fftlen - len / bits_per_word - 1;
	gwsetmulbyconst (a);
	y = gwalloc ();
	iters = 0;
	while (bit < len) {

/* Error check the last 50 iterations, before writing an */
/* intermediate file (either user-requested stop or a */
/* 30 minute interval expired), and every 128th iteration. */

		stopping = stopCheck ();
		echk = stopping || ERRCHK || (bit >= len - 50);
		if ((bit & 127) == 0) {
			echk = 1;
			time (&current_time);
			saving = (current_time - start_time > write_time);
		} else
			saving = 0;

/* Process this bit */

		if (bitval (N, len-bit-1)) {
			gwsetnormroutine (0, echk, 1);
			gwsquare (x);
		} else {
			gwsetnormroutine (0, echk, 0);
			gwsquare (x);
		}

		gwcopyzero (x, y, zerowordslow);
		gwsetnormroutine (zerowordshigh, echk, 0);
		gwfftmul (recip, y);

		gwsetnormroutine (0, echk, 0);
		gwfftmul (n, y);
		gwsub (y, x);

		bit++;
		iters++;
#ifdef GTEST
		if (bit == 100) {
			gen64 (str, x);
			term_giants ();
			gwdone ();
			return (TRUE);
		}
#endif

/* If the sum of the output values is an error (such as infinity) */
/* then raise an error. */

		if (gw_test_illegal_sumout ()) {
			sprintf (buf, ERRMSG0, bit, len, ERRMSG1A);
			OutputBoth (buf);
			goto error;
		}

/* Check that the sum of the input numbers squared is approximately */
/* equal to the sum of unfft results.  Since this check may not */
/* be perfect, check for identical results after a restart. */

		if (gw_test_mismatched_sums ()) {
			static unsigned long last_bit = 0;
			static double last_suminp = 0.0;
			static double last_sumout = 0.0;
			double suminp, sumout;
			suminp = gwsuminp (x) + gwsuminp (y);
			sumout = gwsumout (x) + gwsumout (y);
			if (bit == last_bit &&
			    suminp == last_suminp &&
			    sumout == last_sumout) {
				writeResults (ERROK);
				saving = 1;
			} else {
				char	msg[80];
				sprintf (msg, ERRMSG1B, suminp, sumout);
				sprintf (buf, ERRMSG0, bit, len, msg);
				OutputBoth (buf);
				last_bit = bit;
				last_suminp = suminp;
				last_sumout = sumout;
				goto error;
			}
		}

/* Check for excessive roundoff error  */

		if (echk && MAXERR > 0.40) {
			static unsigned long last_bit = 0;
			static double last_maxerr = 0.0;
			if (MAXERR < 0.48 &&
			    bit == last_bit &&
			    MAXERR == last_maxerr) {
				writeResults (ERROK);
				saving = 1;
			} else {
				char	msg[80];
				sprintf (msg, ERRMSG1C, MAXERR);
				sprintf (buf, ERRMSG0, bit, len, msg);
				OutputBoth (buf);
				last_bit = bit;
				last_maxerr = MAXERR;
				goto error;
			}
		}

		if (ERRCHK) {
			if (MAXERR < reallyminerr && bit > 30)
				reallyminerr = MAXERR;
			if (MAXERR > reallymaxerr)
				reallymaxerr = MAXERR;
		}

/* Print a message every so often */

		if (bit % ITER_OUTPUT == 0) {
			char	fmt_mask[80];
			double	pct;
			pct = trunc_percent (bit * 100.0 / len);
			sprintf (fmt_mask, "%%.%df%%%% of %%s", PRECISION);
			sprintf (buf, fmt_mask, pct, str);
			title (buf);
			ReplaceableLine (2);	/* Replace line */
			sprintf (fmt_mask,
				 "%%s, bit: %%ld / %%ld [%%.%df%%%%]",
				 PRECISION);
			sprintf (buf, fmt_mask, str, bit, len, pct);
			OutputStr (buf);
			if (ERRCHK && bit > 30) {
				OutputStr (".  Round off: ");
				sprintf (buf, "%10.10f", reallyminerr);
				OutputStr (buf);
				sprintf (buf, " to %10.10f", reallymaxerr);
				OutputStr (buf);
			}
			end_timer (0);
			if (CUMULATIVE_TIMING) {
				OutputStr (".  Time thusfar: ");
			} else {
				OutputStr (".  Time per bit: ");
				divide_timer (0, iters);
				iters = 0;
			}
			print_timer (0, TIMER_NL | TIMER_OPT_CLR);
			start_timer (0);
		}

/* Print a results file message every so often */

		if (bit % ITER_OUTPUT_RES == 0 || (NO_GUI && stopping)) {
			sprintf (buf, "Bit %ld / %ld\n", bit, len);
			writeResults (buf);
		}

/* Write results to a file every DISK_WRITE_TIME minutes */
/* On error, retry in 10 minutes (it could be a temporary */
/* disk-full situation) */

		if (saving || stopping) {
			write_time = DISK_WRITE_TIME * 60;
			if (! writeToFile (filename, bit, x)) {
				sprintf (buf, WRITEFILEERR, filename);
				OutputBoth (buf);
				if (write_time > 600) write_time = 600;
			}
			time (&start_time);

/* If an escape key was hit, write out the results and return */

			if (stopping) return (FALSE);
		}
	}
	tmp = popg ((PARG >> 5) + 1);
	specialgwtobinary (x, tmp);
	ulsubg (a, tmp);
	if (!isZero (tmp)) {
		*res = FALSE;	/* Not a prime */
		sprintf (res64, "%08lX%08lX", tmp->n[1], tmp->n[0]);
	}
	pushg (1);
	gwfree (x);
	gwfree (y);

/* Cleanup */

	end_timer (1);
	if (*res)
		sprintf (buf, "%s is a probable prime.\n", str);
	else
		sprintf (buf, "%s is not prime.  Res64: %s\n", str, res64);
	if ((*res && IniGetInt (INI_FILE, "OutputPrimes", 1)) ||
	    (!*res && IniGetInt (INI_FILE, "OutputComposites", 1)))
		writeResults (buf);
	sprintf (buf+strlen(buf)-1, "  Time: ");
	ReplaceableLine (2);	/* Replace line */
	OutputStr (buf);
	print_timer (1, TIMER_CLR | TIMER_NL);
	term_giants ();
	gwdone ();
	_unlink (filename);
	return (TRUE);

/* An error occured, sleep, then try restarting at last save point. */

error:	term_giants ();
	gwdone ();

/* Output a message saying we are restarting */

	OutputBoth (ERRMSG2);
	OutputBoth (ERRMSG3);

/* Sleep five minutes before restarting */

	if (! SleepFive ()) return (FALSE);

/* Restart */

	goto restart;
}

/* Process a number from newpgen output file */

int process_num (
	unsigned long k,
	unsigned long base,
	unsigned long n,
	int	incr,
	int	*res)
{
	if (base == 2 && (incr == -1 || incr == +1)) {
		return (isProthPRP (3, k, n, incr, res));
	} else {
		char	buf[100];
		int	bits, retval;
		if (k == 1)
			sprintf (buf, "%lu^%lu%c%lu", base, n,
				 incr < 0 ? '-' : '+', abs(incr));
		else
			sprintf (buf, "%lu*%lu^%lu%c%lu", k, base, n,
				 incr < 0 ? '-' : '+', abs(incr));
		bits = (int) ((n * log(base) + log(k)) / log(2));
		N = newgiant ((bits >> 4) + 8);
		ultog (base, N);
		power (N, n);
		ulmulg (k, N);
		iaddg (incr, N);
		retval = isPRP (buf, 3, res);
		free (N);
		return (retval);
	}
}



void primeContinue ()
{
	int	work;

/* If we are done, return */

	if (IniGetInt (INI_FILE, "WorkDone", 1)) return;

/* Set appropriate priority */

	SetPriority ();

/* Case off the work type */

	work = IniGetInt (INI_FILE, "Work", 0);

/* Handle a newpgen output file */

	if (work == 0) {
		char	inputfile[80], outputfile[80];
		FILE *fd;
		unsigned long i, chainlen, k, n, base, nfudge, mask;
		int	firstline, line, outfd, res;
		char	c;

		IniGetString (INI_FILE, "PgenInputFile", inputfile, 80, NULL);
		IniGetString (INI_FILE, "PgenOutputFile", outputfile, 80, NULL);
		firstline = IniGetInt (INI_FILE, "PgenLine", 1);

		fd = fopen (inputfile, "r");
		if (fd == NULL) {
			IniWriteInt (INI_FILE, "WorkDone", 1);
			return;
		}

		mask = 0;
		fscanf (fd, "%lu:%c:%lu:%lu:%lu\n", &i, &c, &chainlen, &base, &mask);
		if (mask & 0x40) {
			OutputStr ("Primoral NewPgen files are not supported.\n");
			return;
		}
		if (chainlen == 0) chainlen = 1;

		if (! fileExists (outputfile)) {
			outfd = _open (outputfile, _O_TEXT | _O_RDWR | _O_CREAT, 0666);
			if (outfd) {
				char	buf[100];
				if (mask == 0)
					sprintf (buf, "%lu:%c:%lu:%lu\n", i, c, chainlen, base);
				else
					sprintf (buf, "%lu:%c:%lu:%lu:%lu\n", i, c, chainlen, base, mask);
				_write (outfd, buf, strlen (buf));
				_close (outfd);
			}
		}


/* THIS SECTION IS FOR BACKWARDS COMPATIBILITY WITH PREVIOUS PRP.EXE */
/* That version used the one character code to determine what to do. */
/* The new version uses the mask field. */

		if (mask == 0 || IniGetInt (INI_FILE, "UseCharCode", 0)) {

/* The variable c is a one character code as follows: */
/* P : k.b^n+1 (Plus)
   M : k.b^n-1 (Minus)
   T: k.b^n+-1 (Twin)
   S: k.b^n-1; k.b^(n+1)-1 (SG (CC 1st kind len 2))
   C: k.b^n+1; k.b^(n+1)+1 (CC 2nd kind len 2)
   B: k.b^n+-1; k.b^(n+1)+-1 (BiTwin)
   J: k.b^n+-1; k.b^(n+1)-1 (Twin/SG)
   K: k.b^n+-1; k.b^(n+1)+1 (Twin/CC)
   Y : k.b^n+1 + others (Lucky Plus)
   Z : k.b^n-1 + others (Lucky Minus)
   1: CC 1st kind chain
   2: CC 2nd kind chain
   3: BiTwin chain */
/* Undo the increment of n that newpgen did on types 1, 2, 3 */
/* Map P, M, Y, Z, T, S, C, B to their more generic counterparts */

		nfudge = 0;
		if (c == '1') nfudge = 1;
		if (c == '2') nfudge = 1;
		if (c == '3') nfudge = 1;
		if (c == 'P') c = '2', chainlen = 1;
		if (c == 'M') c = '1', chainlen = 1;
		if (c == 'Y') c = '2', chainlen = 1;
		if (c == 'Z') c = '1', chainlen = 1;
		if (c == 'T') c = '3', chainlen = 1;
		if (c == 'S') c = '1', chainlen = 2;
		if (c == 'C') c = '2', chainlen = 2;
		if (c == 'B') c = '3', chainlen = 2;

/* Process each line in the newpgen output file */

		for (line = 1; ; line++) {

/* Read the line, break at EOF */

			k = 0;
			fscanf (fd, "%lu %lu\n", &k, &n);
			if (k == 0) {
				IniWriteInt (INI_FILE, "WorkDone", 1);
				break;
			}

/* Skip this line if requested (we processed it on an earlier run) */

			if (line < firstline) continue;

/* Test numbers according to the c variable */

			for (i = 0; i < chainlen; i++) {
				if (c == '1' || c == '3') {
					if (! process_num (k, base, n - nfudge + i, -1, &res)) goto done;
					if (!res) break;
				}
				if (c == '2' || c == '3') {
					if (! process_num (k, base, n - nfudge + i, +1, &res)) goto done;
					if (!res) break;
				}
				if (c == 'J') {
					int	res2;
					if (! process_num (k, base, n, -1, &res)) goto done;
					if (!res) break;
					if (! process_num (k, base, n, +1, &res)) goto done;
					if (! process_num (k, base, n+1, -1, &res2)) goto done;
					res |= res2;
					break;
				}
				if (c == 'K') {
					int	res2;
					if (! process_num (k, base, n, +1, &res)) goto done;
					if (!res) break;
					if (! process_num (k, base, n, -1, &res)) goto done;
					if (! process_num (k, base, n+1, +1, &res2)) goto done;
					res |= res2;
					break;
				}
				if (c == 'A') {
					if (! process_num (1, base, n, k+k-1, &res)) goto done;
					if (!res) break;
				}
			}

/* If all numbers tested were probable primes, copy the line to the output file */

			if (res) {
				outfd = _open (outputfile, _O_TEXT | _O_RDWR | _O_APPEND | _O_CREAT, 0666);
				if (outfd) {
					char	buf[100];
					sprintf (buf, "%lu %lu\n", k, n);
					_write (outfd, buf, strlen (buf));
					_close (outfd);
				}
			}

			IniWriteInt (INI_FILE, "PgenLine", line + 1);
		}


/* THIS IS THE NEW SECTION.  It uses both the mask field and the */
/* character code to determine what to do */

		} else {
			unsigned long kk, nn;

/* NEWPGEN output files use the mask as defined below: */
#define MODE_PLUS    0x01	/* k.b^n+1 */
#define MODE_MINUS   0x02	/* k.b^n-1 */
#define MODE_2PLUS   0x04	/* k.b^(n+1)+1 (*) */
#define MODE_2MINUS  0x08	/* k.b^(n+1)-1 (*) */
#define MODE_4PLUS   0x10	/* k.b^(n+2)+1 (*) */
#define MODE_4MINUS  0x20	/* k.b^(n+2)-1 (*) */
#define MODE_PRIMORIAL 0x40	/* PRIMORIAL - can't handle this */
#define MODE_PLUS5  0x80	/* k.b^n+5 */
#define MODE_AP	    0x200	/* 2^n+2k-1 */
#define MODE_PLUS7  0x800	/* k.b^n+7 */
#define MODE_2PLUS3 0x1000	/* 2k.b^n+3 */
#define MODE_NOTGENERALISED 0x400
/* Those entries that have a (*) next to them are modified if the */
/* MODE_NOTGENERALISED flag is set.  If it is set, they are changed */
/* as follows: */
/* MODE_2PLUS      2k.b^n+1 */
/* MODE_2MINUS     2k.b^n-1 */
/* MODE_4PLUS      4k.b^n+1 */
/* MODE_4MINUS     4k.b^n-1 */
/* Similarly, longer chains are affected in the same way (so if the base */
/* is 3 and we are after a CC of the 1st kind of length 4, rather that */
/* looking at k.3^n-1 & k.3^(n+1)-1 & k.3^(n+2)-1 & k.3^(n+3)-1 we look */
/* at k.3^n-1 & 2k.3^n-1 & 4k.3^n-1 & 8k.3^n-1). */

/* Process each line in the newpgen output file */

		for (line = 1; ; line++) {

/* Read the line, break at EOF */

			k = 0;
			fscanf (fd, "%lu %lu\n", &k, &n);
			if (k == 0) {
				IniWriteInt (INI_FILE, "WorkDone", 1);
				break;
			}

/* Skip this line if requested (we processed it on an earlier run) */

			if (line < firstline) continue;

/* Undo the increment of n that newpgen did on types 1, 2, 3 */

			nn = n;
			if (c == '1' || c == '2' || c == '3') nn--;

/* In "lucky" modes, only test the central candidate */

			if (c == 'Y') {
				nn--;
				mask = MODE_2PLUS;
			}
			if (c == 'Z') {
				nn--;
				mask = MODE_2MINUS;
			}

/* Test numbers according to the mask variable */
/* The J and K types (Twin/CC and Twin/SG) are special in that they */
/* are output if either a Twin OR a CC/SG is found */

			kk = k;
			for (i = 0; i < chainlen; i++) {
				if (c == 'J') {
					int	res2;
					if (! process_num (kk, base, nn, -1, &res)) goto done;
					if (!res) break;
					if (! process_num (kk, base, nn, +1, &res)) goto done;
					if (! process_num (kk, base, nn+1, -1, &res2)) goto done;
					res |= res2;
					break;
				}
				if (c == 'K') {
					int	res2;
					if (! process_num (kk, base, nn, +1, &res)) goto done;
					if (!res) break;
					if (! process_num (kk, base, nn, -1, &res)) goto done;
					if (! process_num (kk, base, nn+1, +1, &res2)) goto done;
					res |= res2;
					break;
				}
				if (mask & MODE_MINUS) {
					if (! process_num (kk, base, nn, -1, &res)) goto done;
					if (!res) break;
				}
				if (mask & MODE_PLUS) {
					if (! process_num (kk, base, nn, +1, &res)) goto done;
					if (!res) break;
				}
				if (mask & MODE_PLUS5) {
					if (! process_num (kk, base, nn, +5, &res)) goto done;
					if (!res) break;
				}
				if (mask & MODE_PLUS7) {
					if (! process_num (kk, base, nn, +7, &res)) goto done;
					if (!res) break;
				}
				if (mask & MODE_2PLUS3) {
					if (! process_num (kk+kk, base, nn, +3, &res)) goto done;
					if (!res) break;
				}
				if (mask & MODE_AP) {
					if (! process_num (1, base, nn, kk+kk-1, &res)) goto done;
					if (!res) break;
				}

/* Bump k or n for the next itereation or for the MODE_2PLUS and */
/* MODE_2MINUS flags */

				if (mask & MODE_NOTGENERALISED) kk *= 2;
				else nn += 1;

/* If chainlength is more than 1, then we let the for loop do the work */
/* rather than the MODE_2PLUS, etc. flags */

				if (chainlen > 1) continue;

				if (mask & MODE_2MINUS) {
					if (! process_num (kk, base, nn, -1, &res)) goto done;
					if (!res) break;
				}
				if (mask & MODE_2PLUS) {
					if (! process_num (kk, base, nn, +1, &res)) goto done;
					if (!res) break;
				}

/* Bump k or n for the MODE_4PLUS and MODE_4MINUS flags */

				if (mask & MODE_NOTGENERALISED) kk *= 2;
				else nn += 1;

				if (mask & MODE_4MINUS) {
					if (! process_num (kk, base, nn, -1, &res)) goto done;
					if (!res) break;
				}
				if (mask & MODE_4PLUS) {
					if (! process_num (kk, base, nn, +1, &res)) goto done;
					if (!res) break;
				}
			}

/* If all numbers tested were probable primes, copy the line to the output file */

			if (res) {
				outfd = _open (outputfile, _O_TEXT | _O_RDWR | _O_APPEND | _O_CREAT, 0666);
				if (outfd) {
					char	buf[100];
					sprintf (buf, "%lu %lu\n", k, n);
					_write (outfd, buf, strlen (buf));
					_close (outfd);
				}
			}

			IniWriteInt (INI_FILE, "PgenLine", line + 1);
		}

		}

done:		fclose (fd);
	}

/* Handle an expr */

	else {
		OutputStr ("Expression testing not yet implemented.\n");
		IniWriteInt (INI_FILE, "WorkDone", 1);
	}
}
