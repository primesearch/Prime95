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

/* Global variables */

EXTERNC unsigned long FACPASS = 0;
EXTERNC unsigned long FACHSW = 0;	/* High word of found factor */
EXTERNC unsigned long FACMSW = 0;	/* Middle word of found factor */
EXTERNC unsigned long FACLSW = 0;	/* Low word of found factor */
gwnum LLDATA = NULL;		/* For convienence, the Lucas-Lehmer and */
void *FACDATA = NULL;		/* factor data is kept in a global */
int STOP_REASON = 0;		/* Reason stopCheck stopped processing */
time_t STOP_TIME = 0;		/* Time to stop processing because */
				/* more or less memory is available */
int HIGH_MEMORY_USAGE = 0;	/* Set if we are using a lot of memory */
				/* If user changes the available memory */
				/* settings, then we should stop and */
				/* restart our computations */

char ERRMSG0[] = "Iteration: %ld/%ld, %s";
char ERRMSG1A[] = "ERROR: ILLEGAL SUMOUT\n";
char ERRMSG1B[] = "ERROR: SUM(INPUTS) != SUM(OUTPUTS), %.16g != %.16g\n";
char ERRMSG1C[] = "ERROR: ROUND OFF (%.10g) > 0.40\n";
char ERRMSG2[] = "Possible hardware failure, consult the readme file.\n";
char ERRMSG3[] = "Continuing from last save file.\n";
char ERRMSG4[] = "Waiting five minutes before restarting.\n";
char ERROK[] = "Disregard last error.  Result is reproducible and thus not a hardware problem.\n";
char WRITEFILEERR[] = "Error writing intermediate file: %s\n";
char DONE_MSG1[] = "There are no more exponents to test.\n";
char DONE_MSG2[] = "Please send the results.txt file to woltman@alum.mit.edu\n";
char DONE_MSG3[] = "Contact the PrimeNet server for more exponents.\n";
char RENAME_MSG[] = "Renaming intermediate file %s to %s.\n";

void gwtobinary (gwnum, giant);
void binarytogw (giant,	gwnum);

/* Return true is exponent yields a known Mersenne prime */

int isKnownMersennePrime (
	unsigned long p)
{
	return (p == 2 || p == 3 || p == 5 || p == 7 || p == 13 || p == 17 ||
		p == 19 || p == 31 || p == 61 || p == 89 || p == 107 ||
		p == 127 || p == 521 || p == 607 || p == 1279 || p == 2203 ||
		p == 2281 || p == 3217 || p == 4253 || p == 4423 ||
		p == 9689 || p == 9941 || p == 11213 || p == 19937 ||
		p == 21701 || p == 23209 || p == 44497 || p == 86243 ||
		p == 110503 || p == 132049 || p == 216091 || p == 756839 ||
		p == 859433 || p == 1257787 || p == 1398269 || p == 2976221 ||
		p == 3021377 || p == 6972593);
}

/* Routine to set the stop timer when available memory settings change */

void memSettingsChanged (void) {
	if (HIGH_MEMORY_USAGE || STOP_TIME) {
		OutputStr ("Restarting with new memory settings.\n");
		STOP_TIME = 1;
	}
}

/* Return memory available now */
/* Compute the time where we must halt processing to change the */
/* amount of available memory. */ 

unsigned int avail_mem (void)
{
	time_t	t;
	struct tm *x;
	unsigned int curtime;

/* Set flag indicating we are using a lot of memory */

	HIGH_MEMORY_USAGE = TRUE;

/* If the same memory is available both day and night, then return */
/* that value and note that we won't have to stop in the future. */

	if (DAY_MEMORY == NIGHT_MEMORY) {
		STOP_TIME = 0;
		return (DAY_MEMORY);
	}

/* Determine whether it is daytime or nighttime. */
/* Return corresponding available memory. */
/* Set timer for when daytime or nighttime ends. */

	time (&t);
        x = localtime (&t);
	curtime = x->tm_hour * 60 + x->tm_min;
	if (DAY_START_TIME < DAY_END_TIME) {
		if (curtime < DAY_START_TIME) {
			STOP_TIME = t + (DAY_START_TIME - curtime) * 60;
			return (NIGHT_MEMORY);
		} else if (curtime < DAY_END_TIME) {
			STOP_TIME = t + (DAY_END_TIME - curtime) * 60;
			return (DAY_MEMORY);
		} else {
			STOP_TIME = t + (DAY_START_TIME + 1440 - curtime) * 60;
			return (NIGHT_MEMORY);
		}
	} else {
		if (curtime < DAY_END_TIME) {
			STOP_TIME = t + (DAY_END_TIME - curtime) * 60;
			return (DAY_MEMORY);
		} else if (curtime < DAY_START_TIME) {
			STOP_TIME = t + (DAY_START_TIME - curtime) * 60;
			return (NIGHT_MEMORY);
		} else {
			STOP_TIME = t + (DAY_END_TIME + 1440 - curtime) * 60;
			return (DAY_MEMORY);
		}
	}
}

/* Make a string out of a 96-bit value (a found factor) */

void makestr (
	unsigned long hsw,
	unsigned long msw,
	unsigned long lsw,
	char	*buf)			/* An 80 character output buffer */
{
	int	i, j, k, carry;
	unsigned long x[3];
	char	pow[80];

	x[0] = hsw; x[1] = msw; x[2] = lsw;
	for (i = 0; i < 79; i++) pow[i] = '0', buf[i] = '0';
	pow[78] = '1';
	pow[79] = buf[79] = 0;

	for (i = 3; i--; ) {
		for (j = 0; j < 32; j++) {
			if (x[i] & 1) {
				carry = 0;
				for (k = 79; k--; ) {
					buf[k] = buf[k] - '0' +
						pow[k] - '0' + carry;
					carry = buf[k] / 10;
					buf[k] %= 10;
					buf[k] += '0';
				}
			}
			carry = 0;
			for (k = 79; k--; ) {
				pow[k] = (pow[k] - '0') * 2 + carry;
				carry = pow[k] / 10;
				pow[k] %= 10;
				pow[k] += '0';
			}
			x[i] >>= 1;
		}
	}
	while (buf[0] == '0') strcpy (buf, buf+1);
}

/* Generate the 64-bit residue of a Lucas-Lehmer test */

void generateResidue64 (
	unsigned long units_bit,
	unsigned long *reshi,
	unsigned long *reslo)
{
	unsigned long i, word, bit_in_word, high32, low32;
	long	val;
	int	j, bits, bitsout, carry;

/* Find out where the least significant bit is */

	bitaddr (units_bit, &word, &bit_in_word);

/* Check for a carry out of the previous word */

	for (i = word; ; ) {
		if (i == 0) i = FFTLEN;
		i--;
		get_fft_value (LLDATA, i, &val);
		if (val != 0 || i == word) {
			carry = (val >= 0) ? 0 : -1;
			break;
		}
	}

/* Collect bits until we have 64 of them */

	high32 = low32 = 0;
	bitsout = 0;
	for (i = word; bitsout < 64; i = (i + 1) % FFTLEN) {
		get_fft_value (LLDATA, i, &val);
		val += carry;
		bits = (int) BITS_PER_WORD;
		if (is_big_word (i)) bits++;
		if (bit_in_word) {
			val >>= bit_in_word;
			bits -= (int) bit_in_word;
			bit_in_word = 0;
		}
		for (j = 0; j < bits && bitsout < 64; j++) {
			low32 >>= 1;
			if (high32 & 1) low32 += 0x80000000;
			high32 >>= 1;
			if (val & 1) high32 += 0x80000000;
			val >>= 1;
			bitsout++;
		}
		carry = (int) val;
	}

/* Return the result */

	*reshi = high32;
	*reslo = low32;
}

/* Return TRUE if a continuation file exists.  If one does exist, */
/* make sure it is named pXXXXXXX. */

int continuationFileExists (
	char	*filename)
{
	char	backupname[32];
	char	buf[80];

	if (fileExists (filename)) return (TRUE);
	strcpy (backupname, filename);
	backupname[0] = 'q';
	if (fileExists (backupname)) {
		sprintf (buf, RENAME_MSG, backupname, filename);
		OutputBoth (buf);
		rename (backupname, filename);
		return (TRUE);
	}
	backupname[0] = 'r';
	if (fileExists (backupname)) {
		sprintf (buf, RENAME_MSG, backupname, filename);
		OutputBoth (buf);
		rename (backupname, filename);
		return (TRUE);
	}
	return (FALSE);
}

/* Write intermediate Lucas-Lehmer results to a file */
/* Note: we can't output the bits == 16 case in two bytes because the */
/* rounding code in lucas.asm can create a value of 32768 which does*/
/* not fit in a short. */

int writeToFile (
	char	*filename,
	unsigned long counter,
	unsigned long units_bit,
	unsigned long error_count)
{
	int	fd, bits;
	unsigned long i;
	char	buf[4096];
	char	*bufp;
	short	type;
	unsigned long buggy_error_count = 0;
	long	sum = 0;

/* If we are allowed to create multiple intermediate files, then */
/* write to a file called rXXXXXXX. */

	if (TWO_BACKUP_FILES && strlen (filename) == 8)
		filename[0] = 'r';

/* Now save to the intermediate file */

	fd = _open (filename, _O_BINARY | _O_WRONLY | _O_TRUNC | _O_CREAT, 0666);
	if (fd < 0) return (FALSE);
	if (FFTLEN < 8192)
		type = (short) FFTLEN + 1;
	else
		type = (short) (FFTLEN / 1024);
	if (_write (fd, &type, sizeof (short)) != sizeof (short))
		goto writeerr;
	if (_write (fd, &counter, sizeof (long)) != sizeof (long))
		goto writeerr;

	bits = (int) BITS_PER_WORD + 1;
	bufp = buf;
	for (i = 0; i < FFTLEN; i++) {
		long	x;
		if (bufp - buf >= sizeof (buf) - 20) {
			if (_write (fd, buf, bufp - buf) != bufp - buf)
				goto writeerr;
			bufp = buf;
		}
		get_fft_value (LLDATA, i, &x);
		if (bits <= 15) {
			short y;
			y = (short) x;
			memcpy (bufp, &y, sizeof (y));
			bufp += sizeof (y);
  		} else {
			long y;
			y = (long) x;
			memcpy (bufp, &y, sizeof (y));
			bufp += sizeof (y);
		}
		sum += x;
	}
	sum += units_bit;
	sum += error_count;
	/* Kludge so that buggy v17 save files are rejected */
	if (units_bit != 0) sum ^= 0x1;
	memcpy (bufp, &sum, sizeof (long));
	bufp += sizeof (long);
	memcpy (bufp, &units_bit, sizeof (unsigned long));
	bufp += sizeof (unsigned long);
	memcpy (bufp, &buggy_error_count, sizeof (unsigned long));
	bufp += sizeof (unsigned long);
	memcpy (bufp, &error_count, sizeof (unsigned long));
	bufp += sizeof (unsigned long);
	if (_write (fd, buf, bufp - buf) != bufp - buf) goto writeerr;
	_commit (fd);
	_close (fd);

/* Now rename the intermediate files */

	if (TWO_BACKUP_FILES && strlen (filename) == 8) {
		char	backupname[16];
		strcpy (backupname, filename);
		backupname[0] = 'q'; filename[0] = 'p';
		_unlink (backupname);
		 rename (filename, backupname);
		backupname[0] = 'r';
		rename (backupname, filename);
	}

	return (TRUE);

/* An error occured.  Delete the current file and rename the backup */
/* intermediate file */

writeerr:
	_close (fd);
	_unlink (filename);
	return (FALSE);
}

/* Update the error count at the end of an intermediate file */
/* Note: If a save file hasn't been created yet or the save file is */
/* a pre-V19 save file, then the error will not get recorded. */

void writeNewErrorCount (
	char	*filename,
	unsigned long new_error_count)
{
	int	fd;
	long	offset, units_bit, sum, trash;
	unsigned long old_error_count;

/* Open the intermediate file, position past the FFT data */

	fd = _open (filename, _O_BINARY | _O_RDWR);
	if (fd < 0) return;

	offset = sizeof (short) +			/* type field */
		 sizeof (long) +			/* counter field */
		 FFTLEN * ((BITS_PER_WORD + 1 <= 15) ?	/* FFT data */
			sizeof (short) : sizeof (long));
	
/* Read in the checksum, units_bit, and old error count */

	_lseek (fd, offset, SEEK_SET);
	if (_read (fd, &sum, sizeof (long)) != sizeof (long))
		goto err;
	if (_read (fd, &units_bit, sizeof (long)) != sizeof (long))
		goto err;
	if (_read (fd, &trash, sizeof (long)) != sizeof (long))
		goto err;
	if (_read (fd, &old_error_count, sizeof (long)) != sizeof (long))
		goto err;

/* Update the checksum */

	if (units_bit != 0) sum ^= 0x1;
	sum = sum - old_error_count + new_error_count;
	if (units_bit != 0) sum ^= 0x1;

/* Write out the checksum and new error count */

	_lseek (fd, offset, SEEK_SET);
	if (_write (fd, &sum, sizeof (long)) != sizeof (long))
		goto err;
	_lseek (fd, offset + 3 * sizeof (long), SEEK_SET);
	if (_write (fd, &new_error_count, sizeof (long)) != sizeof (long))
		goto err;

/* Close file and return */

err:	_close (fd);
}

/* Read the data portion of an intermediate Lucas-Lehmer results file */

int readFileData (
	int	fd,
	unsigned long *units_bit,
	unsigned long *error_count)
{
	unsigned long i, buggy_error_count;
	long	sum, filesum;
	int	bits, zero;		/* Guard against a zeroed out file */

	if (fd < 0) return (FALSE);

	bits = (int) BITS_PER_WORD + 1;
	sum = 0;
	zero = TRUE;
	for (i = 0; i < FFTLEN; i++) {
		long	x;
		if (bits <= 15) {
			short y;
			if (_read (fd, &y, sizeof (y)) != sizeof (y))
				goto err;
			x = y;
		} else {
			long y;
			if (_read (fd, &y, sizeof (y)) != sizeof (y))
				goto err;
			x = y;
		}
		sum += x;
		if (x) zero = FALSE;
		set_fft_value (LLDATA, i, x);
	}
	if (_read (fd, &filesum, sizeof (long)) != sizeof (long)) goto err;
	if (_read (fd, units_bit, sizeof (long)) != sizeof (long))
		*units_bit = 0;
	else
		sum += *units_bit;
	/* V18 and earlier did not save a correct error count */
	/* Read in this erroneous error count */
	if (_read (fd, &buggy_error_count, sizeof (long)) != sizeof (long))
		buggy_error_count = 0;
	else
		sum += buggy_error_count;
	/* Now read in the correct V19 error count */
	if (_read (fd, error_count, sizeof (long)) != sizeof (long))
		*error_count = 0x80000000 | buggy_error_count;
	else
		sum += *error_count;
	/* Kludge so that buggy v17 save files are rejected */
	/* V18 and later flip the bottom checksum bit */
	if (*units_bit != 0) sum ^= 0x1;
	/* Kludge so that only large exponent v17 saves are rejected */
	/* Clear bottom checksum bits if v17 would not corrupt this save file */
	if (PARG < 4194304) sum &= 0xFFFFFFFE, filesum &= 0xFFFFFFFE;
	if (filesum != sum) goto err;
	if (zero) goto err;
	_close (fd);
	return (TRUE);
err:	_close (fd);
	return (FALSE);
}

/* Update the work-to-do list.  This isn't perfect code, but should */
/* handle all the cases where the worktodo file wasn't bizarrely edited. */

void updateWorkToDo (
	unsigned long exp,	/* exponent that was worked on */
	int	work_type,	/* type of work performed */
	unsigned long data)	/* how_far_factored OR ECM bound */
{
	unsigned int i;


/* Note if the well-behaved-work-option is on, then just */
/* delete the first line and write the file every half hour */

	if (WELL_BEHAVED_WORK) {
		static time_t last_time_written = 0;
		time_t	current_time;
		IniDeleteLine (WORKTODO_FILE, 1);
		time (&current_time);
		if (current_time > last_time_written + 1800) {
			if (last_time_written) IniFileClose (WORKTODO_FILE);
			last_time_written = current_time;
		}
		return;
	}

/* Read the WorkToDo from disk prior to updating it. */

	IniFileOpen (WORKTODO_FILE, 0);

/* Clear the exponent from the work-to-do-file */

	for (i = 1; ; i++) {
		struct work_unit w;
		int	changed = FALSE, deleteLine = FALSE, done = FALSE;

/* Read the line of the work file */

		if (! parseWorkToDoLine (i, &w)) break;

/* Skip the line if the exponent does not match */

		if (exp != w.p) continue;

/* Switch off the work type just completed */

		switch (work_type) {

/* Delete line if a factor was found or this was a factoring-only work type. */
/* Otherwise update Test, Dblchk, Pfactor lines. */

		case WORK_FACTOR:
			if (data == 0 || w.work_type == WORK_FACTOR)
				deleteLine = TRUE;
			else if (w.work_type == WORK_TEST ||
				 w.work_type == WORK_DBLCHK ||
				 w.work_type == WORK_PFACTOR) {
				w.bits = data;
				changed = TRUE;
			}
			break;

/* Delete line if an LL test completed. */

		case WORK_TEST:
			deleteLine = TRUE;
			break;

/* Delete matching ECM line */

		case WORK_ECM:
			if (w.work_type == WORK_ECM &&
			    w.B1 == data && (int) PLUS1 == w.plus1) {
				deleteLine = TRUE;
				done = TRUE;
			}
			break;

/* Delete matching Pminus1 and Pfactor lines, update Test and DblChk lines. */

		case WORK_PMINUS1:
			if ((w.work_type == WORK_PMINUS1 && (int) PLUS1 == w.plus1) ||
			    w.work_type == WORK_PFACTOR) {
				deleteLine = TRUE;
				done = TRUE;
			} else if (w.work_type == WORK_TEST ||
				   w.work_type == WORK_DBLCHK) {
				w.pminus1ed = TRUE;
				changed = TRUE;
			}
			break;
		}

/* Now implement the decisions made above */

		if (changed) {
			char	buf[20];
			sprintf (buf, "%ld,%d,%d", w.p, w.bits, w.pminus1ed);
			IniReplaceLineAsString (
				WORKTODO_FILE, i,
				(w.work_type == WORK_DBLCHK) ? "DoubleCheck" :
				(w.work_type == WORK_TEST) ? "Test" : "Pfactor",
				buf);
		}

		if (deleteLine) {
			IniDeleteLine (WORKTODO_FILE, i);
			i--;
		}

		if (done) break;
	}

/* Write out the updated worktodo file */

	IniFileClose (WORKTODO_FILE);
}

/* Returns true if this is a priority work item */

int isPriorityWork (
	struct work_unit *w)
{
	if (w->p == EXP_BEING_WORKED_ON) return (FALSE);
	if (w->work_type == WORK_ADVANCEDTEST) return (TRUE);
	if (IniGetInt (INI_FILE, "SequentialWorkToDo", 1)) return (FALSE);
	if ((w->work_type == WORK_TEST ||
	     w->work_type == WORK_DBLCHK) &&
	    (w->bits < factorLimit (w->p, 0) || !w->pminus1ed)) return (TRUE);
	return (FALSE);
}

/* Check if any of the Lucas-Lehmer test lines also require factoring. */
/* This will force factoring to be done first - giving us more accurate */
/* estimates of how much work is queued up. */

int getPriorityWork (void) {
	unsigned int i;
	struct work_unit w;

	IniFileOpen (WORKTODO_FILE, 0);
	for (i = 1; ; i++) {
		if (!parseWorkToDoLine (i, &w)) break;
		if (isPriorityWork (&w)) return (i);
	}

/* No priority work to do */

	return (0);
}

/* Clear rolling average start time in case Factoring, Torture Test, */
/* Advanced/Test, etc. interrupted a Lucas-Lehmer test. */

void clearRollingStart (void)
{
	IniWriteInt (LOCALINI_FILE, "RollingStartTime", 0);
}
		
/* We're starting a prime, save the current time so that */
/* we can maintain the rolling average computations. */

void startRollingAverage (void)
{
	time_t	current_time;
	time (&current_time);
	if (VACATION_END && !ON_DURING_VACATION) current_time = 0;
	IniWriteInt (LOCALINI_FILE, "RollingStartTime", current_time);
}

/* We've completed 65536 Lucas-Lehmer iterations or a factoring. */
/* Examine how long it took compared to our expectaions.  Maintain */
/* a rolling average for more accurate time estimates. */

void updateRollingAverage (
	double	est)			/* estimated time for last work unit */
{
#ifndef SERVER_TESTING
	time_t	start_time, current_time;

/* Get the current time and when the iterations started */

	time (&current_time);
	start_time = IniGetInt (LOCALINI_FILE, "RollingStartTime", 0);
	if (start_time == 0 || current_time <= start_time) return;

/* Compute the rolling average where the existing rolling average
/* accounts for 90% of the new rolling average and the current */
/* data point accounts for 10%.  The rolling average */
/* measures this computer's speed as compared to our estimated */
/* speed for this CPU.  A value of 1200 means the this computer */
/* is 20% faster than expected.  A value of 1000 means the this computer */
/* is as fast as expected.   A value of 800 means this computer */
/* is 20% slower than expected. */

/* The old formula used prior to version 19.1 (ra = 0.9 * ra + 0.1 * est) */
/* produced poor results on machines that were not on 24 hours a day.  As */
/* an example the wild fluctuations of 2000, 667, 2000, 667, etc. averaged */
/* out to well more than the correct 1000.  Brian Beesley suggested a new */
/* formula (ra = 1 / (0.9 / ra + 0.1 / est)) which works much better. */

	est = est / (current_time - start_time) * 1000.0;
	if (est > 50000.0) return;	/* A safeguard against bogus data */
	if (est < 0.5 * ROLLING_AVERAGE) est = 0.5 * ROLLING_AVERAGE;
	if (est > 2.0 * ROLLING_AVERAGE) est = 2.0 * ROLLING_AVERAGE;

	ROLLING_AVERAGE = (unsigned int)
		(1.0 / (0.9 / ROLLING_AVERAGE + 0.1 / est) + 0.5);
	if (ROLLING_AVERAGE < 20) ROLLING_AVERAGE = 20;
	if (ROLLING_AVERAGE > 4000) ROLLING_AVERAGE = 4000;
	IniWriteInt (LOCALINI_FILE, "RollingAverage", ROLLING_AVERAGE);
#endif
}

/* Utility output routines */

void LineFeed (void)
{
	OutputStr ("\n");
}

void OutputTimeStamp ()
{
	time_t	this_time;
	char	tbuf[40], buf[40];

	if (TIMESTAMPING) {
		time (&this_time);
		strcpy (tbuf, ctime (&this_time)+4);
		tbuf[12] = 0;
		sprintf (buf, "[%s] ", tbuf);
		OutputStr (buf);
	}
}

void OutputNum (
	unsigned long num)
{
	char	buf[20];
	sprintf (buf, "%lu", num);
	OutputStr (buf);
}

int stopCheck (void)
{
static	int	time_counter = 0;

/* If the ESC key is hit, stop processing */

	if (escapeCheck ()) {
		STOP_REASON = STOP_ESCAPE;
		return (TRUE);
	}

/* Perform timed checks less often */

	if (++time_counter <= 100) return (FALSE);
	time_counter = 0;

/* If this is a timed run, check if our time has expired */
/* This is a little used feature described in the readme.txt file */

	if (END_TIME) {
		time_t	this_time;
		time (&this_time);
		if (this_time > END_TIME) {
			STOP_REASON = STOP_TIMEOUT;
			return (TRUE);
		}
	}

/* If this is a timed run, check if our time has expired */
/* This is a little used feature described in the readme.txt file */

	if (STOP_TIME) {
		time_t	this_time;
		time (&this_time);
		if (this_time > STOP_TIME) {
			STOP_REASON = STOP_MEM_CHANGED;
			return (TRUE);
		}
	}

/* No need to stop */

	return (FALSE);
}

/* Sleep five minutes before restarting */

int SleepFive (void)
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

void clear_timers (void) {
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
	if (RDTSC_TIMING && (CPU_FLAGS & CPU_RDTSC)) {
		unsigned long hi, lo;
		rdtsc (&hi, &lo);
		timers[i] -= (double) hi * 4294967296.0 + lo;
	} else {
		struct _timeb timeval;
		_ftime (&timeval);
		timers[i] -= (double) timeval.time * 1000.0 + timeval.millitm;
	}
}

void end_timer (
	int	i)
{
	if (RDTSC_TIMING && (CPU_FLAGS & CPU_RDTSC)) {
		unsigned long hi, lo;
		rdtsc (&hi, &lo);
		timers[i] += (double) hi * 4294967296.0 + lo;
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
	if (RDTSC_TIMING && (CPU_FLAGS & CPU_RDTSC))
		return (timers[i] / CPU_SPEED / 1000000.0);
	else
		return (timers[i] / 1000.0);
}

#define TIMER_NL	0x1
#define TIMER_CLR	0x2
#define TIMER_OPT_CLR	0x4
#define TIMER_MS	0x8
#define TIMER_OUT_BOTH  0x1000

void print_timer (
	int	i,
	int	flags)
{
	double	t;
	int	style;
	char	buf[40];

/* Format the timer value in one of several styles */

	t = timer_value (i);
	style = IniGetInt (INI_FILE, "TimingOutput", 0);
	if (style == 0) {
		if (flags & TIMER_MS) style = 4;
		else style = 1;
	}
	if (style == 1)
		sprintf (buf, "%.3f sec.", t);
	else if (style == 2)
		sprintf (buf, "%.1f ms.", t * 1000.0);
	else if (style == 3)
		sprintf (buf, "%.2f ms.", t * 1000.0);
	else
		sprintf (buf, "%.3f ms.", t * 1000.0);
	if (RDTSC_TIMING == 2 && (CPU_FLAGS & CPU_RDTSC)) {
		sprintf (buf+strlen(buf), " (%.0f clocks)", timers[i]);
	}
	if (flags & TIMER_NL) strcat (buf, "\n");

/* Clear the timer */

	if (flags & TIMER_CLR) timers[i] = 0.0;
	if ((flags & TIMER_OPT_CLR) && !CUMULATIVE_TIMING) timers[i] = 0.0;

/* Output the formatted timer value */

	if (flags & TIMER_OUT_BOTH)
		OutputBoth (buf);
	else
		OutputStr (buf);
}

/* More utility routines */

/* Prepare for making a factoring run */

void factorSetup (void)
{

/* Allocate 1MB memory for factoring */

	SRCARG = FACDATA = malloc (1000000);
}

/* Cleanup after making a factoring run */

void factorDone (void)
{

/* Free factoring data */

	free (FACDATA);
}

/* Prepare for running a Lucas-Lehmer test */

void lucasSetup (
	unsigned long p,	/* Exponent to test */
	unsigned long fftlen)	/* Specific FFT length to use, or zero */
{

/* Init the FFT code */

	gwsetup (p, fftlen, 0);

/* Allocate memory for the Lucas-Lehmer data (the number to square) */

	LLDATA = gwalloc ();
}

/* Clean up after running a Lucas-Lehmer test */

void lucasDone (void)
{

/* Free memory for the Lucas-Lehmer data */

	gwfree (LLDATA);

/* Cleanup the FFT code */

	gwdone ();
}

/* Continue factoring/testing Mersenne numbers */

void primeContinue (void) 
{
	struct work_unit w;
	char	filename[20];
	unsigned int line, pass, sequential;
	int	fd;
	short	type;
	unsigned long counter, fftlen;
static	int	first_call = TRUE;

/* Output a message indicating we are starting an LL test */

	if (first_call) {
		char	buf[80];
		first_call = FALSE;
		sprintf (buf, "Mersenne number primality test program version %s\n", VERSION);
		OutputStr (buf);
	}

/* See if we are in a sleep time period */
/* Sleep in two second intervals so we can detect a thread stop */
/* or program exit command. */

again:	if (SLEEP_TIME) {
		time_t	current_time;
		char	buf[30];
		sprintf (buf, "Sleeping until %s\n", ctime (&SLEEP_TIME));
		OutputStr (buf);
		title ("Sleeping");
		ChangeIcon (IDLE_ICON);
		for ( ; ; ) {
			time (&current_time);
			if (current_time >= SLEEP_TIME) break;
			Sleep (2000);
			if (escapeCheck ()) return;
		}
		ChangeIcon (WORKING_ICON);
		readIniFiles ();
		goto again;
	}

/* If the ESC key is hit (or we shouldn't even get started due to */
/* a laptop running on a battery) then return. */

	if (escapeCheck ()) return;

/* Set the process/thread priority */

	SetPriority ();

/* Send new completion dates once a month */

	ConditionallyUpdateEndDates ();

/* Every time the user chooses Test/Continue, clear the timer that */
/* prevents communication for a period of time.  This allows the user */
/* to try something and if it doesn't work, ESC and choose Test/Continue */
/* to try some other system settings (without waiting an hour). */

	next_comm_time = 0; 

/* Loop until the ESC key is hit or the entire work-to-do INI file */
/* is processed and we are not connected to the server. */

	sequential = IniGetInt (INI_FILE, "SequentialWorkToDo", 1);
	for ( ; ; ) {

/* Clear variable indicating reason for stopping */

	STOP_REASON = 0;

/* Send any queued results and get work to do (if necessary) */

	CHECK_WORK_QUEUE = 1;
	if (! communicateWithServer ()) goto check_stop_code;

/* If the ESC key is hit, break out of this loop */

	if (escapeCheck ()) goto check_stop_code;

/* See if there is any work queued up for us to do. */
/* If not, end the thread if we're not connected to primenet. */
/* Otherwise loop until we get some work from primenet. */

	IniFileOpen (WORKTODO_FILE, 0);
	if (IniGetNumLines (WORKTODO_FILE) == 0) {
		ChangeIcon (IDLE_ICON);
		if (! USE_PRIMENET) {
			OutputSomewhere (DONE_MSG1);
			OutputSomewhere (DONE_MSG2);
			OutputSomewhere (DONE_MSG3);
			goto done;
		}
		title ("IDLE");
		Sleep (1000);
		ChangeIcon (WORKING_ICON);
		continue;
	}

/* Clear timer indicating we need to stop and reprocess the worktodo.ini */
/* file because more or less memory is now available. */

	STOP_TIME = 0;

/* Make three passes over the worktodo.ini file looking for the ideal */
/* piece of work to do.  In pass 1, we look for high-priority work.  This */
/* includes trial and P-1 factoring prior to an LL test.  If a factor is */
/* found, it can reduce the amount of work we have queued up, requiring */
/* us to ask the server for more.  In pass 2, we process the file in */
/* order (except for LL tests that are not yet ready because the P-1 */
/* factoring has not completed).  In pass 3, as a last resort we start */
/* LL tests where P-1 factoring is stalled because of low memory. */

/* Skip one pass on large well-behaved work files */

	for (pass = WELL_BEHAVED_WORK ? 2 : 1; pass <= 3; pass++) {

/* Examine each line in the worktodo.ini file */

	for (line = 1; ; line++) {

/* Clear flag indicating we are using a lot of memory */

	HIGH_MEMORY_USAGE = FALSE;

/* Read the line from the work file, break when out of lines */

	if (!parseWorkToDoLine (line, &w)) break;

/* Make sure this line of work from the file makes sense. The exponent */
/* should be a prime number, bounded by values we can handle, and we */
/* should never be asked to factor a number more than we are capable of. */

	if (pass == 1 && !isPrime (w.p) &&
	    w.work_type != WORK_ECM && w.work_type != WORK_PMINUS1 &&
	    w.work_type != WORK_ADVANCEDFACTOR) {
		char	buf[80];
		sprintf (buf, "Error: Work-to-do file contained composite exponent: %ld\n", w.p);
		LogMsg (buf);
		updateWorkToDo (w.p, WORK_FACTOR, 0);
		goto did_some_work;
	}
	if ((w.work_type == WORK_TEST || w.work_type == WORK_DBLCHK) &&
	    (w.p < MIN_PRIME || w.p > MAX_PRIME)) {
		char	buf[80];
		sprintf (buf, "Error: Work-to-do file contained bad LL exponent: %ld\n", w.p);
		LogMsg (buf);
		updateWorkToDo (w.p, WORK_FACTOR, 0);
		goto did_some_work;
	}
	if (w.work_type == WORK_FACTOR &&
	    (w.p < 727 || w.p > MAX_FACTOR || w.bits >= factorLimit (w.p, 1))) {
		char	buf[100];
		sprintf (buf, "Error: Work-to-do file contained bad factoring assignment: %ld,%d\n", w.p, w.bits);
		LogMsg (buf);
		updateWorkToDo (w.p, WORK_FACTOR, 0);
		goto did_some_work;
	}

/* Does a continuation file exist?  If so, open and validate it, setting */
/* fftlen if this is an LL save file.  Otherwise, set fd to zero. */

readloop:
	tempFileName (filename, w.p);
	fd = 0;
	fftlen = 0;
	counter = 0;
	if (continuationFileExists (filename)) {

/* Read the initial data from the file */

		if (! readFileHeader (filename, &fd, &type, &counter)) {
			char	buf[80];
readerr:		sprintf (buf, READFILEERR, filename);
			OutputBoth (buf);
			_unlink (filename);
			goto readloop;
		}

/* Type 2 files are factoring continuation files */
/* Type 3 files are Advanced / Factoring continuation files */

		if (type != 2 && type != 3) {

/* Deduce the fftlen from the type field */

	 		if (type & 1)
				fftlen = (unsigned long) type - 1;
			else	
				fftlen = (unsigned long) type * 1024;

/* Check for corrupt LL continuation files. */

			if ((type & 1 &&
			     fftlen != map_exponent_to_fftlen (w.p, GW_MERSENNE_MOD)) ||
			    counter > w.p) {
				_close (fd);
				goto readerr;
			}
		}
	}

/* See if the exponent needs more factoring.  We treat factoring that is */
/* part of a pfactor or LL test as priority work (done in pass 1).  If the */
/* save file indicates the LL test has begun, then don't do the factoring */
/* (it would destroy the LL save file!) */

	if (((w.work_type == WORK_FACTOR && pass == 2) ||
	     ((w.work_type == WORK_PFACTOR ||
	       w.work_type == WORK_TEST ||
	       w.work_type == WORK_DBLCHK) &&
	      ((pass == 1 && !sequential) || pass == 2) &&
	      ! IniGetInt (INI_FILE, "SkipTrialFactoring", 0))) &&
	    (fd == 0 || fftlen == 0) &&
	    w.bits < factorLimit (w.p, w.work_type == WORK_FACTOR)) {
		int	res;
		if (! primeFactor (w.p, w.bits, &res, w.work_type, fd))
			goto check_stop_code;
		if (res == 999) goto readerr;
		if (w.work_type == WORK_FACTOR && WELL_BEHAVED_WORK) {
			line--;
			continue;
		}
		goto did_some_work;
	}

/* Make sure the first-time user runs a successful self-test. */

	if (((pass == 1 && !sequential && !w.pminus1ed) || pass == 2) &&
	    (w.work_type == WORK_TEST ||
	     w.work_type == WORK_DBLCHK || w.work_type == WORK_PFACTOR)) {
		if (! selfTest (w.p)) goto check_stop_code;
	}

/* See if this exponent needs special P-1 factoring.  We treat P-1 factoring */
/* that is part of an LL test as priority work (done in pass 1).  If the */
/* save file indicates an LL test is more than 50% complete, then don't */
/* do the factoring - it would likely annoy the user. */

	if (((w.work_type == WORK_PFACTOR && pass == 2) ||
	     ((w.work_type == WORK_TEST || w.work_type == WORK_DBLCHK) &&
	      ! w.pminus1ed &&
	      ((pass == 1 && !sequential) || pass == 2) &&
	      (fd == 0 || fftlen == 0 || counter < w.p / 2)))) {
		if (fd) _close (fd);
		if (! pfactor (&w)) {
			if (STOP_REASON != STOP_NOT_ENOUGH_MEM) goto check_stop_code;
		} else
			goto did_some_work;
	}

/* Run the LL test */

	if ((pass == 1 && w.work_type == WORK_ADVANCEDTEST) ||
	    (((pass == 2 && (w.pminus1ed || counter >= w.p / 2)) || pass == 3) &&
	     (w.work_type == WORK_TEST || w.work_type == WORK_DBLCHK))) {
		unsigned long units_bit, err_cnt;

/* Setup */

		lucasSetup (w.p, fftlen);

/* Read the initial data from the file, on failure try the backup */
/* intermediate file. */

		if (fd && fftlen) {
			if (! readFileData (fd, &units_bit, &err_cnt)) {
				lucasDone ();
				goto readerr;
			}

/* Handle case where the save file was for a different FFT length than */
/* we would prefer to use.  This can happen, for example, when upgrading */
/* to a Pentium 4 with its SSE2-based FFT using less precicion. */

			if (fftlen != map_exponent_to_fftlen (w.p, 0)) {
				giant	g;
				g = newgiant ((w.p + 32) / sizeof (short));
				gwtobinary (LLDATA, g);
				lucasDone ();
				lucasSetup (w.p, 0);
				binarytogw (g, LLDATA);
				free (g);
			}
		}

/* Start off with the 1st Lucas number */

		else {
			counter = 2;
			units_bit = 0;
			err_cnt = 0;
		}

/* We read the continuation file successfully, continue LL testing */

		if (! prime (w.p, counter, units_bit, err_cnt)) {
			lucasDone ();
			goto check_stop_code;
		}
		lucasDone ();
		goto did_some_work;
	}

/* Check for Advanced / Factoring work */

	if (w.work_type == WORK_ADVANCEDFACTOR && pass == 2) {
		if (! primeSieve (w.bits, w.B1, (unsigned short) w.B2_start,
				  (unsigned short) w.B2_end, fd))
			goto check_stop_code;
		goto did_some_work;
	}

/* Close the temporary file if we haven't used it by now */

	if (fd) _close (fd);

/* See if this is an ECM factoring line */

	if (w.work_type == WORK_ECM && pass == 2) {
		if (! ecm (w.p, w.B1, w.B2_start, w.B2_end, w.curves_to_do,
			   w.curves_completed, w.curve, w.plus1))
			goto check_stop_code;
		goto did_some_work;
	}

/* See if this is an P-1 factoring line */

	if (w.work_type == WORK_PMINUS1 && pass == 2) {
		if (! pminus1 (w.p, w.B1, w.B2_start, w.B2_end, w.plus1, FALSE))
			goto check_stop_code;
		goto did_some_work;
	}

/* This work-to-do line wasn't processed by any of the cases above */
/* Move onto the next worktodo line */

	}

/* Make another pass over the worktodo.ini file */

	}

/* Ugh, we made three passes over the worktodo file and couldn't find */
/* any work to do.  I think this can only happen if we are low on memory. */
/* If STOP_TIME is set, sleep until the timer expires. */

	if (STOP_TIME) {
		time_t	current_time;
		char	buf[30];
		sprintf (buf, "Sleeping until %s\n", ctime (&STOP_TIME));
		OutputStr (buf);
		title ("Sleeping");
		ChangeIcon (IDLE_ICON);
		for ( ; ; ) {
			time (&current_time);
			if (current_time >= STOP_TIME) break;
			Sleep (2000);
			if (escapeCheck ()) return;
		}
		ChangeIcon (WORKING_ICON);
	}

/* Loop to potentially contact the server for more work and then process */
/* the updated work-to-do file. */

did_some_work:;
	}

/* Examine the stop reason (if any) and either return, reread the INI */
/* files, or examine the work-to-do file for more work. */

check_stop_code:

/* If this was a timed run, read the new INI file settings */
/* and continue processing */

	if (STOP_REASON == STOP_TIMEOUT) {
		readIniFiles ();
		goto again;
	}

/* If we stopped because the available memory changed, then reprocess */
/* the worktodo file from the beginning in case there are items that */
/* were stalled waiting for more memory. */

	if (STOP_REASON == STOP_MEM_CHANGED) goto did_some_work;

/* Otherwise the escape key was hit, close the worktodo file in case */
/* WELL_BEHAVED_WORK caused us to delay writing the file, then return */

done:	IniFileClose (WORKTODO_FILE);
}

/* Increment the error counter.  The error counter is one 32-bit */
/* field that contains 5 values - a flag if this is a contiuation */
/* from a save file that did not track error counts, a count of */
/* errors that were reproducible, a count of ILLEAL SUMOUTs, */
/* a count of convolution errors above 0.4, and a count of */
/* SUMOUTs not close enough to SUMINPs. */

void inc_error_count (
	int	type,
	unsigned long *error_count)
{
	unsigned long addin, maxval, temp;
	
	addin = 1 << (type * 8);
	maxval = ((type == 3) ? 127 : 255) * addin;
	temp = *error_count & maxval;
	if (temp != maxval) temp += addin;
	*error_count = (*error_count & ~maxval) + temp;
}

/* Prepare for subtracting 2 from the squared result.  Also keep track */
/* of the location of the ever changing units bit. */

void lucas_fixup (
	unsigned long *units_bit)
{
	unsigned long sub_bit, word, bit_in_word;
	long	addin;

/* We are about to square the number, the units bit position will double */

	*units_bit <<= 1;
	if (*units_bit >= PARG) *units_bit -= PARG;

/* We will subtract two using the new location of the units bit */
	
	sub_bit = *units_bit + 1;
	if (sub_bit >= PARG) sub_bit -= PARG;
	bitaddr (sub_bit, &word, &bit_in_word);

/* Tell gwnum code the value to add to the squared result. */

	addin = -(1L << bit_in_word);
	gwsetaddin (word, addin);
}

/* Do the Lucas-Lehmer test */

int prime (
	unsigned long p,
	unsigned long counter,
	unsigned long units_bit,
	unsigned long error_count)
{
	long	write_time = DISK_WRITE_TIME * 60;
	unsigned long iters;
	char	filename[20];
	double	reallyminerr = 1.0;
	double	reallymaxerr = 0.0;
	int	priority_work = 0;
	int	escaped = 0;
	int	saving;
	unsigned long i, high32, low32;
	int	isPrime, rc;
	char	buf[160];
	time_t	start_time, current_time;
	unsigned long interimFiles, interimResidues;
static	unsigned long last_counter = 0;		/* Iteration of last error */

/* A new option to create interim save files every N iterations. */
/* This allows two machines to simultanously work on the same exponent */
/* and compare results along the way. */

	interimFiles = IniGetInt (INI_FILE, "InterimFiles", 0);
	interimResidues = IniGetInt (INI_FILE, "InterimResidues", interimFiles);

/* Clear all timers */

	clear_timers ();

/* Get the current time */

	time (&start_time);

/* Init filename */

	tempFileName (filename, p);

/* Init the title */

	sprintf (buf, "%ld / %ld", counter, p);
	title (buf);

/* Init global vars for Test/Status and CommunicateWithServer */

	EXP_BEING_WORKED_ON = p;
	EXP_BEING_FACTORED = 0;
	EXP_PERCENT_COMPLETE = (double) counter / (double) p;

/* Start off with the 1st Lucas number - four */
/* Note we do something a little strange here.  We actually set the */
/* first number to 4 but shifted by a random amount.  This lets two */
/* different machines check the same Mersenne number and operate */
/* on different FFT data - thus greatly reducing the chance that */
/* a CPU or program error corrupts the results. */

	if (counter == 2) {
		unsigned long i, word, bit_in_word, hi, lo;
		srand ((unsigned) time (NULL));
		units_bit = (rand () << 16) + rand ();
		if (CPU_FLAGS & CPU_RDTSC) { rdtsc(&hi,&lo); units_bit += lo; }
		units_bit = units_bit % p;
		bitaddr ((units_bit + 2) % p, &word, &bit_in_word);
		for (i = 0; i < FFTLEN; i++) {
			set_fft_value (LLDATA, i, (i == word) ? (1L << bit_in_word) : 0);
		}
		startRollingAverage ();

/* Immediately create a save file so that writeNewErrorCount can properly */
/* keep track of error counts. */

		if (p > 1500000) start_time = 0;

/* Output a message indicating we are starting an LL test */

		sprintf (buf, "Starting primality test of M%ld\n", p);
		OutputStr (buf);
	}

/* Otherwise, output a message indicating we are resuming an LL test */

	else {
		char	fmt_mask[80];
		double	pct;
		pct = trunc_percent (counter * 100.0 / p);
		sprintf (fmt_mask,
			 "Resuming primality test of M%%ld at iteration %%ld [%%.%df%%%%]\n",
			 PRECISION);
		sprintf (buf, fmt_mask, p, counter, pct);
		OutputStr (buf);
	}

/* Special hack for Emil Steen who got a roundoff error of 0.4915!!! */ 
/* This occured when continuing a v18 exponent that now uses a larger */
/* FFT size.  Since the error is greater than 0.48 he got hung in a loop. */
/* Since we have no idea if the roundoff error was 0.4915 or the deadly */
/* 0.5085, we will rotate the data in hopes of producing friendlier */
/* values in the FFT. */

	if (IniGetInt (INI_FILE, "ShiftHack", 0)) {
		gwadd (LLDATA, LLDATA);
		units_bit = (units_bit + 1) % p;
	}

/* Compute numbers in the lucas series, write out every 30 minutes to a file */

	iters = 0;
	while (counter < p) {
		int	stopping, echk;

/* Use this opportunity to perform other miscellaneous tasks that may */
/* be required by this particular OS port */

		doMiscTasks ();

/* Every so often, communicate with the server. */
/* Even more rarely, set flag to see if we have enough work queued up. */
/* Also, see if we should switch to factoring.  This is sometimes done */
/* when we get a new assignment that hasn't been factored.  To be sure */
/* that we will be able to immediately startup a LL test when this one */
/* completes, we will factor the next exponent and then come back to */
/* finish the LL test of this exponent. */

		if ((counter & 0x3F) == 0) {
			EXP_PERCENT_COMPLETE = (double) counter / (double) p;
			if ((counter & 0xFFFF) == 0) {
				ConditionallyUpdateEndDates ();
				CHECK_WORK_QUEUE = 1;
			}
			if (!communicateWithServer ()) escaped = 1;
			if ((counter & 0xFFFF) == 0) {
				if (getPriorityWork ()) priority_work = 1;
				updateRollingAverage (
					map_fftlen_to_timing (
						FFTLEN, GW_MERSENNE_MOD,
						CPU_TYPE, CPU_SPEED) *
					65536.0 * 24.0 / CPU_HOURS);
				startRollingAverage ();
			}
		}

/* Error check the last 50 iterations, before writing an */
/* intermediate file (either user-requested stop or a */
/* 30 minute interval expired), and every 128th iteration. */

		escaped |= stopCheck ();
		time (&current_time);
		saving = (current_time - start_time > write_time ||
			  counter == last_counter);
		stopping = escaped || priority_work;
		echk = stopping || saving || ERRCHK || (counter >= p - 50) ||
		       (counter & 127) == 0;
		MAXERR = 0.0;

/* Do a Lucas-Lehmer iteration */

		start_timer (0);
#ifndef SERVER_TESTING
		gwsetnormroutine (0, echk, 0);
		gwstartnextfft (!stopping && !saving && !(counter+1 == p) &&
				!(interimResidues && (counter+1) % interimResidues <= 2));
		lucas_fixup (&units_bit);
//gwnum qw = gwalloc ();
//gwfft (LLDATA, qw);
//gwfftfftmul (qw, qw, LLDATA);
//	gwcopy (LLDATA, qw);
//	gwmul (qw, LLDATA);
//gwfree(qw);
		gwsquare (LLDATA);
#endif
		end_timer (0);
		iters++;

/* If the sum of the output values is an error (such as infinity) */
/* then raise an error. */

		if (gw_test_illegal_sumout ()) {
			sprintf (buf, ERRMSG0, counter, p, ERRMSG1A);
			OutputBoth (buf);
			inc_error_count (2, &error_count);
			goto restart;
		}

/* Check that the sum of the input numbers squared is approximately */
/* equal to the sum of unfft results.  Since this check may not */
/* be perfect, check for identical results after a restart. */

		if (gw_test_mismatched_sums ()) {
			static double last_suminp = 0.0;
			static double last_sumout = 0.0;
			if (counter == last_counter &&
			    gwsuminp (LLDATA) == last_suminp &&
			    gwsumout (LLDATA) == last_sumout) {
				OutputBoth (ERROK);
				inc_error_count (3, &error_count);
				GWERROR = 0;
			} else {
				char	msg[100];
				sprintf (msg, ERRMSG1B,
					 gwsuminp (LLDATA),
					 gwsumout (LLDATA));
				sprintf (buf, ERRMSG0, counter, p, msg);
				OutputBoth (buf);
				last_counter = counter;
				last_suminp = gwsuminp (LLDATA);
				last_sumout = gwsumout (LLDATA);
				inc_error_count (0, &error_count);
				goto restart;
			}
		}

/* Check for excessive roundoff error  */

		if (echk && MAXERR > 0.40) {
			static double last_maxerr = 0.0;
			if (MAXERR < 0.48 &&
			    counter == last_counter &&
			    MAXERR == last_maxerr) {
				OutputBoth (ERROK);
				inc_error_count (3, &error_count);
				GWERROR = 0;
			} else {
				char	msg[100];
				sprintf (msg, ERRMSG1C, MAXERR);
				sprintf (buf, ERRMSG0, counter, p, msg);
				OutputBoth (buf);
				last_counter = counter;
				last_maxerr = MAXERR;
				inc_error_count (1, &error_count);
				goto restart;
			}
		}

/* Some special debugging code to dump out FFT results */

#ifdef XXX
xxx:{
unsigned long i;
char	buf[1024];
int	fd = _open (RESFILE, _O_TEXT | _O_RDWR | _O_CREAT, 0666);
while (_read (fd, buf, sizeof (buf)));
for (i = 0; i < FFTLEN; i++) {
	double *dblp = addr (LLDATA, i);
	if (*dblp < -0.0001 || *dblp > 0.0001) {
		sprintf (buf, "i: %ld, %10.6f\n", i, *dblp);
		_write (fd, buf, strlen (buf));
	}
}
_close (fd);
}
#endif

/* Update counter and maximum round-off error */

		counter++;
		if (ERRCHK) {
			if (MAXERR < reallyminerr && counter > 30)
				reallyminerr = MAXERR;
			if (MAXERR > reallymaxerr)
				reallymaxerr = MAXERR;
		}

/* Print a message every so often */

		if (counter % ITER_OUTPUT == 0) {
			char	fmt_mask[80];
			double	pct;
			pct = trunc_percent (counter * 100.0 / p);
			sprintf (fmt_mask, "%%.%df%%%% of %%ld", PRECISION);
			sprintf (buf, fmt_mask, pct, p);
			title (buf);
			sprintf (fmt_mask,
				 "Iteration: %%ld / %%ld [%%.%df%%%%]",
				 PRECISION);
			sprintf (buf, fmt_mask, counter, p, pct);
			OutputTimeStamp ();
			OutputStr (buf);
			if (ERRCHK && counter > 30) {
				OutputStr (".  Round off: ");
				sprintf (buf, "%10.10f", reallyminerr);
				OutputStr (buf);
				sprintf (buf, " to %10.10f", reallymaxerr);
				OutputStr (buf);
			}
			OutputStr (".  Per iteration time: ");
			divide_timer (0, iters);
			print_timer (0, TIMER_NL | TIMER_OPT_CLR);
			if (!CUMULATIVE_TIMING) iters = 0;
		}

/* Print a results file message every so often */

		if (counter % ITER_OUTPUT_RES == 0 || (NO_GUI && stopping)) {
			sprintf (buf, "Iteration %ld / %ld\n", counter, p);
			writeResults (buf);
		}

/* If an escape key was hit, write out the results and return */

		if (stopping) {
			if (! writeToFile (filename, counter,
					   units_bit, error_count)) {
				sprintf (buf, WRITEFILEERR, filename);
				OutputBoth (buf);
				return (FALSE);
			}
			if (escaped) return (FALSE);
			return (TRUE);
		}

/* Write results to a file every DISK_WRITE_TIME minutes */
/* On error, retry in 10 minutes (it could be a temporary */
/* disk-full situation) */

		if (saving) {
			write_time = DISK_WRITE_TIME * 60;
			if (! writeToFile (filename, counter,
					   units_bit, error_count)) {
				sprintf (buf, WRITEFILEERR, filename);
				OutputBoth (buf);
				if (write_time > 600) write_time = 600;
			}
			time (&start_time);
		}

/* Output the 64-bit residue at specified interims.  Also output the */
/* residues for the next two iterations so that we can compare our */
/* residues to programs that start counter at zero or one. */

		if (interimResidues && counter % interimResidues <= 2) {
			generateResidue64 (units_bit, &high32, &low32);
			sprintf (buf, 
				 "M%ld interim WX%d residue %08lX%08lX at iteration %ld\n",
				 p, PORT, high32, low32, counter);
			OutputBoth (buf);
		}

/* Write a save file every "interimFiles" iterations. */

		if (interimFiles && counter % interimFiles == 0) {
			char	interimfile[20];
			sprintf (interimfile, "%.8s.%03d",
				 filename, counter / interimFiles);
			writeToFile (interimfile, counter,
				     units_bit, error_count);
		}

/* This is a kludge to handle exponents smaller than the FFT length. */
/* The multiply normalize code can leave the second FFT value non-zero */
/* even though it should contain no data. */

//		if (p < FFTLEN) {
//			void gwtobinary (gwnum, giant);
//			void binarytogw (giant, gwnum);
//			giant g = popg(1);
//			gwtobinary (LLDATA, g);
//			binarytogw (g, LLDATA);
//			pushg (1);
//		}
	}

/* Check for a successful completion */
/* We found a prime if result is zero */
/* Note that all values of -1 is the same as zero */

	for (i = 0, isPrime = 1; isPrime && i < FFTLEN; i++) {
		long	val;
		get_fft_value (LLDATA, i, &val);
		if (isPrime <= 2 && val == 0)
			isPrime = 2;
		else if ((isPrime & 1) &&
			 (val == -1 ||
			  (p < FFTLEN && !is_big_word (i))))
			isPrime = 3;
		else
			isPrime = FALSE;
	}

/* Format the output message */

	if (isPrime) {
		sprintf (buf, "M%ld is prime! WX%d: %08lX\n",
			 p, PORT, SEC1 (p));
		high32 = low32 = 0;
	} else {
		generateResidue64 (units_bit, &high32, &low32);
		sprintf (buf,
			 "M%ld is not prime. Res64: %08lX%08lX. WX%d: %08lX,%ld,%08lX\n",
			 p, high32, low32, PORT,
			 SEC2 (p, high32, low32, units_bit, error_count),
			 units_bit, error_count);
	}

/* Output results to the screen, results file, and server */

	{
		struct primenetAssignmentResult pkt;
		memset (&pkt, 0, sizeof (pkt));
		pkt.exponent = p;
		pkt.resultType =
			isPrime ? PRIMENET_RESULT_PRIME : PRIMENET_RESULT_TEST;
		sprintf (pkt.resultInfo.residue,
			 "%08lX%08lX", high32, low32);
		spoolMessage (PRIMENET_ASSIGNMENT_RESULT, &pkt);
	}
	OutputStr (buf);
	rc = writeResults (buf);
	spoolMessage (PRIMENET_RESULT_MESSAGE, buf);

/* Delete the continuation files - assuming the results file write */
/* was successful. */

	if (rc) _unlink (filename);
	filename[0] = 'q';
	_unlink (filename);

/* Clear prime from the to-do list */

	updateWorkToDo (p, WORK_TEST, 0);

/* Clear rolling average start time in case Advanced/Test */
/* interrupted a Lucas-Lehmer test. */

	clearRollingStart ();
		
/* Output good news to the screen in an infinite loop */

	if (isPrime && !SILENT_VICTORY && !isKnownMersennePrime (p)) {
		title ("New Prime!!!");
		for ( ; ; ) {
			OutputStr ("New Mersenne Prime!!!!  ");
			OutputStr (buf);
			flashWindowAndBeep ();
			if (escapeCheck ()) return (FALSE);
			if (!communicateWithServer ()) return (FALSE);
		}
	}

/* All done */

	return (TRUE);

/* An error occured, sleep, then try restarting at last save point. */
/* Clear rolling average in case we will be backtracking past a */
/* 64K iteration boundary. */

restart:clearRollingStart ();

/* Output a message saying we are restarting */

	OutputBoth (ERRMSG2);
	OutputBoth (ERRMSG3);

/* Update the error count in the save file */

	writeNewErrorCount (filename, error_count);

/* Sleep five minutes before restarting */

	if (! SleepFive ()) return (FALSE);

/* Return so that last continuation file is read in */

	return (TRUE);
}

#define TORTURE1 "Beginning a continuous self-test to check your computer.\n"
#if defined (__linux__) || defined (__FreeBSD__) || defined (__EMX__)
#define TORTURE2 "Please read stress.txt.  Hit ^C to end this test.\n"
#else
#define TORTURE2 "Please read stress.txt.  Choose Test/Stop to end this test.\n"
#endif
#define SELFMSG1A "The program will now perform a self-test to make sure the\n"
#define SELFMSG1B "Lucas-Lehmer code is working properly on your computer.\n"
#define SELFMSG1C "This will take about an hour.\n"
#define SELFMSG5C "This will take about 29 hours.\n"
#define SELF1 "Test %i, %i Lucas-Lehmer iterations of M%ld using %ldK FFT length.\n"
#define SELFFAIL "FATAL ERROR: Final result was %08lX, expected: %08lX.\n"
#define SELFFAILW "FATAL ERROR: Writing to temp file.\n"
#define SELFFAILR "FATAL ERROR: Reading from temp file.\n"
char SELFFAIL1[] = "ERROR: ILLEGAL SUMOUT\n";
char SELFFAIL2[] = "FATAL ERROR: Resulting sum was %.16g, expected: %.16g\n";
char SELFFAIL3[] = "FATAL ERROR: Rounding was %.10g, expected less than 0.4\n";
char SELFFAIL4[] = "Possible hardware failure, consult readme.txt file, restarting test.\n";
char SELFFAIL5[] = "Hardware failure detected, consult stress.txt file.\n";

#define SELFPASS "Self-test %iK passed!\n"
char SelfTestIniMask[] = "SelfTest%iPassed";

#define MAX_SELF_TEST_ITERS	399

int SELF_TEST_ERRORS = 0;
int SELF_TEST_WARNINGS = 0;

struct self_test_info {
	unsigned long p;
	unsigned long reshi;
};

struct self_test_info SELF_TEST_DATA[MAX_SELF_TEST_ITERS] = {
{78643201, 0x2D9C8904}, {78643199, 0x7D469182}, {75497473, 0x052C7FD8},
{75497471, 0xCCE7495D}, {71303169, 0x467A9338}, {71303167, 0xBBF8B37D},
{68157441, 0xBE71E616}, {68157439, 0x93A71CC2}, {66060289, 0xF296BB99},
{66060287, 0x649EEF2A}, {62390273, 0xBC8DFC27}, {62390271, 0xDE7D5B5E},
{56623105, 0x0AEBF972}, {56623103, 0x1BA96297}, {53477377, 0x5455F347},
{53477375, 0xCE1C7F78}, {50331649, 0x3D746AC8}, {50331647, 0xE23F2DE6},
{49807361, 0xB43EF4C5}, {49807359, 0xA8BEB02D}, {47185921, 0xD862563C},
{47185919, 0x17281086}, {41943041, 0x0EDA1F92}, {41943039, 0xDE6911AE},
{39845889, 0x43D8A96A}, {39845887, 0x3D118E8F}, {37748737, 0x38261154},
{37748735, 0x22B34CD2}, {35651585, 0xB0E48D2E}, {35651583, 0xCC3340C6},
{34865153, 0xD2C00E6C}, {34865151, 0xFA644F69}, {33030145, 0x83E5738D},
{33030143, 0x6EDBC5B5}, {31195137, 0xFF9591CF}, {31195135, 0x04577C70},
{29884417, 0xACC36457}, {29884415, 0xC0FE7B1E}, {28311553, 0x780EB8F5},
{28311551, 0xE6D128C3}, {26738689, 0x09DC45B0}, {26738687, 0xDC7C074A},
{24903681, 0xA482CF1E}, {24903679, 0x4B3F5121}, {23592961, 0xAFE3C198},
{23592959, 0xCF9AD48C}, {20971521, 0x304EC13B}, {20971519, 0x9C4E157E},
{19922945, 0x83FE36D9}, {19922943, 0x9C60E7A2}, {18874369, 0x83A9F8CB},
{18874367, 0x5A6E22E0}, {17825793, 0xF3A90A5E}, {17825791, 0x6477CA76},
{17432577, 0xCAB36E6A}, {17432575, 0xB8F814C6}, {16515073, 0x91EFCB1C},
{16515071, 0xA0C35CD9}, {15597569, 0x12E057AD}, {15597567, 0xC4EFAEFD},
{14942209, 0x1C912A7B}, {14942207, 0xABA9EA6E}, {14155777, 0x4A943A4E},
{14155775, 0x00789FB9}, {13369345, 0x27A041EE}, {13369343, 0xA8B01A41},
{12451841, 0x4DC891F6}, {12451839, 0xA75BF824}, {11796481, 0xFDD67368},
{11796479, 0xE0237D19}, {10485761, 0x15419597}, {10485759, 0x154D473B},
{10223617, 0x26039EB7}, {10223615, 0xC9DFB1A4}, {9961473, 0x3EB29644},
{9961471, 0xE2AB9CB2}, {9437185, 0x42609D65}, {9437183, 0x77ED0792},
{8716289, 0xCCA0C17B}, {8716287, 0xD47E0E85}, {8257537, 0x80B5C05F},
{8257535, 0x278AE556}, {7798785, 0x55A2468D}, {7798783, 0xCF62032E},
{7471105, 0x0AE03D3A}, {7471103, 0xD8AB333B}, {7077889, 0xC516359D},
{7077887, 0xA23EA7B3}, {6684673, 0xA7576F00}, {6684671, 0x057E57F4},
{6422529, 0xC779D2C3}, {6422527, 0xA8263D37}, {6225921, 0xB46AEB2F},
{6225919, 0xD0A5FD5F}, {5898241, 0xE46E76F9}, {5898239, 0x29ED63B2},
{5505025, 0x83566CC3}, {5505023, 0x0B9CBE64}, {5242881, 0x3CC408F6},
{5242879, 0x0EA4D112}, {4980737, 0x6A2056EF}, {4980735, 0xE03CC669},
{4718593, 0x87622D6B}, {4718591, 0xF79922E2}, {4587521, 0xE189A38A},
{4587519, 0x930FF36C}, {4358145, 0xDFEBF850}, {4358143, 0xBB63D330},
{4128769, 0xC0844AD1}, {4128767, 0x25BDBFC3}, {3932161, 0x7A525A7E},
{3932159, 0xF30C9045}, {3735553, 0xFAD79E97}, {3735551, 0x005ED15A},
{3538945, 0xDDE5BA46}, {3538943, 0x15ED5982}, {3342337, 0x1A6E87E9},
{3342335, 0xECEEA390}, {3276801, 0x3341C77F}, {3276799, 0xACA2EE28},
{3112961, 0x2BDF9D2B}, {3112959, 0xA0AC8635}, {2949121, 0x36EDB768},
{2949119, 0x53FD5473}, {2785281, 0x66816C94}, {2785279, 0x059E8D6B},
{2654999, 0x07EE900D},{2621441, 0x2BC1DACD},{2621439, 0xBCBA58F1},
{2653987, 0xB005CACC},{2651879, 0x38DCD06B},{2654003, 0x1ED556E7},
{2620317, 0x09DB64F8},{2539613, 0x4146EECA},{2573917, 0x939DA3B3},
{2359297, 0x73A131F0},{2359295, 0x53A92203},{2646917, 0x71D4E5A2},
{2605473, 0xE11637FC},{2495213, 0x89D80370},{2540831, 0x2CF01FBB},
{2654557, 0x4106F46F},{2388831, 0xA508B5A7},{2654777, 0x9E744AA3},
{2584313, 0x800E9A61},{2408447, 0x8C91E8AA},{2408449, 0x437ECC01},
{2345677, 0x60AEE9C2},{2332451, 0xAB209667},{2330097, 0x3FB88055},
{2333851, 0xFE4ECF19},{2444819, 0x56BF33C5},{2555671, 0x9DC03527},
{2654333, 0xE81BCF40},{2543123, 0x379CA95D},{2432123, 0x5952676A},
{2321123, 0x24DCD25F},{2654227, 0xAC3B7F2B},
{2329999, 0xF5E902A5},{2293761, 0x9E4BBB8A},{2293759, 0x1901F07B},
{2236671, 0x45EB162A},{2193011, 0x382B6E4B},{2329001, 0x4FF052BB},
{2327763, 0x3B315213},{2325483, 0x0DC5165A},{2323869, 0xD220E27F},
{2315679, 0xF650BE33},{2004817, 0xC2FF3440},{2130357, 0xC25804D8},
{2288753, 0xA4DD9AAD},{2266413, 0x675257DB},{2244765, 0xC08FF487},
{2222517, 0x1A128B22},{2200339, 0x0EB0E827},{2328117, 0x0A24673A},
{2329557, 0x2E267692},{2188001, 0xD012AF6A},{2166567, 0x509BA41A},
{2144651, 0x54CFC0E6},{2122923, 0xA47068E6},{2100559, 0xACFAB4E1},
{2088461, 0xEA01E860},{2066543, 0x847DF0D0},{2044767, 0x04225888},
{2022823, 0x6EA34B32},{2328527, 0xC55E3E05},{2327441, 0x207C8CEC},
{2326991, 0x0A4F2ACD},{2009987, 0xE6A59DEF},
{1999999, 0xD645A18F},{1966081, 0xB88828A1},{1966079, 0x5BD87C45},
{1998973, 0xCBDD74F7},{1997651, 0x666B0CB1},{1675001, 0x50A94DB7},
{1977987, 0x30D1CD1F},{1955087, 0x5B9426A4},{1933071, 0x23C1AF0B},
{1911957, 0xF7699248},{1899247, 0x11C76E04},{1877431, 0xA3299B39},
{1855067, 0x35243683},{1833457, 0xCF630DC0},{1811987, 0x7C7022EC},
{1799789, 0xEFEC47B7},{1777773, 0x0F16E2D6},{1755321, 0x1AC5D492},
{1733333, 0x5DA0555E},{1711983, 0xDC19DA8B},{1699779, 0x2B44914E},
{1677323, 0x03D3980B},{1995091, 0x922E555B},{1993041, 0x0CA8451B},
{1991991, 0xDFFB212D},{1679779, 0x51D75E0F},{1684993, 0x048BBCE8},
{1970009, 0x646E0DFA},{1957445, 0xC8D244ED},{1999997, 0x5FC899D0},
{1998983, 0x1CD518AA},{1999007, 0xA9DD8591},
{1674999, 0xDB0169D8},{1638401, 0xD3F8A8C5},{1638399, 0xF270D8DD},
{1674997, 0xC824EF15},{1674551, 0xD844AEAD},{1674001, 0x8F5EFA50},
{1345001, 0x18EE2E2D},{1655083, 0x09B30DEE},{1633941, 0x0B87C8B1},
{1611557, 0x6B57E48D},{1599549, 0x48EA38B2},{1577771, 0xCE84D9DC},
{1555947, 0x6797EEF4},{1533349, 0xD6897409},{1511861, 0x8A8177AC},
{1499625, 0x56BB6FB3},{1477941, 0xF3DD8ED3},{1455931, 0x31A222C7},
{1433069, 0x28F01E1B},{1411747, 0x680C6E39},{1399449, 0xB7F01A54},
{1377247, 0xE656F652},{1355991, 0xB2AA2819},{1350061, 0x31F9A728},
{1673881, 0xA51D38E4},{1672771, 0x5474B6F9},{1671221, 0x2710DDEA},
{1670551, 0x31FC3838},{1660881, 0x4C5B22C5},{1650771, 0x998F747B},
{1655001, 0x164659A6},{1674339, 0xED2D23E2},
{1344999, 0x158AA064},{1310721, 0x5694A427},{1310719, 0x258BDDE3},
{1344997, 0x1D059D4F},{1344551, 0x60606AA3},{1344001, 0x9AC6AB36},
{1322851, 0x3A000D0A},{1300993, 0x77CB0184},
{1288771, 0x7431D9E2},{1266711, 0xB4BC4E8D},{1244881, 0x48BC9FF9},
{1222991, 0x3F5FC39E},{1200881, 0xD5DF4944},{1188441, 0xD9D8968B},
{1166661, 0xD4AB97F4},{1144221, 0x9940943B},{1122001, 0x647406B8},
{1100881, 0x3AD40CE0},{1088511, 0xD578BB51},{1066837, 0x2F82BFBB},
{1044811, 0x7C6EDDD1},{1022991, 0x6A1C2DD4},{1000001, 0x2879748F},
{1343881, 0xB59E8006},{1342771, 0x87563FFE},{1341221, 0x29AD6127},
{1340551, 0x17DB4ACB},{1330881, 0x9642F068},
{942079, 0xE528A9B0},{974849, 0x79791EDB},{983041, 0x29216C43},
{901121, 0x26C4E660},{917503, 0x5F244685},{933889, 0x62490F57},
{851967, 0x331AA906},{860161, 0x41185F27},{884735, 0x7BC7A661},
{802817, 0xA9645693},{819199, 0x48AFB0A5},{835585, 0x706437D3},
{753663, 0x99C43F31},{778241, 0x1729A6C4},{786431, 0x61080929},
{720897, 0x1E96863D},{737279, 0x1B07A764},{745473, 0x7BCE80AA},
{655359, 0x1107F161},{659457, 0x589C16A4},{688127, 0xD01E5A85},
{622593, 0x26F6FC8C},{630783, 0x4DD2E603},{638977, 0xC88F34B4},
{589823, 0x0290B60B},{602113, 0xEFCD5BA8},{614399, 0x6408F880},
{557057, 0xC30FE589},{565247, 0xF4CA3679},{573441, 0xF8F039AA},
{532479, 0x0072FE03},{540673, 0xDA0E0D99},{544767, 0x62443C6B},
{491521, 0x3F520DFA},{516095, 0xA6BD9423},{524289, 0xCD591388},
{466943, 0xE10EE929},{471041, 0x18752F40},{487423, 0x933FFF17},
{442369, 0xC22471C3},{450559, 0x025B1320},{458753, 0xE296CC00},
{417791, 0x080C803C},{425985, 0xB2095F04},{430079, 0x98B1EC61},
{393217, 0x26DD79ED},{401407, 0x2F0F75F9},{409601, 0xAEFAC2F8},
{372735, 0xCB6D00A2},{376833, 0x915D5458},{389119, 0x6188E38D},
{344065, 0x4D0C5089},{360447, 0x84AC5CFD},{368641, 0x72414364},
{319487, 0x24ED1BE9},{327681, 0x3101106A},{329727, 0x5BDB69AF},
{307201, 0x68536CD1},{311295, 0x69778074},{315393, 0x429D4950},
{286719, 0x1A31A686},{294913, 0xF55727C6},{301055, 0x33BDB242},
{272385, 0xEF6EC4B4},{278527, 0x05530FD5},{282625, 0x34A4E699},
{262143, 0xA9638844},{266241, 0xE0969CED},{270335, 0x14AD54BE},
{243713, 0xC19AEA91},{245759, 0x7538BF0B},{258049, 0x73F541AD},
{229375, 0x6E42B26A},{233473, 0x1964F897},{235519, 0x661BBC3F},
{215041, 0x04D5D2F0},{221183, 0xA89E7764},{225281, 0x20876BED},
{204799, 0xD20C2126},{208897, 0x9D4DCF0E},{212991, 0x1FF00E2A},
{194561, 0x6ED1CB70},{196607, 0x3190D5F5},{200705, 0xFAD28F5A},
{184319, 0x360EF08E},{186369, 0x0F001482},{188415, 0x86FCE4D6},
{164865, 0x4942B002},{172031, 0xC5AF29DB},{180225, 0x35D49D74},
{157695, 0x5422FACF},{159745, 0xB5CD03A1},{163839, 0x1CA6048E},
{150529, 0x7412F09C},{153599, 0xA9FAAE69},{155649, 0xA7B736AF},
{141311, 0x7A5D0730},{143361, 0x580F4DC4},{147455, 0x176B299A},
{135169, 0x65AC10A4},{136191, 0xC4591D37},{139265, 0xBCE1FC80},
{129023, 0xAFE1E7A8},{131073, 0xC5AAB12F},{133119, 0xDE51C35A},
{117761, 0x054A26F6},{121855, 0x55AF2385},{122881, 0x652827AC},
{112639, 0x6FA4DB24},{114689, 0x0BBAF161},{116735, 0xB85F0E8E},
{106497, 0xF833D925},{107519, 0x80F177D8},{110593, 0x1A56AA86},
{100351, 0x1DE12CE6},{102401, 0x19F967B4},{104447, 0xF9F3CDFD}
};

int selfTest (
	unsigned long pArg)
{
	unsigned long fftlen;
	int	i;

/* Set the process/thread priority */

	SetPriority ();

/* Clear rolling average start time in case we've */
/* interrupted a Lucas-Lehmer test. */

	clearRollingStart ();

/* Clear counters */

	SELF_TEST_ERRORS = 0;
	SELF_TEST_WARNINGS = 0;

/* Are we to self-test all fft lengths? (pArg = 0) */
/* Check for the torture-test case (pArg = 1) */

	if (pArg <= 1) {
		int	lengths[29] = {1024,8,10,896,768,12,14,640,512,16,20,448,384,24,28,320,256,32,40,224,192,48,56,160,128,64,80,112,96};
		int	min_fft, max_fft, test_time;
		time_t	start_time, current_time;

/* Make sure the user really wants to spend many hours doing this now */

		if (pArg == 0) {
			OutputStr (SELFMSG1A);
			OutputStr (SELFMSG1B);
			OutputStr (SELFMSG5C);
			test_time = 60;
		} else {
			OutputStr (TORTURE1);
			OutputStr (TORTURE2);
			test_time = IniGetInt (INI_FILE, "TortureTime", 15);
			time (&start_time);
		}

/* Now self-test each fft length */

		min_fft = IniGetInt (INI_FILE, "MinTortureFFT", 8);
		max_fft = IniGetInt (INI_FILE, "MaxTortureFFT", 4096);
		for ( ; ; ) {
		    for (i = 0; i < 29; i++) {
			if (lengths[i] < min_fft || lengths[i] > max_fft)
				continue;
			if (! selfTestInternal (lengths[i]*1024, test_time, i)) {
				if (pArg == 1) {
					char buf[100];
					int	hours, minutes;
					time (&current_time);
					minutes = (current_time - start_time) / 60;
					hours = minutes / 60;
					minutes = minutes % 60;
					OutputStr ("Torture Test ran ");
					if (hours) {
						sprintf (buf, "%d hours, ", hours);
						OutputStr (buf);
					}
					sprintf (buf, "%d minutes - %d errors, %d warnings.\n", minutes, SELF_TEST_ERRORS, SELF_TEST_WARNINGS);
					OutputStr (buf);
				}
				return (FALSE);
			}
		    }
		    if (pArg == 0) break;
		}
	}

/* Otherwise we are self-testing a specific word length. */

	else {
		char	iniName[32];

/* What fft length are we running? */

		fftlen = map_exponent_to_fftlen (pArg, GW_MERSENNE_MOD);

/* If fftlength is less than 64K return (we don't have any small exponents */
/* in our self test data) */

		if (fftlen < 65536) return (TRUE);

/* If so, make sure we haven't done so already. */

		sprintf (iniName, SelfTestIniMask, (int) (fftlen/1024));
		if (IniGetInt (LOCALINI_FILE, iniName, 0)) return (TRUE);

/* Make sure the user really wants to spend an hour doing this now */

		OutputStr (SELFMSG1A);
		OutputStr (SELFMSG1B);
		OutputStr (SELFMSG1C);

/* Do the self test */

		if (! selfTestInternal (fftlen, 60, -1)) return (FALSE);
	}

/* Self test completed! */

	return (TRUE);
}

int selfTestInternal (
	unsigned long fftlen,
	unsigned int test_time,	/* Number of minutes to self-test */
	int	torture_count)	/* Index into array of fft sizes */
{
	unsigned long k, limit;
	unsigned int i, iter, countdown;
	char	filename[16];
	char	buf[120];
	char	iniName[32];
	time_t	start_time, current_time;
static	int	data_index[29] = {0};

/* Set the title */

	title ("Self-Test");

/* Generate file name for temp files */

	strcpy (filename, "ptemp");
	strcat (filename, EXTENSION);

/* Determine the range from which we'll choose an exponent to test. */

	limit = map_fftlen_to_max_exponent (fftlen, GW_MERSENNE_MOD);

/* Get the current time */

	time (&start_time);

/* Start in the self test data array where we left off the last time */
/* torture test executed this FFT length. */

	i = (torture_count < 0) ? 0 : data_index[torture_count];

/* Loop testing various exponents from self test data array until */
/* time runs out */

	countdown = 1;
	for (iter = 1; ; iter++) {
		unsigned long p, reshi, reslo, units_bit;
		unsigned int ll_iters;

/* Find next self test data entry to work on */

		for ( ; ; i++) {

/* Wrap in the self test data array */

			if (i == MAX_SELF_TEST_ITERS) i = 0;

/* Now select the actual exponent */

			p = SELF_TEST_DATA[i].p;
			if (p > limit) continue;

/* The SSE2 carry propagation code gets into trouble if there are too */
/* few bits per FFT word!  Thus, we'll require at least 8 bits per */
/* word here. */

			if ((CPU_FLAGS & CPU_SSE2) && p / fftlen < 8) continue;

/* We've found an exponent to test! */

			break;
		}

/* Output start message */

		ll_iters = (p < 1000000) ? 1000 : 400;
		sprintf (buf, SELF1, iter, ll_iters, p, fftlen / 1024);
		OutputStr (buf);

/* Now run Lucas setup, for extra safety double the maximum allowable */
/* sum(inputs) vs. sum(outputs) difference. */

		lucasSetup (p, fftlen);
		MAXDIFF *= 2.0;

/* Init data area with a pre-determined value */

restart_test:	units_bit = 0;
		dbltogw (4.0, LLDATA);

/* Do Lucas-Lehmer iterations */

		for (k = 0; k < ll_iters; k++) {
			int	fd;
			short	type;
			unsigned long trash;

/* One Lucas-Lehmer test with error checking */

			gwsetnormroutine (0, 1, 0);
			gwstartnextfft (k != 100 && k != ll_iters - 1);
			lucas_fixup (&units_bit);
			gwsquare (LLDATA);

/* If the sum of the output values is an error (such as infinity) */
/* then raise an error. */

			if (gw_test_illegal_sumout ()) {
				OutputBoth (SELFFAIL1);
				OutputBoth (SELFFAIL4);
				SELF_TEST_WARNINGS++;
				goto restart_test;
			}

/* Check that the sum of the input numbers squared is approximately */
/* equal to the sum of unfft results. */

			if (gw_test_mismatched_sums ()) {
				sprintf (buf, SELFFAIL2,
					 gwsumout (LLDATA),
					 gwsuminp (LLDATA));
				OutputBoth (buf);
				OutputBoth (SELFFAIL5);
				SELF_TEST_ERRORS++;
				lucasDone ();
				return (FALSE);
			}

/* Make sure round off error is tolerable */

			if (MAXERR > 0.45) {
				sprintf (buf, SELFFAIL3, MAXERR);
				OutputBoth (buf);
				OutputBoth (SELFFAIL5);
				SELF_TEST_ERRORS++;
				lucasDone ();
				return (FALSE);
			}

/* Abort if user demands it */

			if (escapeCheck ()) {
				lucasDone ();
				return (FALSE);
			}

/* Test our ability to read and write files too. */

			if (k != 100) continue;
			if (--countdown) continue;
			if (! writeToFile (filename, 0, 0, 0)) {
				OutputBoth (SELFFAILW);
				SELF_TEST_ERRORS++;
				lucasDone ();
				return (FALSE);
			}
			countdown = 9;
			dbltogw (0.0, LLDATA); /* clear memory */
			if (! readFileHeader (filename, &fd, &type, &trash) ||
			    ! readFileData (fd, &trash, &trash)) {
				OutputBoth (SELFFAILR);
				SELF_TEST_ERRORS++;
				lucasDone ();
				return (FALSE);
			}
			_unlink (filename);
		}

/* Compare final 32 bits with the pre-computed array of correct residues */

		generateResidue64 (units_bit, &reshi, &reslo);
		lucasDone ();
		if (reshi != SELF_TEST_DATA[i].reshi) {
			sprintf (buf, SELFFAIL, reshi, SELF_TEST_DATA[i].reshi);
			OutputBoth (buf);
			OutputBoth (SELFFAIL5);
			SELF_TEST_ERRORS++;
			return (FALSE);
		}

/* Bump index into self test data array */

		i++;

/* Has time expired? */

		time (&current_time);
		if ((unsigned int) (current_time - start_time) >= test_time * 60) break;
	}

/* Save our position in self test data array for next time torture test */
/* executes this FFT length */

	if (torture_count >= 0) data_index[torture_count] = i;

/* We've passed the self-test.  Remember this in the .INI file */
/* so that we do not need to do this again. */

	sprintf (buf, SELFPASS, (int) (fftlen/1024));
	OutputBoth (buf);
	sprintf (iniName, SelfTestIniMask, (int) (fftlen/1024));
	IniWriteInt (LOCALINI_FILE, iniName, 1);
	return (TRUE);
}

/* Read a file of exponents to run LL iterations on as part of a QA process */
/* The format of this file is: */
/*	exponent,optional fft length,num iters,optional shift count,residue */

int lucas_QA (
	int	type)
{
	FILE	*fd;

/* Set the title */

	title ("QA");

/* Open QA file */

	fd = fopen ("qa", "r");
	if (fd == NULL) {
		OutputStr ("File named 'qa' could not be opened.\n");
		return (TRUE);
	}

/* Loop until the entire file is processed */

	for ( ; ; ) {
		unsigned long p, fftlen, iters, units_bit;
		char	buf[500], res[80];
		unsigned long reshi, reslo;
		unsigned long i, word, bit_in_word;
		double	maxsumdiff, maxerr, toterr, M, S;
/*		unsigned long bins[501]; */
		gwnum	t1, t2;
		unsigned int iters_unchecked;

/* Read a line from the file */

		p = 0;
		fscanf (fd, "%lu,%lu,%lu,%lu,%s\n",
			&p, &fftlen, &iters, &units_bit, &res);
		if (p == 0) break;

/* Now run Lucas setup */

		lucasSetup (p, fftlen);
		maxsumdiff = 0.0;
		maxerr = 0.0; toterr = 0.0;
/*		memset (bins, 0, sizeof (bins)); */
		iters_unchecked = (type >= 2) ? 12 : (type >= 1) ? 5 : 40;

/* Check for a randomized units bit */

		if (units_bit >= p) {
			unsigned long hi, lo;
			srand ((unsigned) time (NULL));
			units_bit = (rand () << 16) + rand ();
			if (CPU_FLAGS & CPU_RDTSC) { rdtsc (&hi,&lo); units_bit += lo; }
			units_bit = units_bit % p;
			sprintf (buf, "Units bit = %lu\n", units_bit);
			OutputBoth (buf);
		}

/* Init data area with a pre-determined value */

		bitaddr ((units_bit + 2) % p, &word, &bit_in_word);
		for (i = 0; i < FFTLEN; i++) {
			set_fft_value (LLDATA, i,
				       (type == 3) ? i % 97 :
				       (i == word) ? (1L << bit_in_word) : 0);
		}

/* The thorough, P-1, and ECM tests use more than one number */

		if (type >= 2) {
			t1 = gwalloc ();
			dbltogw (234872639921.0, t1);
			gwfft (t1, t1);
			t2 = gwalloc ();
			dbltogw (1982387192367.0, t2);
			gwfft (t2, t2);
			MAXDIFF *= 16;
		}

/* Do Lucas-Lehmer iterations */

		for (i = 0; i < iters; i++) {

/* One Lucas-Lehmer iteration with error checking */

			gwsetnormroutine (0, ERRCHK || (i & 127) == 64, 0);
			gwstartnextfft (i < iters / 2);
			if (type <= 1) {		/* Typical LL test */
				lucas_fixup (&units_bit);
				gwsquare (LLDATA);
			} else if (type == 2) {		/* Thorough test */
				unsigned long j;
				for (j = 0; j < (i & 7); j++) {
					gwadd (LLDATA, LLDATA);
					units_bit = (units_bit+1) % p;
				}
				if ((i & 15) == 13) {
					gwadd3quick (LLDATA, LLDATA, t1);
					gwsub3quick (t1, LLDATA, LLDATA);
					gwadd3 (LLDATA, LLDATA, t1);
					gwsub3 (t1, LLDATA, LLDATA);
					gwaddsub4 (LLDATA, LLDATA, t1, t2);
					gwaddsub (t1, LLDATA);
					gwadd (t2, LLDATA);
				}
				lucas_fixup (&units_bit);
				if ((i & 3) == 0) {
					gwsquare (LLDATA);
				} else if ((i & 3) == 1) {
					gwfft (LLDATA, LLDATA);
					gwfftfftmul (LLDATA, LLDATA, LLDATA);
				} else {
					gwfft (LLDATA, t1);
					gwfftmul (t1, LLDATA);
				}
			} else if (type == 3) {		/* Typical ECM run */
				lucas_fixup (&units_bit);
				gwfftsub3 (t1, t2, t2);
				gwfft (LLDATA, LLDATA);
				gwfftfftmul (t2, LLDATA, t2);
				gwswap (t1, LLDATA);
				gwswap (t2, LLDATA);
			}

/* Keep track of the standard deviation - see Knuth vol 2 */

			if ((type == 1 || type == 3) && i > iters_unchecked) {
/*				bins[(int)(MAXERR*1000.0+0.5)]++; */
				toterr += MAXERR;
				if (i == iters_unchecked + 1) {
					M = MAXERR;
					S = 0.0;
				} else {
					double	newM;
					newM = M + (MAXERR - M) /
						   (i - iters_unchecked);
					S = S + (MAXERR - M) * (MAXERR - newM);
					M = newM;
				}
			}

/* Maintain maximum error info */

			if (MAXERR > maxerr) maxerr = MAXERR;
			MAXERR = 0.0;

/* Maintain maximum suminp/sumout difference */

			if (fabs (gwsuminp (LLDATA) - gwsumout (LLDATA)) >
								maxsumdiff) {
				maxsumdiff = fabs (gwsuminp (LLDATA) -
						  gwsumout (LLDATA));
			}

/* If the sum of the output values is an error (such as infinity) */
/* then raise an error.  For some reason these bad values are treated */
/* as zero by the C compiler.  There is probably a better way to */
/* check for this error condition. */

			if (gw_test_illegal_sumout ()) {
				OutputBoth ("Warning: ILLEGAL SUMOUT\n");
				dbltogw (11.0, LLDATA);
				GWERROR = 0;
			}

/* Check that the sum of the input numbers squared is approximately */
/* equal to the sum of unfft results. */

			if (gw_test_mismatched_sums ()) {
				OutputBoth ("Warning: SUMOUT MISMATCH\n");
				GWERROR = 0;
			}

/* Abort if user demands it */

			if (escapeCheck ()) {
				lucasDone ();
				fclose (fd);
				return (FALSE);
			}

/* Test our ability to read and write files too. */

			if (type == 2 && (i & 1023) == 255) {
				char	filename[16];
				int	fd;
				short	type;
				unsigned long trash;
				strcpy (filename, "ptemp");
				strcat (filename, EXTENSION);
				if (! writeToFile (filename, 0, 0, 0))
					OutputBoth ("Warning: File write failed\n");
				dbltogw (0.0, LLDATA); /* clear memory */
				if (! readFileHeader (filename, &fd, &type, &trash) ||
				    ! readFileData (fd, &trash, &trash))
					OutputBoth ("Warning: File read failed\n");
				_unlink (filename);
			}
		}

/* Generate residue and cleanup */

		generateResidue64 (units_bit, &reshi, &reslo);
		lucasDone ();

/* Output array of distributions of MAXERR */

		if (type == 1 || type == 3) {
			S = sqrt (S / (iters - iters_unchecked - 1));
			toterr /= iters - iters_unchecked;
			sprintf (buf, "avg: %6.6f, stddev: %6.6f, #stdev to 0.5: %6.6f\n",
				 toterr, S, (0.50 - toterr) / S);
			OutputBoth (buf);
/*			for (i = 0; i < 501; i++) {
				sprintf (buf, "%lu,%lu\n", i, bins[i]);
				OutputBoth (buf);
			} */
		}

/* Compare residue with correct residue from the input file */

		sprintf (buf, "%08X%08X", reshi, reslo);
		if (stricmp (res, buf)) {
			sprintf (buf, "Warning: Residue mismatch. Expected %s\n", res);
			OutputBoth (buf);
		}

/* Output message */

		sprintf (buf, "Exp/iters: %lu/%lu, res: %08X%08X, maxerr: %6.6f, maxdiff: %9.9f/%9.9f\n",
			 p, iters, reshi, reslo, maxerr, maxsumdiff, MAXDIFF);
		OutputBoth (buf);
	}
	fclose (fd);

	return (TRUE);
}

/* Generate random FFT data for timing the Lucas-Lehmer code */

void generateRandomData (void)
{
	unsigned long i;

/* Fill data space with random values. */

	srand ((unsigned int) PARG);
	for (i = 0; i < FFTLEN; i++) {
		set_fft_value (LLDATA, i, rand() & 0xFF);
	}
}

/* Time a few iterations of an LL test on a given exponent */

void primeTime (
	unsigned long p,
	unsigned long iterations)
{
#define SAVED_LIMIT	10
	unsigned long i, saved, save_limit;
	char	buf[80];
	double	time, saved_times[SAVED_LIMIT];
	int	days, hours, minutes;

/* Set the process/thread priority */

	SetPriority ();

/* Clear all timers */

	clear_timers ();

/* Look for special values to run QA suite */

	if (p >= 9994 && p <= 9999) {
		lucas_QA (9999 - p);
		return;
	}

/* Init the FFT code */

	lucasSetup (p, 0);

/* Fill data space with random values. */

	generateRandomData ();

/* Do one squaring untimed, to prime the caches and start the */
/* post-FFT process going. */

	gwsetnormroutine (0, ERRCHK != 0, 0);
	gwstartnextfft (TRUE);
	gwsquare (LLDATA);

/* Compute numbers in the lucas series */
/* Note that for reasons unknown, we've seen cases where printing out
/* the times on each iteration greatly impacts P4 timings. */

	save_limit = (CPU_TYPE >= 12 && p <= 4000000) ? SAVED_LIMIT : 1;
	for (i = 0, saved = 0; i < iterations; i++) {
		start_timer (0);
		gwsquare (LLDATA);
		end_timer (0);
		timers[1] += timers[0];
		saved_times[saved++] = timers[0];
		timers[0] = 0;
		if (saved == save_limit || i == iterations - 1) {
			unsigned long j;
			for (j = 0; j < saved; j++) {
				OutputStr ("p: ");
				OutputNum (p);
				OutputStr (".  Time: ");
				timers[0] = saved_times[j];
				print_timer (0, TIMER_MS | TIMER_NL | TIMER_CLR);
			}
			saved = 0;
		}
		if (escapeCheck ()) {
			lucasDone ();
			return;
		}
	}
	lucasDone ();
	time = timer_value (1);

/* Print an estimate for how long it would take to test this number */

	OutputStr ("Iterations: ");
	OutputNum (iterations);
	OutputStr (".  Total time: ");
	print_timer (1, TIMER_NL | TIMER_CLR);
	time = time * p / iterations;
	days = (int) (time / 86400.0); time -= (double) days * 86400.0;
	hours = (int) (time / 3600.0); time -= (double) hours * 3600.0;
	minutes = (int) (time / 60.0);
	OutputStr ("Estimated time to complete this exponent: ");
	sprintf (buf, days == 1 ? "%d day, " : "%d days, ", days);
	OutputStr (buf);
	sprintf (buf, hours == 1 ? "%d hour, " : "%d hours, ", hours);
	OutputStr (buf);
	sprintf (buf, minutes == 1 ? "%d minute, " : "%d minutes.\n", minutes);
	OutputStr (buf);
//unsigned long *z;
//z = (unsigned long *) COPYZERO[0];
//if (z != NULL) for (i = 0; i < 32; i++) {
//sprintf (buf, "timer %d: %d\n", i, z[i]);
//if(z[i]) OutputBoth (buf);}
}

#define BENCH1 "\nYour timings will be written to the results.txt file.\n"
#define BENCH2 "Compare your results to other computers at http:\\www.mersenne.org\\bench.htm\n"
#define BENCH3 "That web page also contains instructions on how your results can be included.\n\n"

/* Time a few iterations of many FFT lengths */

void primeBench (void)
{
	unsigned long i, j, iterations;
	double	best_time;
	char	buf[80];
	int	fft_lengths[12] = {256, 320, 384, 448, 512, 640, 768, 892, 1024, 1280, 1536, 1792};

/* Set the process/thread priority */

	SetPriority ();

/* Output startup message */

	OutputStr (BENCH1);
	OutputBoth (BENCH2);
	OutputBoth (BENCH3);

/* Loop over all 12 FFT lengths */

	for (i = 0; i < 12; i++) {

/* Initialize for this FFT length.  Compute the number of iterations to time. */
/* This is based on the fact that it doesn't take too long for my 1400 MHz P4 */
/* to run 10 iterations of a 1792K FFT. */

		lucasSetup (5000000, fft_lengths[i] * 1024);
		iterations = 10 * 1792 * CPU_SPEED / 1400 / fft_lengths[i];
		if (iterations < 10) iterations = 10;

/* Output start message for this FFT length */

		sprintf (buf, "Timing %d iterations at %dK FFT length.  ", iterations, fft_lengths[i]);
		OutputStr (buf);

/* Fill data space with random values. */

		generateRandomData ();

/* Do one squaring untimed, to prime the caches and start the */
/* POSTFFT optimization going. */

		gwsetnormroutine (0, 0, 0);
		gwstartnextfft (TRUE);
		gwsquare (LLDATA);

/* Compute numbers in the lucas series */
/* Note that for reasons unknown, we've seen cases where printing out
/* the times on each iteration greatly impacts P4 timings. */

		for (j = 0; j < iterations; j++) {
			if (escapeCheck ()) {
				lucasDone ();
				return;
			}
			clear_timers ();
			start_timer (0);
			gwsquare (LLDATA);
			end_timer (0);
			if (j == 0 || timers[0] < best_time) best_time = timers[0];
		}
		lucasDone ();

/* Print the best time for this FFT length */

		OutputStr ("Best time: ");
		sprintf (buf, "Best time for %dK FFT length: ", fft_lengths[i]);
		writeResults (buf);
		timers[0] = best_time;
		print_timer (0, TIMER_NL | TIMER_MS | TIMER_OUT_BOTH);
	}
	OutputStr ("\n");
}

/* Wrapper code that verifies any factors found by the assembly code */

int factorAndVerify (
	unsigned long p)
{
	unsigned long hsw, msw;
	int	res;

/* Remember starting point in case of an error */

	hsw = FACHSW;
	msw = FACMSW;

/* Call assembly code */

loop:	res = factor64 ();

/* If a factor was not found, return. */

	if (res == 2) return (2);

/* Otherwise verify the factor. */

	if (FACHSW || FACMSW || FACLSW > 1) {
		giant	f, x;

		f = newgiant (100);
		itog ((int) FACHSW, f);
		gshiftleft (32, f);
		uladdg (FACMSW, f);
		gshiftleft (32, f);
		uladdg (FACLSW, f);

		x = newgiant (100);
		itog (2, x);
		powermod (x, p, f);
		res = isone (x);
	
		free (f);
		free (x);

		if (res) return (1);
	}

/* If factor is no good, print an error message, sleep, and */
/* restart the factoring code. */

	OutputBoth ("ERROR: Incorrect factor found.\n");
	if (! SleepFive ()) return (FALSE);
	FACLSW = p;
	SRCARG = FACDATA;
	setupf ();
	FACHSW = hsw;
	FACMSW = msw;
	goto loop;
}

/* Compute percent completion of a factoring job */

double facpct (
	short	pass,
	unsigned int bits,
	unsigned long endpthi,
	unsigned long endptlo)
{
	double	startpt, endpt, current;

	current = FACHSW * 4294967296.0 + FACMSW;
	endpt = endpthi * 4294967296.0 + endptlo;
	if (current > endpt) current = endpt;
        if (bits < 32) bits = 32;
        startpt = pow (2.0, bits-32);
	return ((pass + (current - startpt) / (endpt - startpt)) / 16.0 * 100.0);
}

/* Trial factor a Mersenne number prior to running a Lucas-Lehmer test */

char FACMSG[] = "Factoring M%%ld to 2^%%d is %%.%df%%%% complete.";

int primeFactor (
	unsigned long p,		/* Exponent to factor */
	unsigned int bits,		/* How far already factored in bits */
	int	*result,		/* Returns true if factor found */
	int	work_type,		/* Work type from worktodo.ini file */
	int	fd)			/* Continuation file handle or zero */
{
	int	continuation, old_style;
	long	write_time = DISK_WRITE_TIME * 60;
	unsigned int test_bits;
	unsigned long endpthi, endptlo;
	short	pass;
	time_t	start_time, current_time;
	struct work_unit w;
	char	filename[20];
	char	buf[80], str[80];

/* Clear all timers */

	clear_timers ();

/* Get the current time */

	time (&start_time);

/* Init temporary file name */

	tempFileName (filename, p);

/* Create a work packet for later updating of the rolling average */
 
	w.work_type = WORK_FACTOR;
	w.p = p;
	w.bits = bits;

/* Determine how much we should factor (in bits) */

	test_bits = factorLimit (p, work_type == WORK_FACTOR);
	if (test_bits < 32) test_bits = 32;

/* By default we do not do the old style of factoring which was 16 passes */
/* each pass testing from bits to test_bits */

	old_style = FALSE;

/* Is this a continuation?  If so, read continuation file. */
/* There are two types of continuation files, handle both */

	if (fd) {
		short	shortdummy;
		if (_read (fd, &pass, sizeof (short)) != sizeof (short))
			goto readerr;
		pass--;
		if (_read (fd, &shortdummy, sizeof (short)) != sizeof (short))
			goto readerr;
		*result = shortdummy;
		if (pass < 900) {
			unsigned long startpt;
			FACHSW = 0;
			if (_read (fd, &FACMSW, sizeof (long)) != sizeof (long))
				goto readerr;
			if (_read (fd, &startpt, sizeof (long)) != sizeof (long))
				goto readerr;
			endpthi = 0;
			if (_read (fd, &endptlo, sizeof (long)) != sizeof (long))
				goto readerr;
			old_style = TRUE;
		} else {
			if (_read (fd, &shortdummy, sizeof (short)) != sizeof (short))
				goto readerr;
			old_style = shortdummy;
			if (_read (fd, &shortdummy, sizeof (short)) != sizeof (short))
				goto readerr;
			bits = shortdummy;
			if (_read (fd, &pass, sizeof (short)) != sizeof (short))
				goto readerr;
			if (_read (fd, &FACHSW, sizeof (long)) != sizeof (long))
				goto readerr;
			if (_read (fd, &FACMSW, sizeof (long)) != sizeof (long))
				goto readerr;
			if (_read (fd, &endpthi, sizeof (long)) != sizeof (long))
				goto readerr;
			if (_read (fd, &endptlo, sizeof (long)) != sizeof (long))
				goto readerr;
		}
		_close (fd);
		continuation = TRUE;
	} else {
		startRollingAverage ();
		*result = FALSE;
		continuation = FALSE;
	}

/* Is prime already factored enough?  If so, return.  However, if we are */
/* processing a continuation file, then return the result from the file. */
/* This bizarre case happens when someone lowers the FactorOverride value */
/* in the middle of a factoring job. */

	if (bits >= test_bits) {
		if (fd) goto done;
		return (TRUE);
	}

/* Init the title */

	sprintf (buf, "Factoring M%ld", p);
	title (buf);
	sprintf (buf, "%s factoring M%ld to 2^%d\n",
		 fd ? "Resuming" : "Starting", p, test_bits);
	OutputStr (buf);

/* Load the factoring code */

	factorSetup ();

/* Loop testing larger and larger factors until we've tested to the */
/* appropriate number of bits.  Advance one bit at a time to minimize wasted */
/* time looking for a second factor after a first factor is found. */

	while (test_bits > bits) {
	    unsigned int end_bits;
	    unsigned long iters, iters_r;
	    int	stopping = 0;

/* Advance one bit at a time to minimize wasted time looking for a */
/* second factor after a first factor is found. */

	    end_bits = (bits < 50) ? 50 : bits + 1;
	    if (old_style) {
		    unsigned long x;
		    if (endpthi || endptlo == 0xFFFFFFFF) {
			    end_bits = 64;
			    endpthi = 1;
			    endptlo = 0;
		    }
		    else for (x = endptlo, end_bits = 31; x; x >>= 1) end_bits++;
	    }
	    if (end_bits > test_bits) end_bits = test_bits;

/* Compute the ending point for each pass */

	    if (!continuation) {
		if (end_bits < 64) {
			endpthi = 0;
			endptlo = 1L << (end_bits-32);
		} else {
			endpthi = 1L << (end_bits-64);
			endptlo = 0;
		}
	    }

/* Sixteen passes.  Two for the 1 or 7 mod 8 factors times two for the */
/* 1 or 2 mod 3 factors times four for the 1, 2, 3, or 4 mod 5 factors. */

	    iters_r = 0;
	    iters = 0;
	    if (! continuation) pass = 0;
	    for ( ; pass < 16; pass++) {

/* Set the starting point only if we are not resuming from */
/* a continuation file.  For no particularly good reason we */
/* quickly redo trial factoring for factors below 2^50. */

		if (continuation)
			continuation = FALSE;
		else {
			if (bits < 50) {
				FACHSW = 0;
				FACMSW = 0;
			} else if (bits < 64) {
				FACHSW = 0;
				FACMSW = 1L << (bits-32);
			} else {
				FACHSW = 1L << (bits-64);
				FACMSW = 0;
			}
		}

/* Only test for factors less than 2^32 on the first pass */

		if (FACHSW == 0 && FACMSW == 0 && pass != 0) FACMSW = 1;

/* Setup the factoring program */

		FACPASS = pass;
		FACLSW = p;
		SRCARG = FACDATA;
		setupf ();

/* Loop until all factors tested or factor found */

		for ( ; ; ) {
			int	res;

/* Use this opportunity to perform other miscellaneous tasks that may */
/* be required by this particular OS port */

			doMiscTasks ();

/* Periodically set global vars for Test/Status and CommunicateWithServer */

			if ((iters & 0xFFF) == 0) {
				EXP_BEING_WORKED_ON = p;
				EXP_BEING_FACTORED = 1;
				if (*result)
					EXP_PERCENT_COMPLETE = 999.0;
				else
					EXP_PERCENT_COMPLETE =
						pow (0.5, test_bits - end_bits + 1) *
						(1.0 + facpct (pass, bits, endpthi, endptlo) / 100.0);
			}

/* Do a chunk of factoring */

#ifdef SERVER_TESTING
			if (FACMSW & 0x80000000) FACHSW++;
			FACMSW=FACMSW*2+1000;
			if ((rand() & 0x1FF) == 0x111) break;
#else
			start_timer (0);
			res = factorAndVerify (p);
			end_timer (0);
			if (res != 2) break;
#endif

/* Send queued messages to the server every so often */

			if ((iters_r & 0x7F) == 0) {
				if (!communicateWithServer ()) stopping = 1;
			}

/* Set flag if we are stopping */

			stopping |= stopCheck ();

/* Output informative message */

			if (++iters >= ITER_OUTPUT) {
				char	fmt_mask[80];
				double	percent;
				percent = facpct (pass, bits, endpthi, endptlo);
				percent = trunc_percent (percent);
				sprintf (fmt_mask, FACMSG, PRECISION);
				sprintf (buf, fmt_mask, p, end_bits, percent);
				OutputTimeStamp ();
				OutputStr (buf);
				OutputStr ("  Time: ");
				print_timer (0, TIMER_NL | TIMER_OPT_CLR);
				sprintf (fmt_mask, "%%.%df%%%%", PRECISION);
				sprintf (buf, fmt_mask, percent);
				title (buf);
				iters = 0;
			}

/* Output informative message */

			if (++iters_r >= ITER_OUTPUT_RES ||
			    (NO_GUI && stopping)) {
				char	fmt_mask[80];
				double	percent;
				percent = facpct (pass, bits, endpthi, endptlo);
				percent = trunc_percent (percent);
				sprintf (fmt_mask, FACMSG, PRECISION);
				sprintf (buf, fmt_mask, p, end_bits, percent);
				strcat (buf, "\n");
				writeResults (buf);
				iters_r = 0;
			}

/* Test for completion */

			if (FACHSW > endpthi ||
			    (FACHSW == endpthi && FACMSW >= endptlo))
				goto nextpass;

/* If an escape key was hit, write out the results and return */

			time (&current_time);
			if (stopping || current_time-start_time > write_time) {
				short	shortdummy;
				long	longdummy;
				fd = _creat (filename, 0666);
				_close (fd);
				fd = _open (filename, _O_BINARY | _O_RDWR);
				shortdummy = 2;
				_write (fd, &shortdummy, sizeof (short));
				longdummy = 0;
				_write (fd, &longdummy, sizeof (long));
				shortdummy = 999;
				_write (fd, &shortdummy, sizeof (short));
				shortdummy = *result;
				_write (fd, &shortdummy, sizeof (short));
				shortdummy = old_style;
				_write (fd, &shortdummy, sizeof (short));
				shortdummy = bits;
				_write (fd, &shortdummy, sizeof (short));
				_write (fd, &pass, sizeof (short));
				_write (fd, &FACHSW, sizeof (long));
				_write (fd, &FACMSW, sizeof (long));
				_write (fd, &endpthi, sizeof (long));
				_write (fd, &endptlo, sizeof (long));
				_commit (fd);
				_close (fd);
				if (stopping) {
					factorDone ();
					return (FALSE);
				}
				start_time = current_time;
			}
		}

/* Format the output message */

		makestr (FACHSW, FACMSW, FACLSW, str);

/* Output results to the screen, results file, and server */

		{
			struct primenetAssignmentResult pkt;
			memset (&pkt, 0, sizeof (pkt));
			pkt.exponent = p;
			pkt.resultType = PRIMENET_RESULT_FACTOR;
			strcpy (pkt.resultInfo.factor, str);
			spoolMessage (PRIMENET_ASSIGNMENT_RESULT, &pkt);
		}
		sprintf (buf, "M%ld has a factor: %s\n", p, str);
		OutputBoth (buf);
		spoolMessage (PRIMENET_RESULT_MESSAGE, buf);

/* Set flag indicating a factor has been found */

		*result = TRUE;

/* We used to continue factoring to find a smaller factor in a later pass. */
/* However, there was a bug - restarting from the save file skipped the */
/* further factoring AND the time it takes to search for smaller factors */
/* is getting longer and longer as we factor deeper and deeper.  Therefore, */
/* in version 20 I've elected to no longer search for smaller factors. */
/* The one exception is users that are using FactorOverride to locate */
/* small factors. */

		if (work_type != WORK_FACTOR ||
		    IniGetInt (INI_FILE, "FactorOverride", 0) <= 60) break;

		if (FACMSW != 0xFFFFFFFF) {
			endpthi = FACHSW;
			endptlo = FACMSW+1;
		} else {
			endpthi = FACHSW+1;
			endptlo = 0;
		}

/* Every so often, set flag to see if we have enough work queued up */
/* Also, update the percent complete */

nextpass:	ConditionallyUpdateEndDates ();
		CHECK_WORK_QUEUE = 1;
	    }

/* If we've found a factor, then we are done factoring */

	    if (*result) break;

/* Advance the how far factored variable */

	    bits = end_bits;
	    old_style = FALSE;
	}

/* Clear prime from the to-do list */

done:	if (*result)
		updateWorkToDo (p, WORK_FACTOR, 0);

/* Output message if no factor found */

	else {
		struct primenetAssignmentResult pkt;

		memset (&pkt, 0, sizeof (pkt));
		pkt.exponent = p;
		pkt.resultType = PRIMENET_RESULT_NOFACTOR;
		pkt.resultInfo.how_far_factored = bits;
		spoolMessage (PRIMENET_ASSIGNMENT_RESULT, &pkt);

		sprintf (buf, "M%ld no factor to 2^%d, WX%d: %08lX\n",
			 p, bits, PORT, SEC3 (p));
		OutputBoth (buf);
		spoolMessage (PRIMENET_RESULT_MESSAGE, buf);

		updateWorkToDo (p, WORK_FACTOR, bits);

/* Update the rolling average when no factor is found */

		updateRollingAverage (raw_work_estimate (&w) * 24.0 / CPU_HOURS);
	}

/* Delete the continuation file */

	_unlink (filename);

/* Clear rolling average start time in case we've */
/* interrupted a Lucas-Lehmer test. */

	clearRollingStart ();

/* All done */

	factorDone ();
	return (TRUE);

/* Return a kludgy error code to indicate an error reading intermediate file */

readerr:_close (fd);
	*result = 999;
	return (TRUE);
}


/* Factor a range of primes using factors of the specified size */

char NOFAC[] = "M%ld no factor from 2^%d to 2^%d, WX%d: %08lX\n";

int primeSieve (
	unsigned long startp,
	unsigned long endp,
	unsigned short minbits,
	unsigned short maxbits,
	int	fd)			/* Continuation file */
{
	long	write_time = DISK_WRITE_TIME * 60;
	char	filename[20], buf[80], str[80];
	int	increasing, continuation;
	unsigned long p, endpthi, endptlo, iters, iters_r;
	short	pass, found_factor;
	int	stopping = 0;
	time_t	start_time, current_time;

/* Clear all timers */

	clear_timers ();

/* Get the current time */

	time (&start_time);

/* Clear rolling average start time in case we've */
/* interrupted a Lucas-Lehmer test. */

	clearRollingStart ();

/* Special value indicates check all factors in a file named "factors" */

	if (startp == 8888) {
		primeSieveTest ();
		return (FALSE);
	}

/* Is this a continuation?  If so, read continuation file. */

	if (fd) {
		short	oldfromfile;
		_read (fd, &startp, sizeof (long));
		_read (fd, &endp, sizeof (long));
		_read (fd, &minbits, sizeof (short));
		_read (fd, &maxbits, sizeof (short));
		_read (fd, &oldfromfile, sizeof (short));
		_read (fd, &p, sizeof (long));
		_read (fd, &FACMSW, sizeof (long));
		_read (fd, &pass, sizeof (short));
		_read (fd, &endptlo, sizeof (long));
		_read (fd, &found_factor, sizeof (short));
		FACHSW = 0;
		_read (fd, &FACHSW, sizeof (long));
		endpthi = 0;
		_read (fd, &endpthi, sizeof (long));
		_close (fd);
		continuation = TRUE;
	} else
		continuation = FALSE;

/* Init filename */

	tempFileName (filename, 0);

/* Load the factoring code */

	factorSetup ();

/* Loop until all the entire range is factored */

	increasing = (startp <= endp);
	if (!continuation) p = startp;
	if ((p & 1) == 0) p = increasing ? p + 1 : p - 1;
	for ( ; ; p = increasing ? p + 2 : p - 2) {

/* Reached endp? */

		if (increasing) {
			if (p > endp) break;
		} else {
			if (p < endp) break;
		}

/* Is p a prime? */

		if (!isPrime (p)) goto nextp;

/* Determine how much we should factor */

		if (!continuation) {
			if (maxbits < 64) {
				endpthi = 0;
				endptlo = 1L << (maxbits-32);
			} else {
				endpthi = 1L << (maxbits-64);
				endptlo = 0;
			}
		}

/* Sixteen passes! two for the 1 or 7 mod 8 factors times two for the */
/* 1 or 2 mod 3 factors times four for the 1, 2, 3, or 4 mod 5 factors. */

		if (!continuation) {
			pass = 1;
			found_factor = FALSE;
		}
		for ( ; pass <= 16; pass++) {

/* Setup the factoring program */

		if (!continuation) {
			if (minbits <= 32) {
				FACHSW = 0;
				FACMSW = 0;
			} else if (minbits <= 64) {
				FACHSW = 0;
				FACMSW = 1L << (minbits-33);
			} else {
				FACHSW = 1L << (minbits-65);
				FACMSW = 0;
			}
			if (pass > 1 && FACMSW == 0) FACMSW = 1;
		}
		FACPASS = pass-1;
		FACLSW = p;
		SRCARG = FACDATA;
		setupf ();

/* Loop until all factors tested or factor found */

		iters = 0;
		iters_r = 0;
		for ( ; ; ) {
			int	res;

/* Test for completion */

			if (FACHSW > endpthi ||
			    (FACHSW == endpthi && FACMSW >= endptlo))
				goto nextpass;

/* Factor some more */

			start_timer (0);
			res = factorAndVerify (p);
			end_timer (0);
			if (res != 2) break;

/* Send queued messages to the server every so often */

			if ((iters & 0x7F) == 0) {
				if (!communicateWithServer ()) stopping = 1;
			}

/* Set flag if we are stopping */

			stopping |= stopCheck ();

/* Output informative message */

			if (++iters >= ITER_OUTPUT) {
				char	fmt_mask[80];
				double	percent;
				percent = facpct ((short) (pass-1), minbits-1, endpthi, endptlo);
				percent = trunc_percent (percent);
				sprintf (fmt_mask, FACMSG, PRECISION);
				sprintf (buf, fmt_mask, p, maxbits, percent);
				OutputTimeStamp ();
				OutputStr (buf);
				OutputStr ("  Time: ");
				print_timer (0, TIMER_NL | TIMER_OPT_CLR);
				sprintf (fmt_mask, "%%.%df%%%%", PRECISION);
				sprintf (buf, fmt_mask, percent);
				title (buf);
				iters = 0;
			}

/* Output informative message */

			if (++iters_r >= ITER_OUTPUT_RES) {
				char	fmt_mask[80];
				double	percent;
				percent = facpct ((short) (pass-1), minbits-1, endpthi, endptlo);
				percent = trunc_percent (percent);
				sprintf (fmt_mask, FACMSG, PRECISION);
				sprintf (buf, fmt_mask, p, maxbits, percent);
				strcat (buf, "\n");
				writeResults (buf);
				iters_r = 0;
			}

/* If an escape key was hit, write out the results and return */

			time (&current_time);
			if (stopping || current_time-start_time > write_time) {
				short	three = 3;
				fd = _open (filename, _O_BINARY | _O_WRONLY | _O_TRUNC | _O_CREAT, 0666);
				_write (fd, &three, sizeof (short));
				_write (fd, &p, sizeof (long)); /* dummy */
				_write (fd, &startp, sizeof (long));
				_write (fd, &endp, sizeof (long));
				_write (fd, &minbits, sizeof (short));
				_write (fd, &maxbits, sizeof (short));
				_write (fd, &three, sizeof (short));
				_write (fd, &p, sizeof (long));
				_write (fd, &FACMSW, sizeof (long));
				_write (fd, &pass, sizeof (short));
				_write (fd, &endptlo, sizeof (long));
				_write (fd, &found_factor, sizeof (short));
				_write (fd, &FACHSW, sizeof (long));
				_write (fd, &endpthi, sizeof (long));
				_commit (fd);
				_close (fd);
				if (stopping) {
					factorDone ();
					return (FALSE);
				}
				start_time = current_time;
			}
		}

/* Format the output message */

		makestr (FACHSW, FACMSW, FACLSW, str);

/* Output results */

		{
			struct primenetAssignmentResult pkt;
			memset (&pkt, 0, sizeof (pkt));
			pkt.exponent = p;
			pkt.resultType = PRIMENET_RESULT_FACTOR;
			strcpy (pkt.resultInfo.factor, str);
			spoolMessage (PRIMENET_ASSIGNMENT_RESULT, &pkt);
		}
		sprintf (buf, "M%ld has a factor: %s\n", p, str);
		OutputBoth (buf);
		spoolMessage (PRIMENET_RESULT_MESSAGE, buf);

/* Next pass */

		if (FACMSW != 0xFFFFFFFF) {
			endpthi = FACHSW;
			endptlo = FACMSW+1;
		} else {
			endpthi = FACHSW+1;
			endptlo = 0;
		}
		found_factor = TRUE;
nextpass:	continuation = FALSE;
		}

/* Output message if no factor found */

		if (! found_factor) {
			sprintf (buf, NOFAC, p, minbits-1, maxbits, PORT, SEC4 (p));
			OutputBoth (buf);
			spoolMessage (PRIMENET_RESULT_MESSAGE, buf);
		}

/* Factor next prime */

nextp:		continuation = FALSE;
	}

/* Delete the continuation file */

	_unlink (filename);

/* Delete the entry from worktodo.ini */

	updateWorkToDo (0, WORK_FACTOR, 0);

/* All done */

	factorDone ();
	return (TRUE);
}


/* Test the factoring program */

void primeSieveTest (void)
{
	char	buf[500];
	FILE	*fd;
	unsigned long p;

/* Load the factoring code */

	factorSetup ();

/* Open factors file */

	fd = fopen ("factors", "r");

/* Loop until all the entire range is factored */

	while (fscanf (fd, "%ld", &p) && p) {
		unsigned long fachi, facmid, faclo;
		unsigned long i;
		char fac[480];
		char *f;

/* What is the factor? */

		fscanf (fd, "%s", fac);
		fachi = facmid = faclo = 0;
		for (f = fac; *f; f++) {
			if (*f < '0' || *f > '9') continue;
			RES = *f - '0';
			CARRYL = 0;
			muladdhlp (faclo, 10);
			faclo = RES;
			RES = CARRYL;
			CARRYL = 0;
			muladdhlp (facmid, 10);
			facmid = RES;
			fachi = fachi * 10 + CARRYL;
			if (fachi >= 4096) {
				sprintf (buf, "%s factor too big.\n", fac);
				OutputBoth (buf);
				goto nextp;
			}
		}

/* See if p is a prime */

		if (! isPrime (p)) {
			sprintf (buf, "%ld not a prime.\n", p);
			OutputBoth (buf);
			goto nextp;
		}

/* Setup the factoring program */

		i = (fachi % 120 * 16 + facmid % 120 * 16 + faclo % 120) % 120;
		if (i == 1) FACPASS = 0;
		else if (i == 7) FACPASS = 1;
		else if (i == 17) FACPASS = 2;
		else if (i == 23) FACPASS = 3;
		else if (i == 31) FACPASS = 4;
		else if (i == 41) FACPASS = 5;
		else if (i == 47) FACPASS = 6;
		else if (i == 49) FACPASS = 7;
		else if (i == 71) FACPASS = 8;
		else if (i == 73) FACPASS = 9;
		else if (i == 79) FACPASS = 10;
		else if (i == 89) FACPASS = 11;
		else if (i == 97) FACPASS = 12;
		else if (i == 103) FACPASS = 13;
		else if (i == 113) FACPASS = 14;
		else if (i == 119) FACPASS = 15;
		else goto bad;
		FACHSW = fachi;
		FACMSW = facmid;
		FACLSW = p;
		SRCARG = FACDATA;
		setupf ();

/* Factor found, is it a match? */

		do {
			if (factor64 () != 2 &&
			    FACHSW == fachi &&
			    FACMSW == facmid &&
			    FACLSW == faclo) {
				sprintf (buf, "%ld factored OK.\n", p);
				OutputSomewhere (buf);
				goto nextp;
			}
		} while (FACMSW == facmid);

/* Uh oh. */

bad:		sprintf (buf, "%ld factor not found.\n", p);
		OutputBoth (buf);

/* If an escape key was hit, write out the results and return */

nextp:		if (stopCheck ()) break;
		p = 0;
	}

/* All done */

	fclose (fd);
	factorDone ();
}

/* Do the P-1 factoring step prior to a Lucas-Lehmer test */

int pfactor (
	struct work_unit *w)
{
	unsigned long bound1, bound2, squarings;
	double	prob;
	char	buf[120];

/* Deduce the proper P-1 bounds */

	guess_pminus1_bounds (w->p, w->bits,
			      (w->work_type == WORK_DBLCHK ||
			       (w->work_type == WORK_PFACTOR && w->pminus1ed)),
			      &bound1, &bound2, &squarings, &prob);
	if (bound1 == 0) {
		updateWorkToDo (w->p, WORK_PMINUS1, 0);
		return (TRUE);
	}

/* Output a message that P-1 factoring is about to begin */

	sprintf (buf, "Starting P-1 factoring on M%ld with B1=%ld, B2=%ld\n",
		 w->p, bound1, bound2);
	OutputStr (buf);
	sprintf (buf, "Chance of finding a factor is an estimated %.3g%%\n",
		 prob * 100.0);
	OutputStr (buf);

/* Set flag indicating we are using a lot of memory.  Actually, we aren't */
/* in stage 1, but if the memory settings change we want to recompute the */
/* memory bounds. */

	HIGH_MEMORY_USAGE = TRUE;

/* Call the P-1 factoring code */

	return (pminus1 (w->p, bound1, bound1, bound2, FALSE, w->bits));
}
