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
char ERRMSG1D[] = "ERROR: Shift counter corrupt.\n";
char ERRMSG1E[] = "ERROR: Illegal double encountered.\n";
char ERRMSG1F[] = "ERROR: FFT data has been zeroed!\n";
char ERRMSG2[] = "Possible hardware failure, consult the readme.txt file.\n";
char ERRMSG3[] = "Continuing from last save file.\n";
char ERRMSG4[] = "Waiting five minutes before restarting.\n";
char ERRMSG5[] = "For added safety, redoing iteration using a slower, more reliable method.\n";
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
		p == 3021377 || p == 6972593 || p == 13466917 ||
		p == 20996011 || p == 24036583);
}

/* Routine to set the stop timer when available memory settings change */

void memSettingsChanged (void) {
	if (HIGH_MEMORY_USAGE || STOP_TIME) STOP_TIME = 1;
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
		if ((x & 0xFF000000) != 0 && (x & 0xFF000000) != 0xFF000000)
			goto err;
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
	if (*units_bit >= PARG) goto err;
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

	if (WELL_BEHAVED_WORK && work_type == WORK_FACTOR) {
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
	    (w->bits < factorLimit (w->p, w->work_type) || !w->pminus1ed))
		return (TRUE);
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

/* If we're checking for the day/night memory change, check the timer now. */

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
	if (RDTSC_TIMING < 10) {
		timers[i] -= getHighResTimer ();
	} else if (RDTSC_TIMING > 10 && (CPU_FLAGS & CPU_RDTSC)) {
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
	if (RDTSC_TIMING < 10) {
		timers[i] += getHighResTimer ();
	} else if (RDTSC_TIMING > 10 && (CPU_FLAGS & CPU_RDTSC)) {
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
	if (RDTSC_TIMING < 10)
		return (timers[i] / getHighResTimerFrequency ());
	else if (RDTSC_TIMING > 10 && (CPU_FLAGS & CPU_RDTSC))
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
	char	buf[40];

/* The timer could be less than zero if the computer went into hibernation. */
/* Hibernation is where the memory image is saved to disk and the computer */
/* shut off.  Upon power up the memory image is restored but the RDTSC */
/* timestamp counter has been reset to zero. */

	t = timer_value (i);
	if (t < 0.0) {
		strcpy (buf, "Unknown");
		timers[i] = 0.0;
	}

/* Format the timer value in one of several styles */

	else {
		int	style;

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
		if (RDTSC_TIMING == 12 && (CPU_FLAGS & CPU_RDTSC)) {
			sprintf (buf+strlen(buf), " (%.0f clocks)", timers[i]);
		}
	}

/* Append optional newline */

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

static int PAUSED = FALSE;

int pauseWhileRunning ()
{
	int	i;
static	time_t	last_check = 0;
static	int	how_often = -1;

/* In most installations, there is no pause list.  Return quickly. */

	if (PAUSE_WHILE_RUNNING == NULL) return (TRUE);

/* Only check the paused list every 10 seconds */

	if (how_often < 0)
		how_often = IniGetInt (INI_FILE, "PauseCheckInterval", 10);
	if (how_often) {
		time_t	current_time;
		time (&current_time);
		if (current_time < last_check + how_often) return (TRUE);
		last_check = current_time;
	}

/* Check if a process is running such that we should pause */

	if (! checkPauseList ()) return (TRUE);
	PAUSED = TRUE;

/* Every 10 seconds see if we should resume processing.  More frequently */
/* see if we should write a save file and stop processing - most likely */
/* due to a shutdown. */

	for ( ; ; ) {
		for (i = 0; i < 20; i++) {
			Sleep (500);
			if (escapeCheck ()) {
				PAUSED = FALSE;
				return (FALSE);
			}
		}
		if (! checkPauseList ()) {
			OutputStr ("Resuming processing.\n");
			PAUSED = FALSE;
			return (TRUE);
		}
	}
}

int isInPauseList (
	char	*program_name)
{
	char	buf[512];
	char	**p;

	strcpy (buf, program_name);
	strupper (buf);
	for (p = PAUSE_WHILE_RUNNING; *p; p++) {
		if (strstr (buf, *p) != NULL) {
			if (!PAUSED) {
				sprintf (buf,
					 "Pausing because %s is running.\n",
					 program_name);
				OutputStr (buf);
			}
			return (TRUE);
		}
	}
	return (FALSE);
}

/* Prepare for making a factoring run */

void factorSetup (unsigned long p)
{

/* Allocate 1MB memory for factoring */

	if (FACDATA == NULL) FACDATA = malloc (1000000);

/* Call the factoring setup assembly code */

	FACLSW = p;
	SRCARG = FACDATA;
	setupf ();

/* If using the SSE2 factoring code, do more initialization */
/* We need to initialize much of the following data: */
/*	XMM_BITS30		DD	0,3FFFFFFFh,0,3FFFFFFFh
	XMM_INITVAL		DD	0,0,0,0
	XMM_INVFAC		DD	0,0,0,0
	XMM_I1			DD	0,0,0,0
	XMM_I2			DD	0,0,0,0
	XMM_F1			DD	0,0,0,0
	XMM_F2			DD	0,0,0,0
	XMM_F3			DD	0,0,0,0
	XMM_TWO_120_MODF1	DD	0,0,0,0
	XMM_TWO_120_MODF2	DD	0,0,0,0
	XMM_TWO_120_MODF3	DD	0,0,0,0
	XMM_INIT120BS		DD	0,0
	XMM_INITBS		DD	0,0
	XMM_BS			DD	0,0
	XMM_SHIFTER		DD	48 DUP (0)
	TWO_TO_FACSIZE_PLUS_62	DQ	0.0
	SSE2_LOOP_COUNTER	DD	0 */
/* The address to XMM_BITS30 was returned in SRCARG. */

	if (CPU_FLAGS & CPU_SSE2) {
		unsigned long i, bits_in_factor;
		unsigned long *xmm_data;

/* Compute the number of bits in the factors we will be testing */

		if (FACHSW) bits_in_factor = 64, i = FACHSW;
		else if (FACMSW) bits_in_factor = 32, i = FACMSW;
		else return;
		while (i) bits_in_factor++, i >>= 1;

/* Factors 63 bits and below use the non-SSE2 code */

		if (bits_in_factor <= 63) return;

/* Set XMM_SHIFTER values (the first shifter value is not used). */
/* Also compute the initial value. */

		xmm_data = (unsigned long *) SRCARG;
		for (i = 0; p > bits_in_factor + 59; i++) {
			xmm_data[52+i*2] = (p & 1) ? 1 : 0;
			p >>= 1;
		}
		xmm_data[4] =			/* XMM_INITVAL */
		xmm_data[6] = p >= 90 ? 0 : (1 << (p - 60));
		xmm_data[44] = 62 - (120 - bits_in_factor);/* XMM_INIT120BS */
		xmm_data[46] = 62 - (p - bits_in_factor);/* XMM_INITBS */
		xmm_data[100] = i;		/* SSE2_LOOP_COUNTER */
		*(double *)(&xmm_data[98]) = pow (2.0, bits_in_factor + 62);
						/* TWO_TO_FACSIZE_PLUS_62 */

/* Set XMM_BS to 60 - (120 - fac_size + 1) as defined in factor64.mac */

		xmm_data[48] = bits_in_factor - 61;
	}
}

/* Cleanup after making a factoring run */

void factorDone (void)
{

/* Free factoring data */

	free (FACDATA);
	FACDATA = NULL;
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
		char	buf[80];
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

	if (!isPrime (w.p) &&
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
	    (w.p < 727 || w.p > MAX_FACTOR || w.bits >= factorLimit (w.p, w.work_type))) {
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

/* Type 3 files are obsolete Advanced / Factoring continuation files */

		if (type == 3) goto readerr;

/* Type 2 files are factoring continuation files */
/* Type 4 files are Advanced / Factoring continuation files */

		if (type < 2 || type > 4) {

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
	    w.bits < factorLimit (w.p,
			          w.work_type != WORK_PFACTOR ? w.work_type :
			          w.pminus1ed ? WORK_DBLCHK : WORK_TEST)) {
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
		if (! pick_fft_size (w.p)) goto check_stop_code;
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

		if (! pick_fft_size (w.p)) goto check_stop_code;
		lucasSetup (w.p, fftlen ? fftlen : advanced_map_exponent_to_fftlen (w.p));

/* Read the initial data from the file, on failure try the backup */
/* intermediate file. */

		if (fd && fftlen) {
			if (! readFileData (fd, &units_bit, &err_cnt)) {
				lucasDone ();
				goto readerr;
			}

/* Handle case where the save file was for a different FFT length than */
/* we would prefer to use.  This can happen, for example, when upgrading */
/* to a Pentium 4 with its SSE2-based FFT using less precision. */

			if (fftlen != advanced_map_exponent_to_fftlen (w.p)) {
				giant	g;
				g = newgiant ((w.p + 32) / sizeof (short));
				gwtobinary (LLDATA, g);
				lucasDone ();
				lucasSetup (w.p, advanced_map_exponent_to_fftlen (w.p));
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
		if (! pick_fft_size (w.p)) goto check_stop_code;
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
		char	buf[80];
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

	if (STOP_REASON == STOP_MEM_CHANGED) {
		OutputStr ("Restarting with new memory settings.\n");
		goto did_some_work;
	}

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

/* For exponents that are near an FFT limit, do 1000 sample iterations */
/* to see if we should use the smaller or larger FFT size.  We examine */
/* the average roundoff error to determine which FFT size to use. */

int pick_fft_size (
	unsigned long p)
{
	char	buf[120];
	double	softpct, total_error, avg_error, max_avg_error;
	unsigned long small_fftlen, large_fftlen, fftlen;
	int	i;

/* We don't do this for small exponents.  We've not studied the average */
/* error enough on smaller FFT sizes to intelligently pick the FFT size. */
/* Also, for really large exponents there is no larger FFT size to use! */

	if (p <= 5000000) return (TRUE);
	if (p >= 75000000) return (TRUE);

/* Get the info on how what percentage of exponents on either side of */
/* an FFT crossover we will do this 1000 iteration test. */

	IniGetString (INI_FILE, "SoftCrossover", buf, sizeof (buf), "0.2");
	softpct = atof (buf) / 100.0;

/* If this exponent is not close to an FFT crossover, then we are done */

	small_fftlen = map_exponent_to_fftlen (
			(unsigned long) ((1.0 - softpct) * p),
			GW_MERSENNE_MOD);
	large_fftlen = map_exponent_to_fftlen (
			(unsigned long) ((1.0 + softpct) * p),
			GW_MERSENNE_MOD);
	if (small_fftlen == large_fftlen) return (TRUE);

/* If we've already picked an FFT length, then return */

	if (fftlen_from_ini_file (p)) return (TRUE);

/* Let the user be more conservative or more aggressive in picking the */
/* acceptable average error.  By default, we accept an average error */
/* between 0.241 and 0.243 depending on the FFT size. */

	max_avg_error = 0.241 + 0.002 *
		(log (small_fftlen) - log (262144.0)) /
		(log (4194304.0) - log (262144.0));
	IniGetString (INI_FILE, "SoftCrossoverAdjust", buf, sizeof (buf), "0");
	max_avg_error += atof (buf);

/* Print message to let user know what is going on */

	sprintf (buf,
		 "Trying 1000 iterations for exponent %ld using %dK FFT.\n",
		 p, small_fftlen / 1024);
	OutputBoth (buf);
	sprintf (buf,
		 "If average roundoff error is above %.5g, then a larger FFT will be used.\n",
		 max_avg_error);
	OutputBoth (buf);

/* Init the FFT code using the smaller FFT size */

	lucasSetup (p, small_fftlen);

/* Fill data space with random values then do one squaring to make */
/* the data truly random. */

	generateRandomData ();
	gwsetnormroutine (0, TRUE, 0);
	gwstartnextfft (TRUE);
	gwsquare (LLDATA);

/* Average the roundoff error over a 1000 iterations. */

	for (i = 0, total_error = 0.0; ; ) {
		MAXERR = 0.0;
		gwsquare (LLDATA);
		total_error += MAXERR;
		if (escapeCheck ()) {
			lucasDone ();
			return (FALSE);
		}
		if (++i == 1000) break;
		if (i % 100 == 0) {
			sprintf (buf,
				 "After %d iterations average roundoff error is %.5g.\n",
				 i, total_error / (double) i);
			OutputStr (buf);
		}
	}
	avg_error = total_error / 1000.0;
	lucasDone ();

/* Now decide which FFT size to use based on the average error. */
/* Save this info in local.ini so that we don't need to do this again. */
/* We write the SSE2 flag to the INI file so that it won't cause a problem */
/* if the local.ini file is copied between P3 and P4 machines. */

	fftlen = (avg_error <= max_avg_error) ? small_fftlen : large_fftlen;
	sprintf (buf, "%ld,%ld,%d", p, fftlen, CPU_FLAGS & CPU_SSE2 ? 1 : 0);
	IniWriteString (LOCALINI_FILE, "SoftCrossoverData", buf);

/* Output message to user informing him of the outcome. */

	sprintf (buf,
		 "Final average roundoff error is %.5g, using %dK FFT for exponent %ld.\n",
		 avg_error, fftlen / 1024, p);
	OutputBoth (buf);
	return (TRUE);
}

/* Test if we are near the maximum exponent this fft length can test */

int exponent_near_fft_limit ()
{
	unsigned long max_exponent;
	char	pct[30];

	max_exponent = map_fftlen_to_max_exponent (FFTLEN, GW_MERSENNE_MOD);
	IniGetString (INI_FILE, "NearFFTLimitPct", pct, sizeof(pct), "0.5");
	return (PARG > (100.0 - atof (pct)) / 100.0 * max_exponent);
}

/* Do an LL iteration very carefully.  This is done after a normal */
/* iteration gets a roundoff error above 0.40.  This careful iteration */
/* will not generate a roundoff error. */

void careful_iteration (
	gwnum	x,			/* Number to square */
	unsigned long *units_bit)	/* Units bit (if subtracting two) */
{
	gwnum	hi, lo;
	unsigned long i;

/* Copy the data to hi and lo.  Zero out half the FFT data in each. */

	hi = gwalloc ();
	lo = gwalloc ();
	gwcopy (x, hi);
	gwcopy (x, lo);
	for (i = 0; i < FFTLEN/2; i++) set_fft_value (hi, i, 0);
	for ( ; i < FFTLEN; i++) set_fft_value (lo, i, 0);

/* Now do the squaring using three multiplies and adds */

	gwsetnormroutine (0, 0, 0);
	gwstartnextfft (0);
	gwsetaddin (0, 0);
	gwfft (hi, hi);
	gwfft (lo, lo);
	gwfftfftmul (lo, hi, x);
	gwfftfftmul (hi, hi, hi);
	if (units_bit != NULL) lucas_fixup (units_bit);
	gwfftfftmul (lo, lo, lo);
	gwaddquick (x, x);
	gwaddquick (hi, x);
	gwadd (lo, x);

/* Free memory and return */

	gwfree (hi);
	gwfree (lo);
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
	double	*addr1;
	int	priority_work = 0;
	int	escaped, saving, near_fft_limit, sleep5;
	unsigned long i, high32, low32;
	int	isPrime, rc;
	char	buf[160];
	time_t	start_time, current_time;
	unsigned long interimFiles, interimResidues, throttle;
static	unsigned long last_counter = 0;		/* Iteration of last error */
static	int	maxerr_recovery_mode = 0;	/* Big roundoff err rerun */

/* A new option to create interim save files every N iterations. */
/* This allows two machines to simultanously work on the same exponent */
/* and compare results along the way. */

	interimFiles = IniGetInt (INI_FILE, "InterimFiles", 0);
	interimResidues = IniGetInt (INI_FILE, "InterimResidues", interimFiles);

/* Option to slow down the program by sleeping after every iteration.  You */
/* might use this on a laptop or a computer running in a hot room to keep */
/* temperatures down and thus reduce the chance of a hardware error.  */

	throttle = IniGetInt (INI_FILE, "Throttle", 0);

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

/* If we are near the maximum exponent this fft length can test, then we */
/* will error check all iterations */

	near_fft_limit = exponent_near_fft_limit ();

/* Get address of second FFT data element.  We'll use this for very */
/* quickly checking for zeroed FFT data. */

	addr1 = addr (LLDATA, 1);

/* Compute numbers in the lucas series, write out every 30 minutes to a file */

	iters = 0;
	escaped = 0;
	while (counter < p) {
		int	stopping, echk;

/* Use this opportunity to perform other miscellaneous tasks that may */
/* be required by this particular OS port */

		doMiscTasks ();

/* Every so often, communicate with the server. */
/* Even more rarely, set flag to see if we have enough work queued up. */

		if ((counter & 0x3F) == 0) {
			EXP_PERCENT_COMPLETE = (double) counter / (double) p;
			if (!communicateWithServer ()) escaped = 1;
			if (!pauseWhileRunning ()) escaped = 1;

/* See if we should switch to factoring.  This is sometimes done */
/* when we get a new assignment that hasn't been factored.  To be sure */
/* that we will be able to immediately startup a LL test when this one */
/* completes, we will factor the next exponent and then come back to */
/* finish the LL test of this exponent. */

			if ((counter & 0xFFFF) == 0) {
				unsigned long sets;
				double	timing, iters_per_day;
				if (getPriorityWork ()) priority_work = 1;

/* In an effort to reduce wild fluctuations in the rolling average */
/* we will try to update the rolling average only twice a day. */

				timing = map_fftlen_to_timing (
						FFTLEN, GW_MERSENNE_MOD,
						CPU_TYPE, CPU_SPEED);
				iters_per_day = CPU_HOURS * 3600.0 / timing;
				sets = (unsigned long) (iters_per_day / 2.0 / 65536.0);
				if (sets == 0) sets = 1;
				if ((counter >> 16) % sets == 0) {
					updateRollingAverage (
						timing * sets * 65536.0 *
						24.0 / CPU_HOURS);
					startRollingAverage ();
				}
			}
		}

/* Error check the last 50 iterations, before writing an */
/* intermediate file (either user-requested stop or a */
/* 30 minute interval expired), and every 128th iteration. */
/* Also save right after a we pass an errored iteration and several */
/* iterations before retesting an errored iteration so that we don't */
/* have to backtrack very far to do a careful_iteration	(we don't do the */
/* iteration immediately before because on the P4 a save operation will */
/* change the FFT data and make the error non-reproducible. */

		escaped |= stopCheck ();
		saving = (counter == last_counter-8 || counter == last_counter);
		stopping = escaped || priority_work;
		echk = stopping || saving || near_fft_limit || ERRCHK || (counter >= p - 50);
		if ((counter & 127) == 0) {
			echk = 1;
			time (&current_time);
			saving |= (current_time - start_time > write_time);
		}
		MAXERR = 0.0;

/* Do a Lucas-Lehmer iteration */

		start_timer (0);

/* If we are recovering from a big roundoff error, then run one */
/* iteration using three multiplies where half the data is zeroed. */
/* This won't run into any roundoff problems and will protect from */
/* roundoff errors up to 0.6. */

		if (maxerr_recovery_mode && counter == last_counter) {
			careful_iteration (LLDATA, &units_bit);
			maxerr_recovery_mode = 0;
			echk = 0;
		}

/* Otherwise, do a normal iteration */

#ifndef SERVER_TESTING
		else {
			gwsetnormroutine (0, echk, 0);
			gwstartnextfft (!stopping && !saving &&
					!maxerr_recovery_mode &&
					counter+1 != p &&
					(interimResidues == 0 ||
					 (counter+1) % interimResidues > 2));
			lucas_fixup (&units_bit);
//gwnum qw = gwalloc ();
//gwfft (LLDATA, qw);
//gwfftfftmul (qw, qw, LLDATA);
//	gwcopy (LLDATA, qw);
//	gwmul (qw, LLDATA);
//gwfree(qw);
			gwsquare (LLDATA);
			if (throttle) Sleep (throttle);
		}
#endif

/* End iteration timing and increase count of iteration completed */

		end_timer (0);
		iters++;

/* If the sum of the output values is an error (such as infinity) */
/* then raise an error. */

		if (gw_test_illegal_sumout ()) {
			sprintf (buf, ERRMSG0, counter, p, ERRMSG1A);
			OutputBoth (buf);
			inc_error_count (2, &error_count);
			sleep5 = TRUE;
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
				sleep5 = TRUE;
				goto restart;
			}
		}

/* Check for excessive roundoff error.  If round off is too large, repeat */
/* the iteration to see if this was a hardware error.  If it was repeatable */
/* then repeat the iteration using a safer, slower method.  This can */
/* happen when operating near the limit of an FFT. */

		if (echk && MAXERR >= 0.40625) {
			static double last_maxerr = 0.0;
			if (counter == last_counter && MAXERR == last_maxerr) {
				OutputBoth (ERROK);
				inc_error_count (3, &error_count);
				GWERROR = 0;
				OutputBoth (ERRMSG5);
				maxerr_recovery_mode = 1;
				sleep5 = FALSE;
				goto restart;
			} else {
				char	msg[100];
				sprintf (msg, ERRMSG1C, MAXERR);
				sprintf (buf, ERRMSG0, counter, p, msg);
				OutputBoth (buf);
				last_counter = counter;
				last_maxerr = MAXERR;
				inc_error_count (1, &error_count);
				sleep5 = FALSE;
				goto restart;
			}
		}

/* Check if the units_bit is corrupt.  This will make sure we are always */
/* subtracting 2 from the FFT data.  If the FFT data was mysteriously zeroed */
/* and the units_bit value was corrupt then we could get a false positive */
/* result.  With this fix we should get into a safe -2, 2, 2, 2 loop. */

		if (units_bit >= p) {
			sprintf (buf, ERRMSG0, counter, p, ERRMSG1D);
			OutputBoth (buf);
			inc_error_count (2, &error_count);
			sleep5 = TRUE;
			goto restart;
		}

/* Check if the FFT data has been zeroed. This will help reduce the chances */
/* of another false positive being reported. */

		if (*addr1 == 0.0 && p > 1000 &&
		    counter > 50 && counter < p-2 && counter != last_counter) {
			for (i = 2; ; i++) {
				if (*addr (LLDATA, i) != 0.0) break;
				if (i == 50) {
					sprintf (buf, ERRMSG0, counter, p, ERRMSG1F);
					OutputBoth (buf);
					inc_error_count (2, &error_count);
					last_counter = counter;
					sleep5 = TRUE;
					goto restart;
				}
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
				 "M%ld interim WZ%d residue %08lX%08lX at iteration %ld\n",
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
		if (! is_valid_fft_value (LLDATA, i)) {
			sprintf (buf, ERRMSG0, counter, p, ERRMSG1E);
			OutputBoth (buf);
			inc_error_count (2, &error_count);
			sleep5 = TRUE;
			goto restart;
		}
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
		sprintf (buf, "M%ld is prime! WZ%d: %08lX,%08lX\n",
			 p, PORT, SEC1 (p), error_count);
		high32 = low32 = 0;
	} else {
		generateResidue64 (units_bit, &high32, &low32);
		sprintf (buf,
			 "M%ld is not prime. Res64: %08lX%08lX. WZ%d: %08lX,%ld,%08lX\n",
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

	if (!isPrime || isKnownMersennePrime (p)) {
		if (rc) _unlink (filename);
		filename[0] = 'q';
		_unlink (filename);
	}

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

	if (sleep5) OutputBoth (ERRMSG2);
	OutputBoth (ERRMSG3);

/* Update the error count in the save file */

	writeNewErrorCount (filename, error_count);

/* Sleep five minutes before restarting */

	if (sleep5 && ! SleepFive ()) return (FALSE);

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
char SELFFAIL6[] = "Maximum number of warnings exceeded.\n";

#define SELFPASS "Self-test %iK passed!\n"
char SelfTestIniMask[] = "SelfTest%iPassed";

int SELF_TEST_ERRORS = 0;
int SELF_TEST_WARNINGS = 0;

struct self_test_info {
	unsigned long p;
	unsigned long iters;
	unsigned long reshi;
};

#define MAX_SELF_TEST_ITERS	399
struct self_test_info SELF_TEST_DATA[MAX_SELF_TEST_ITERS] = {
{78643201, 400, 0x2D9C8904}, {78643199, 400, 0x7D469182},
{75497473, 400, 0x052C7FD8}, {75497471, 400, 0xCCE7495D},
{71303169, 400, 0x467A9338}, {71303167, 400, 0xBBF8B37D},
{68157441, 400, 0xBE71E616}, {68157439, 400, 0x93A71CC2},
{66060289, 400, 0xF296BB99}, {66060287, 400, 0x649EEF2A},
{62390273, 400, 0xBC8DFC27}, {62390271, 400, 0xDE7D5B5E},
{56623105, 400, 0x0AEBF972}, {56623103, 400, 0x1BA96297},
{53477377, 400, 0x5455F347}, {53477375, 400, 0xCE1C7F78},
{50331649, 400, 0x3D746AC8}, {50331647, 400, 0xE23F2DE6},
{49807361, 400, 0xB43EF4C5}, {49807359, 400, 0xA8BEB02D},
{47185921, 400, 0xD862563C}, {47185919, 400, 0x17281086},
{41943041, 400, 0x0EDA1F92}, {41943039, 400, 0xDE6911AE},
{39845889, 400, 0x43D8A96A}, {39845887, 400, 0x3D118E8F},
{37748737, 400, 0x38261154}, {37748735, 400, 0x22B34CD2},
{35651585, 400, 0xB0E48D2E}, {35651583, 400, 0xCC3340C6},
{34865153, 400, 0xD2C00E6C}, {34865151, 400, 0xFA644F69},
{33030145, 400, 0x83E5738D}, {33030143, 400, 0x6EDBC5B5},
{31195137, 400, 0xFF9591CF}, {31195135, 400, 0x04577C70},
{29884417, 400, 0xACC36457}, {29884415, 400, 0xC0FE7B1E},
{28311553, 400, 0x780EB8F5}, {28311551, 400, 0xE6D128C3},
{26738689, 400, 0x09DC45B0}, {26738687, 400, 0xDC7C074A},
{24903681, 400, 0xA482CF1E}, {24903679, 400, 0x4B3F5121},
{23592961, 400, 0xAFE3C198}, {23592959, 400, 0xCF9AD48C},
{20971521, 400, 0x304EC13B}, {20971519, 400, 0x9C4E157E},
{19922945, 400, 0x83FE36D9}, {19922943, 400, 0x9C60E7A2},
{18874369, 400, 0x83A9F8CB}, {18874367, 400, 0x5A6E22E0},
{17825793, 400, 0xF3A90A5E}, {17825791, 400, 0x6477CA76},
{17432577, 400, 0xCAB36E6A}, {17432575, 400, 0xB8F814C6},
{16515073, 400, 0x91EFCB1C}, {16515071, 400, 0xA0C35CD9},
{15597569, 400, 0x12E057AD}, {15597567, 400, 0xC4EFAEFD},
{14942209, 400, 0x1C912A7B}, {14942207, 400, 0xABA9EA6E},
{14155777, 400, 0x4A943A4E}, {14155775, 400, 0x00789FB9},
{13369345, 400, 0x27A041EE}, {13369343, 400, 0xA8B01A41},
{12451841, 400, 0x4DC891F6}, {12451839, 400, 0xA75BF824},
{11796481, 400, 0xFDD67368}, {11796479, 400, 0xE0237D19},
{10485761, 400, 0x15419597}, {10485759, 400, 0x154D473B},
{10223617, 400, 0x26039EB7}, {10223615, 400, 0xC9DFB1A4},
{9961473, 400, 0x3EB29644}, {9961471, 400, 0xE2AB9CB2},
{9437185, 400, 0x42609D65}, {9437183, 400, 0x77ED0792},
{8716289, 400, 0xCCA0C17B}, {8716287, 400, 0xD47E0E85},
{8257537, 400, 0x80B5C05F}, {8257535, 400, 0x278AE556},
{7798785, 400, 0x55A2468D}, {7798783, 400, 0xCF62032E},
{7471105, 400, 0x0AE03D3A}, {7471103, 400, 0xD8AB333B},
{7077889, 400, 0xC516359D}, {7077887, 400, 0xA23EA7B3},
{6684673, 400, 0xA7576F00}, {6684671, 400, 0x057E57F4},
{6422529, 400, 0xC779D2C3}, {6422527, 400, 0xA8263D37},
{6225921, 400, 0xB46AEB2F}, {6225919, 400, 0xD0A5FD5F},
{5898241, 400, 0xE46E76F9}, {5898239, 400, 0x29ED63B2},
{5505025, 400, 0x83566CC3}, {5505023, 400, 0x0B9CBE64},
{5242881, 400, 0x3CC408F6}, {5242879, 400, 0x0EA4D112},
{4980737, 400, 0x6A2056EF}, {4980735, 400, 0xE03CC669},
{4718593, 400, 0x87622D6B}, {4718591, 400, 0xF79922E2},
{4587521, 400, 0xE189A38A}, {4587519, 400, 0x930FF36C},
{4358145, 400, 0xDFEBF850}, {4358143, 400, 0xBB63D330},
{4128769, 400, 0xC0844AD1}, {4128767, 400, 0x25BDBFC3},
{3932161, 400, 0x7A525A7E}, {3932159, 400, 0xF30C9045},
{3735553, 400, 0xFAD79E97}, {3735551, 400, 0x005ED15A},
{3538945, 400, 0xDDE5BA46}, {3538943, 400, 0x15ED5982},
{3342337, 400, 0x1A6E87E9}, {3342335, 400, 0xECEEA390},
{3276801, 400, 0x3341C77F}, {3276799, 400, 0xACA2EE28},
{3112961, 400, 0x2BDF9D2B}, {3112959, 400, 0xA0AC8635},
{2949121, 400, 0x36EDB768}, {2949119, 400, 0x53FD5473},
{2785281, 400, 0x66816C94}, {2785279, 400, 0x059E8D6B},
{2654999, 400, 0x07EE900D}, {2621441, 400, 0x2BC1DACD},
{2621439, 400, 0xBCBA58F1}, {2653987, 400, 0xB005CACC},
{2651879, 400, 0x38DCD06B}, {2654003, 400, 0x1ED556E7},
{2620317, 400, 0x09DB64F8}, {2539613, 400, 0x4146EECA},
{2573917, 400, 0x939DA3B3}, {2359297, 400, 0x73A131F0},
{2359295, 400, 0x53A92203}, {2646917, 400, 0x71D4E5A2},
{2605473, 400, 0xE11637FC}, {2495213, 400, 0x89D80370},
{2540831, 400, 0x2CF01FBB}, {2654557, 400, 0x4106F46F},
{2388831, 400, 0xA508B5A7}, {2654777, 400, 0x9E744AA3},
{2584313, 400, 0x800E9A61}, {2408447, 400, 0x8C91E8AA},
{2408449, 400, 0x437ECC01}, {2345677, 400, 0x60AEE9C2},
{2332451, 400, 0xAB209667}, {2330097, 400, 0x3FB88055},
{2333851, 400, 0xFE4ECF19}, {2444819, 400, 0x56BF33C5},
{2555671, 400, 0x9DC03527}, {2654333, 400, 0xE81BCF40},
{2543123, 400, 0x379CA95D}, {2432123, 400, 0x5952676A},
{2321123, 400, 0x24DCD25F}, {2654227, 400, 0xAC3B7F2B},
{2329999, 400, 0xF5E902A5}, {2293761, 400, 0x9E4BBB8A},
{2293759, 400, 0x1901F07B}, {2236671, 400, 0x45EB162A},
{2193011, 400, 0x382B6E4B}, {2329001, 400, 0x4FF052BB},
{2327763, 400, 0x3B315213}, {2325483, 400, 0x0DC5165A},
{2323869, 400, 0xD220E27F}, {2315679, 400, 0xF650BE33},
{2004817, 400, 0xC2FF3440}, {2130357, 400, 0xC25804D8},
{2288753, 400, 0xA4DD9AAD}, {2266413, 400, 0x675257DB},
{2244765, 400, 0xC08FF487}, {2222517, 400, 0x1A128B22},
{2200339, 400, 0x0EB0E827}, {2328117, 400, 0x0A24673A},
{2329557, 400, 0x2E267692}, {2188001, 400, 0xD012AF6A},
{2166567, 400, 0x509BA41A}, {2144651, 400, 0x54CFC0E6},
{2122923, 400, 0xA47068E6}, {2100559, 400, 0xACFAB4E1},
{2088461, 400, 0xEA01E860}, {2066543, 400, 0x847DF0D0},
{2044767, 400, 0x04225888}, {2022823, 400, 0x6EA34B32},
{2328527, 400, 0xC55E3E05}, {2327441, 400, 0x207C8CEC},
{2326991, 400, 0x0A4F2ACD}, {2009987, 400, 0xE6A59DEF},
{1999999, 400, 0xD645A18F}, {1966081, 400, 0xB88828A1},
{1966079, 400, 0x5BD87C45}, {1998973, 400, 0xCBDD74F7},
{1997651, 400, 0x666B0CB1}, {1675001, 400, 0x50A94DB7},
{1977987, 400, 0x30D1CD1F}, {1955087, 400, 0x5B9426A4},
{1933071, 400, 0x23C1AF0B}, {1911957, 400, 0xF7699248},
{1899247, 400, 0x11C76E04}, {1877431, 400, 0xA3299B39},
{1855067, 400, 0x35243683}, {1833457, 400, 0xCF630DC0},
{1811987, 400, 0x7C7022EC}, {1799789, 400, 0xEFEC47B7},
{1777773, 400, 0x0F16E2D6}, {1755321, 400, 0x1AC5D492},
{1733333, 400, 0x5DA0555E}, {1711983, 400, 0xDC19DA8B},
{1699779, 400, 0x2B44914E}, {1677323, 400, 0x03D3980B},
{1995091, 400, 0x922E555B}, {1993041, 400, 0x0CA8451B},
{1991991, 400, 0xDFFB212D}, {1679779, 400, 0x51D75E0F},
{1684993, 400, 0x048BBCE8}, {1970009, 400, 0x646E0DFA},
{1957445, 400, 0xC8D244ED}, {1999997, 400, 0x5FC899D0},
{1998983, 400, 0x1CD518AA}, {1999007, 400, 0xA9DD8591},
{1674999, 400, 0xDB0169D8}, {1638401, 400, 0xD3F8A8C5},
{1638399, 400, 0xF270D8DD}, {1674997, 400, 0xC824EF15},
{1674551, 400, 0xD844AEAD}, {1674001, 400, 0x8F5EFA50},
{1345001, 400, 0x18EE2E2D}, {1655083, 400, 0x09B30DEE},
{1633941, 400, 0x0B87C8B1}, {1611557, 400, 0x6B57E48D},
{1599549, 400, 0x48EA38B2}, {1577771, 400, 0xCE84D9DC},
{1555947, 400, 0x6797EEF4}, {1533349, 400, 0xD6897409},
{1511861, 400, 0x8A8177AC}, {1499625, 400, 0x56BB6FB3},
{1477941, 400, 0xF3DD8ED3}, {1455931, 400, 0x31A222C7},
{1433069, 400, 0x28F01E1B}, {1411747, 400, 0x680C6E39},
{1399449, 400, 0xB7F01A54}, {1377247, 400, 0xE656F652},
{1355991, 400, 0xB2AA2819}, {1350061, 400, 0x31F9A728},
{1673881, 400, 0xA51D38E4}, {1672771, 400, 0x5474B6F9},
{1671221, 400, 0x2710DDEA}, {1670551, 400, 0x31FC3838},
{1660881, 400, 0x4C5B22C5}, {1650771, 400, 0x998F747B},
{1655001, 400, 0x164659A6}, {1674339, 400, 0xED2D23E2},
{1344999, 400, 0x158AA064}, {1310721, 400, 0x5694A427},
{1310719, 400, 0x258BDDE3}, {1344997, 400, 0x1D059D4F},
{1344551, 400, 0x60606AA3}, {1344001, 400, 0x9AC6AB36},
{1322851, 400, 0x3A000D0A}, {1300993, 400, 0x77CB0184},
{1288771, 400, 0x7431D9E2}, {1266711, 400, 0xB4BC4E8D},
{1244881, 400, 0x48BC9FF9}, {1222991, 400, 0x3F5FC39E},
{1200881, 400, 0xD5DF4944}, {1188441, 400, 0xD9D8968B},
{1166661, 400, 0xD4AB97F4}, {1144221, 400, 0x9940943B},
{1122001, 400, 0x647406B8}, {1100881, 400, 0x3AD40CE0},
{1088511, 400, 0xD578BB51}, {1066837, 400, 0x2F82BFBB},
{1044811, 400, 0x7C6EDDD1}, {1022991, 400, 0x6A1C2DD4},
{1000001, 400, 0x2879748F}, {1343881, 400, 0xB59E8006},
{1342771, 400, 0x87563FFE}, {1341221, 400, 0x29AD6127},
{1340551, 400, 0x17DB4ACB}, {1330881, 400, 0x9642F068},
{942079, 1000, 0xE528A9B0}, {974849, 1000, 0x79791EDB},
{983041, 1000, 0x29216C43}, {901121, 1000, 0x26C4E660},
{917503, 1000, 0x5F244685}, {933889, 1000, 0x62490F57},
{851967, 1000, 0x331AA906}, {860161, 1000, 0x41185F27},
{884735, 1000, 0x7BC7A661}, {802817, 1000, 0xA9645693},
{819199, 1000, 0x48AFB0A5}, {835585, 1000, 0x706437D3},
{753663, 1000, 0x99C43F31}, {778241, 1000, 0x1729A6C4},
{786431, 1000, 0x61080929}, {720897, 1000, 0x1E96863D},
{737279, 1000, 0x1B07A764}, {745473, 1000, 0x7BCE80AA},
{655359, 1000, 0x1107F161}, {659457, 1000, 0x589C16A4},
{688127, 1000, 0xD01E5A85}, {622593, 1000, 0x26F6FC8C},
{630783, 1000, 0x4DD2E603}, {638977, 1000, 0xC88F34B4},
{589823, 1000, 0x0290B60B}, {602113, 1000, 0xEFCD5BA8},
{614399, 1000, 0x6408F880}, {557057, 1000, 0xC30FE589},
{565247, 1000, 0xF4CA3679}, {573441, 1000, 0xF8F039AA},
{532479, 1000, 0x0072FE03}, {540673, 1000, 0xDA0E0D99},
{544767, 1000, 0x62443C6B}, {491521, 1000, 0x3F520DFA},
{516095, 1000, 0xA6BD9423}, {524289, 1000, 0xCD591388},
{466943, 1000, 0xE10EE929}, {471041, 1000, 0x18752F40},
{487423, 1000, 0x933FFF17}, {442369, 1000, 0xC22471C3},
{450559, 1000, 0x025B1320}, {458753, 1000, 0xE296CC00},
{417791, 1000, 0x080C803C}, {425985, 1000, 0xB2095F04},
{430079, 1000, 0x98B1EC61}, {393217, 1000, 0x26DD79ED},
{401407, 1000, 0x2F0F75F9}, {409601, 1000, 0xAEFAC2F8},
{372735, 1000, 0xCB6D00A2}, {376833, 1000, 0x915D5458},
{389119, 1000, 0x6188E38D}, {344065, 1000, 0x4D0C5089},
{360447, 1000, 0x84AC5CFD}, {368641, 1000, 0x72414364},
{319487, 1000, 0x24ED1BE9}, {327681, 1000, 0x3101106A},
{329727, 1000, 0x5BDB69AF}, {307201, 1000, 0x68536CD1},
{311295, 1000, 0x69778074}, {315393, 1000, 0x429D4950},
{286719, 1000, 0x1A31A686}, {294913, 1000, 0xF55727C6},
{301055, 1000, 0x33BDB242}, {272385, 1000, 0xEF6EC4B4},
{278527, 1000, 0x05530FD5}, {282625, 1000, 0x34A4E699},
{262143, 1000, 0xA9638844}, {266241, 1000, 0xE0969CED},
{270335, 1000, 0x14AD54BE}, {243713, 1000, 0xC19AEA91},
{245759, 1000, 0x7538BF0B}, {258049, 1000, 0x73F541AD},
{229375, 1000, 0x6E42B26A}, {233473, 1000, 0x1964F897},
{235519, 1000, 0x661BBC3F}, {215041, 1000, 0x04D5D2F0},
{221183, 1000, 0xA89E7764}, {225281, 1000, 0x20876BED},
{204799, 1000, 0xD20C2126}, {208897, 1000, 0x9D4DCF0E},
{212991, 1000, 0x1FF00E2A}, {194561, 1000, 0x6ED1CB70},
{196607, 1000, 0x3190D5F5}, {200705, 1000, 0xFAD28F5A},
{184319, 1000, 0x360EF08E}, {186369, 1000, 0x0F001482},
{188415, 1000, 0x86FCE4D6}, {164865, 1000, 0x4942B002},
{172031, 1000, 0xC5AF29DB}, {180225, 1000, 0x35D49D74},
{157695, 1000, 0x5422FACF}, {159745, 1000, 0xB5CD03A1},
{163839, 1000, 0x1CA6048E}, {150529, 1000, 0x7412F09C},
{153599, 1000, 0xA9FAAE69}, {155649, 1000, 0xA7B736AF},
{141311, 1000, 0x7A5D0730}, {143361, 1000, 0x580F4DC4},
{147455, 1000, 0x176B299A}, {135169, 1000, 0x65AC10A4},
{136191, 1000, 0xC4591D37}, {139265, 1000, 0xBCE1FC80},
{129023, 1000, 0xAFE1E7A8}, {131073, 1000, 0xC5AAB12F},
{133119, 1000, 0xDE51C35A}, {117761, 1000, 0x054A26F6},
{121855, 1000, 0x55AF2385}, {122881, 1000, 0x652827AC},
{112639, 1000, 0x6FA4DB24}, {114689, 1000, 0x0BBAF161},
{116735, 1000, 0xB85F0E8E}, {106497, 1000, 0xF833D925},
{107519, 1000, 0x80F177D8}, {110593, 1000, 0x1A56AA86},
{100351, 1000, 0x1DE12CE6}, {102401, 1000, 0x19F967B4},
{104447, 1000, 0xF9F3CDFD}
};

#define MAX_SELF_TEST_ITERS2	370
struct self_test_info SELF_TEST_DATA2[MAX_SELF_TEST_ITERS2] = {
{77497473, 900, 0xF0B43F54}, {76497471, 900, 0xF30AFA95},
{75497473, 900, 0x32D8D3A7}, {75497471, 900, 0x9E689331},
{74497473, 900, 0xD43166A4}, {73497471, 900, 0x639E4F0C},
{72303169, 900, 0x74BDED5C}, {71303169, 900, 0xA2147B5C},
{71303167, 900, 0x717525AB}, {70303167, 900, 0xD716B4F0},
{68060289, 1000, 0xF90C7BFF}, {67060287, 1000, 0xFE9BF47C},
{66060289, 1000, 0x057C60F5}, {66060287, 1000, 0x2ECC97CE},
{65390273, 1000, 0xC55C6369}, {64390271, 1000, 0x48552448},
{63390273, 1000, 0x6FF8CD84}, {62390273, 1000, 0x42ACEB15},
{62390271, 1000, 0x48764DF8}, {61390271, 1000, 0xD5408698},
{57623105, 1200, 0x098B4491}, {56623105, 1200, 0x5E720717},
{56623103, 1200, 0x1980D8BC}, {55623103, 1200, 0xEDD592B6},
{53477377, 1200, 0xBAEF5CCC}, {53477375, 1200, 0x2F296FC8},
{52331647, 1200, 0xA1EAE85D}, {51331649, 1200, 0xE3B39845},
{50331649, 1200, 0x53543DF2}, {50331647, 1200, 0x0049E54B},
{48185921, 1500, 0x78F4AEAA}, {47185921, 1500, 0x4D7FFDDC},
{47185919, 1500, 0x059D196F}, {46185919, 1500, 0x38B1D9AD},
{45943041, 1500, 0x7670FDDF}, {44943039, 1500, 0xA859BBD7},
{43943041, 1500, 0xD673E000}, {42943039, 1500, 0x6B69D8CE},
{41943041, 1500, 0x6E92CE47}, {41943039, 1500, 0x888BEE79},
{39151585, 1900, 0x3B06496C}, {38748737, 1900, 0x6429E0FD},
{38251583, 1900, 0x04AD7F99}, {37748737, 1900, 0x47659BC5},
{37748735, 1900, 0x2DFA41B0}, {36748735, 1900, 0x1A1DA557},
{36251585, 1900, 0x83F23FA8}, {35651585, 1900, 0x3598B4B9},
{35651583, 1900, 0x7E443962}, {35251583, 1900, 0x1CE4D084},
{34230145, 2100, 0x0FDE9717}, {33730143, 2100, 0x54EB5333},
{33030145, 2100, 0xF37897B8}, {33030143, 2100, 0x52B3981B},
{32595137, 2100, 0xA76D0805}, {32095135, 2100, 0xCF443ACD},
{31595137, 2100, 0xA6DEA70A}, {31195137, 2100, 0x0777442D},
{31195135, 2100, 0x9B265F8F}, {30695135, 2100, 0xA3BC760F},
{29311553, 2500, 0xFD1D6D74}, {28811551, 2500, 0xE720BFD3},
{28311553, 2500, 0xA11F75AB}, {28311551, 2500, 0x7E0471E5},
{27738689, 2500, 0xD246DC55}, {27238687, 2500, 0x806A3A62},
{26738689, 2500, 0x8E8450B1}, {26738687, 2500, 0xD4A0DBC9},
{26138689, 2500, 0x47C47755}, {25638687, 2500, 0x7E9C7E8E},
{24903681, 3100, 0x50835AB8}, {24903679, 3100, 0xAE3D2F94},
{24092961, 3100, 0x7B540B4D}, {23892959, 3100, 0xA0D4EC50},
{23592961, 3100, 0x47FBD6FE}, {23592959, 3100, 0x09FD89AB},
{22971521, 3100, 0x99DFEDB9}, {21871519, 3100, 0x35A8B46A},
{20971521, 3100, 0x94C12572}, {20971519, 3100, 0x1F6D3003},
{19922945, 4000, 0x86B106EB}, {19922943, 4000, 0xE1CE3C1A},
{19374367, 4000, 0xD1045A66}, {19174369, 4000, 0x3247CE82},
{18874369, 4000, 0x33BB2689}, {18874367, 4000, 0x6856F21F},
{18474367, 4000, 0x95E2F6FA}, {18274367, 4000, 0x61182009},
{18274369, 4000, 0xB2FD8175}, {18074369, 4000, 0x7F242A6E},
{17432577, 4500, 0x632CAD0B}, {17432575, 4500, 0xC9C79F07},
{17115073, 4500, 0xF2B70D4B}, {16815071, 4500, 0x71B22529},
{16515073, 4500, 0xAB1CC854}, {16515071, 4500, 0xF54D05D7},
{16297569, 4500, 0x6B5F72DA}, {15997567, 4500, 0x9669F188},
{15597569, 4500, 0x352BFCCF}, {15597567, 4500, 0x36B164ED},
{14942209, 5300, 0xEA5DB53B}, {14942207, 5300, 0x6CC650A2},
{14155777, 5300, 0xEB7C125D}, {14155775, 5300, 0xB4C8B09B},
{13969343, 5300, 0x832359A5}, {13669345, 5300, 0x7EE99140},
{13369345, 5300, 0xCDF43471}, {13369343, 5300, 0x343FEA12},
{13069345, 5300, 0x65B17A9B}, {12969343, 5300, 0x063F492B},
{12451841, 6500, 0xCB168E5D}, {12451839, 6500, 0xE91EEB5A},
{12196481, 6500, 0x0A261B7E}, {11796481, 6500, 0x38100A5F},
{11796479, 6500, 0x78FCF8C5}, {11596479, 6500, 0x8C481635},
{11285761, 6500, 0x2580BC8D}, {10885759, 6500, 0x54030992},
{10485761, 6500, 0x054660AA}, {10485759, 6500, 0x50F74AF0},
{9961473, 7800, 0x7991161C}, {9961471, 7800, 0x627F3BEE},
{9837183, 7800, 0xBC67A608}, {9737185, 7800, 0x9A0CBC59},
{9537183, 7800, 0xA6A509A6}, {9437185, 7800, 0x877C09B6},
{9437183, 7800, 0x1D259540}, {9337185, 7800, 0x5EF3F14C},
{9237183, 7800, 0x5780245F}, {9137185, 7800, 0x6C1162A9},
{8716289, 9000, 0x2011133F}, {8716287, 9000, 0xEEEC1181},
{8516289, 9000, 0xF1D93A69}, {8316287, 9000, 0x53D6E3CB},
{8257537, 9000, 0x38DB98D6}, {8257535, 9000, 0x7D1BECA7},
{8098785, 9000, 0x51E9FA27}, {7998783, 9000, 0xF7F14FF2},
{7798785, 9000, 0x8437BC4D}, {7798783, 9000, 0x9E28D8E1},
{7471105, 11000, 0xEFDA89EA}, {7471103, 11000, 0x4061C4BF},
{7377889, 11000, 0x65ABE846}, {7277887, 11000, 0x02B0EBD7},
{7077889, 11000, 0x336E1030}, {7077887, 11000, 0x685B792E},
{6984673, 11000, 0x3AE19FAF}, {6884671, 11000, 0x2A0ED16A},
{6684673, 11000, 0x206A3512}, {6684671, 11000, 0x4FD9980A},
{6225921, 13000, 0x1A922371}, {6225919, 13000, 0xC0F63BD8},
{6198241, 13000, 0xDA664501}, {6098239, 13000, 0xB92015CD},
{5898241, 13000, 0xDA384BD9}, {5898239, 13000, 0x20B59AC8},
{5705025, 13000, 0x941A2DA0}, {5605023, 13000, 0xCFDF5835},
{5505025, 13000, 0x37A6C972}, {5505023, 13000, 0x6252AB5C},
{5120737, 17000, 0x512705D0}, {5030735, 17000, 0x633E3E74},
{4980737, 17000, 0xD8245D49}, {4980735, 17000, 0xFB2C3530},
{4888593, 17000, 0xE3C6EDBC}, {4818591, 17000, 0x89E7FE48},
{4718593, 17000, 0xA23C713D}, {4718591, 17000, 0xC7BA41D6},
{4698593, 17000, 0xA0194103}, {4648591, 17000, 0xD5A50A23},
{4501145, 19000, 0x7BAF4344}, {4458143, 19000, 0x686F6B13},
{4358145, 19000, 0x682E6643}, {4358143, 19000, 0x974DA6CC},
{4298769, 19000, 0x1FC0E577}, {4228767, 19000, 0x46B5F3CD},
{4128769, 19000, 0x59332478}, {4128767, 19000, 0x4AF5C8B8},
{4028769, 19000, 0x542C17CB}, {3978767, 19000, 0x76E41351},
{3835553, 22000, 0x9058FE40}, {3785551, 22000, 0x45EF5C15},
{3735553, 22000, 0x2700B350}, {3735551, 22000, 0x09EDCEAD},
{3688945, 22000, 0x626C29D3}, {3618943, 22000, 0x82B1D4D1},
{3538945, 22000, 0x70331CC6}, {3538943, 22000, 0x00FEB746},
{3342337, 22000, 0x7CEE24AE}, {3342335, 22000, 0x1802D072},
{3242961, 27000, 0xE877F863}, {3172959, 27000, 0x04C9F1F7},
{3112961, 27000, 0x241E93DB}, {3112959, 27000, 0x8D359307},
{2949121, 27000, 0x6B545E09}, {2949119, 27000, 0xAFD6F417},
{2885281, 27000, 0x439E57E6}, {2785281, 27000, 0xB4E40DFE},
{2785279, 27000, 0x3787D3FA}, {2685279, 27000, 0x902967B7},
{2605473, 34000, 0xE21C344E}, {2584313, 34000, 0xFDBCFCB2},
{2573917, 34000, 0x89B5012C}, {2540831, 34000, 0x201BAA90},
{2539613, 34000, 0x2226BA6B}, {2495213, 34000, 0xE3577D9F},
{2408447, 34000, 0x594C9155}, {2388831, 34000, 0x55CE9F16},
{2359297, 34000, 0x09A72A40}, {2359295, 34000, 0x621E8BF9},
{2244765, 39000, 0xEC2F362D}, {2236671, 39000, 0x4B50CA20},
{2222517, 39000, 0x8DA427C0}, {2193011, 39000, 0xD1DE8993},
{2130357, 39000, 0x4B5EBB90}, {2122923, 39000, 0x5F9110FC},
{2100559, 39000, 0xE0CF8904}, {2088461, 39000, 0x26AD1DEA},
{2066543, 39000, 0xB78C9237}, {2004817, 39000, 0x3D7838F8},
{1933071, 46000, 0x86323D21}, {1911957, 46000, 0x500CFEAD},
{1899247, 46000, 0x128667DF}, {1877431, 46000, 0x2A59B6B5},
{1855067, 46000, 0xBE9AABF5}, {1833457, 46000, 0xB84D7929},
{1777773, 46000, 0x771E0A9D}, {1755321, 46000, 0xF93334E3},
{1699779, 46000, 0x07B46DEE}, {1677323, 46000, 0x910E0320},
{1633941, 56000, 0x455509CD}, {1611557, 56000, 0x0F51FA1E},
{1599549, 56000, 0x646A96B0}, {1577771, 56000, 0xA4A21303},
{1555947, 56000, 0x80B84725}, {1533349, 56000, 0x23E9F7B1},
{1477941, 56000, 0x593F208F}, {1455931, 56000, 0x11002C52},
{1433069, 56000, 0x5B641D8B}, {1411747, 56000, 0x5EAE18A8},
{1322851, 75000, 0xD5C50F2E}, {1310721, 75000, 0x855E44A2},
{1310719, 75000, 0xC0836C1F}, {1300993, 75000, 0xF62263D6},
{1288771, 75000, 0x867EBBAB}, {1266711, 75000, 0xBA1FF3BE},
{1244881, 75000, 0xCE8199EB}, {1222991, 75000, 0xCDE49EF5},
{1200881, 75000, 0xC8610F6C}, {1188441, 75000, 0xFC772495},
{1150221, 84000, 0xA3334541}, {1144221, 84000, 0x44307B03},
{1122001, 84000, 0x9B937DCF}, {1108511, 84000, 0x9F3D191E},
{1100881, 84000, 0xBAF4EA2D}, {1096837, 84000, 0xAA9396F1},
{1088511, 84000, 0xB0CB2704}, {1066837, 84000, 0x031F202C},
{1044811, 84000, 0x7EA89CFE}, {1022991, 84000, 0xD42294C8},
{983041, 100000, 0x4052BBC0}, {974849, 100000, 0xB0E9EB07},
{942079, 100000, 0xEE230987}, {933889, 100000, 0x58FA63B0},
{917503, 100000, 0x8B457209}, {901121, 100000, 0xD2325FC4},
{884735, 100000, 0xCBB5A603}, {860161, 100000, 0xBC240C77},
{854735, 100000, 0xE8BE766D}, {851967, 100000, 0x09AD9B74},
{827279, 120000, 0x64B01894}, {819199, 120000, 0xF97F1E2B},
{802817, 120000, 0xC4EDBC3C}, {795473, 120000, 0x046584E0},
{786431, 120000, 0xC6BA553D}, {778241, 120000, 0x856A5147},
{753663, 120000, 0xC7895B4A}, {745473, 120000, 0x42B47EA2},
{737279, 120000, 0x29E477B8}, {720897, 120000, 0x97111FA7},
{662593, 160000, 0x32472A99}, {659457, 160000, 0xEF49D340},
{655359, 160000, 0x75C12C38}, {644399, 160000, 0xDE632783},
{638977, 160000, 0xDCDB98B4}, {630783, 160000, 0x6B8F0706},
{622593, 160000, 0xD732286D}, {614399, 160000, 0x2489EFB3},
{612113, 160000, 0xCAE00EC6}, {602113, 160000, 0x792AD67D},
{580673, 180000, 0xC508CAFA}, {573441, 180000, 0xB0680C2B},
{565247, 180000, 0xF1DBB762}, {557057, 180000, 0x374F647B},
{544767, 180000, 0x3DC41F49}, {540673, 180000, 0x949A4CB7},
{532479, 180000, 0xEA06DC97}, {524289, 180000, 0xA76CE14A},
{522479, 180000, 0xAA8EAC14}, {516095, 180000, 0x04F0CC23},
{501041, 210000, 0xD9F72F62}, {496943, 210000, 0xD62D5380},
{487423, 210000, 0x55ACB2FD}, {471041, 210000, 0xB6AEAB0E},
{466943, 210000, 0x251CDE78}, {458753, 210000, 0xDC40CADB},
{450559, 210000, 0x2AD0CF72}, {442369, 210000, 0x5FF2E46E},
{441041, 210000, 0x1194CC23}, {436943, 210000, 0x0272AF35},
{420217, 270000, 0xD233852A}, {409601, 270000, 0x6F89825C},
{401407, 270000, 0x3D9DE818}, {393217, 270000, 0xDE8E6FF0},
{392119, 270000, 0x30CA58B7}, {389119, 270000, 0x80975797},
{376833, 270000, 0xC75824DB}, {372735, 270000, 0xF8BE0932},
{368641, 270000, 0xA48AC5E3}, {360447, 270000, 0x7DD29C13},
{339487, 340000, 0xA7311A6D}, {335393, 340000, 0xD9704DF2},
{331681, 340000, 0x3316A003}, {329727, 340000, 0xE46D5991},
{327681, 340000, 0xBEDA4A7B}, {319487, 340000, 0xB25C84FF},
{315393, 340000, 0xF5AD1DDA}, {311295, 340000, 0xFE41A12A},
{308295, 340000, 0x03AAC47E}, {307201, 340000, 0xFC08ACCC},
{291913, 380000, 0xC56AB884}, {286719, 380000, 0x248EF622},
{282625, 380000, 0x50A98488}, {280335, 380000, 0x9B64A843},
{278527, 380000, 0x39D5B7DB}, {274335, 380000, 0x48623B41},
{270335, 380000, 0xC04B857A}, {266241, 380000, 0xFE4475F6},
{262143, 380000, 0xADC3ECE9}, {260335, 380000, 0x15B8F9EF},
{250519, 460000, 0xA2FE3B50}, {245759, 460000, 0xC6D800D6},
{245281, 460000, 0x4F23AA34}, {243713, 460000, 0xB30EC823},
{235519, 460000, 0x31FD709E}, {233473, 460000, 0x8FCC69C2},
{231183, 460000, 0xD59255CC}, {229375, 460000, 0x788520D0},
{225281, 460000, 0xD669C8BC}, {221183, 460000, 0x9B915F4B},
{212991, 560000, 0x0555250D}, {210415, 560000, 0x3FC3CCD7},
{208897, 560000, 0x9FF8F462}, {204799, 560000, 0x294EB549},
{200705, 560000, 0x80B1222F}, {196607, 560000, 0x8AB8D945},
{194561, 560000, 0x4140E623}, {188415, 560000, 0xFA0A3453},
{186369, 560000, 0xAC17EAB6}, {184319, 560000, 0x835F341B},
{172031, 800000, 0xF6BD0728}, {163839, 800000, 0x26C78657},
{159745, 800000, 0x6ACBB961}, {157695, 800000, 0x3EA979F3},
{155649, 800000, 0x09C7ADE4}, {153599, 800000, 0xF601EB92},
{147455, 800000, 0x0AA97D21}, {143361, 800000, 0xEA6A01F1},
{141311, 800000, 0x9BB8A6A3}, {135169, 800000, 0xECA55A45}
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
#define SELF_FFT_LENGTHS	37
		int	lengths[SELF_FFT_LENGTHS] = {1024,8,10,896,768,12,14,640,512,16,20,448,384,24,28,320,256,32,40,224,192,48,56,160,128,64,80,112,96,1280,1536,1792,2048,2560,3072,3584,4096};
		int	min_fft, max_fft, test_time;
		time_t	start_time, current_time;
		unsigned int memory;	/* Memory to use during torture test */
		void	*bigbuf = NULL;

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
		memory = IniGetInt (INI_FILE, "TortureMem", DAY_MEMORY > NIGHT_MEMORY ? DAY_MEMORY : NIGHT_MEMORY);
		while (memory > 8 && bigbuf == NULL) {
			bigbuf = malloc (memory * 1000000);
			if (bigbuf == NULL) memory--;
		}
		for ( ; ; ) {
		    for (i = 0; i < SELF_FFT_LENGTHS; i++) {
			if (lengths[i] < min_fft || lengths[i] > max_fft)
				continue;
			if (! selfTestInternal (lengths[i]*1024, test_time, i, memory, bigbuf)) {
				if (pArg == 1) {
					char	buf[100];
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
				GW_BIGBUF = NULL;
				GW_BIGBUF_SIZE = 0;
				free (bigbuf);
				return (FALSE);
			}
		    }
		    if (pArg == 0) break;
		}
		GW_BIGBUF = NULL;
		GW_BIGBUF_SIZE = 0;
		free (bigbuf);
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

		if (! selfTestInternal (fftlen, 60, -1, 0, NULL))
			return (FALSE);
	}

/* Self test completed! */

	return (TRUE);
}

int selfTestInternal (
	unsigned long fftlen,
	unsigned int test_time,	/* Number of minutes to self-test */
	int	torture_count,	/* Index into array of fft sizes */
	unsigned int memory,	/* MB of memory the torture test can use */
	void	*bigbuf)	/* Memory block for the torture test */
{
	unsigned long k, limit;
	unsigned int i, iter, countdown;
	char	filename[16];
	char	buf[120];
	char	iniName[32];
	time_t	start_time, current_time;
static	int	data_index[SELF_FFT_LENGTHS] = {0};
	struct self_test_info *test_data;
	unsigned int test_data_count;

/* Set the title */

	title ("Self-Test");

/* Generate file name for temp files */

	strcpy (filename, "ptemp");
	strcat (filename, EXTENSION);

/* Pick which self test data array to use.  Machines are much faster now */
/* compared to when the torture test was introduced.  This new self test */
/* data will run more iterations and thus stress the cpu more by spending */
/* less time in the initialization code. */

	if (CPU_SPEED < 1000.0) {
		test_data = SELF_TEST_DATA;
		test_data_count = MAX_SELF_TEST_ITERS;
	} else {
		test_data = SELF_TEST_DATA2;
		test_data_count = MAX_SELF_TEST_ITERS2;
	}

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
		unsigned int ll_iters, num_gwnums;
		gwnum	*gwarray, g;

/* Find next self test data entry to work on */

		for ( ; ; i++) {

/* Wrap in the self test data array */

			if (i == test_data_count) i = 0;

/* Now select the actual exponent */

			p = test_data[i].p;
			if (p > limit) continue;

/* The SSE2 carry propagation code gets into trouble if there are too */
/* few bits per FFT word!  Thus, we'll require at least 8 bits per */
/* word here.  Now that the number of iterations changes for each FFT */
/* length I'm raising the requirement to 10 bits to keep timings roughly */
/* equal. */

			if (p / fftlen < 10) continue;

/* We've found an exponent to test! */

			break;
		}

/* Output start message */

		ll_iters = test_data[i].iters;
		sprintf (buf, SELF1, iter, ll_iters, p, fftlen / 1024);
		OutputStr (buf);

/* Now run Lucas setup, for extra safety double the maximum allowable */
/* sum(inputs) vs. sum(outputs) difference. */

		GW_BIGBUF = bigbuf;
		GW_BIGBUF_SIZE = (bigbuf != NULL) ? memory * 1000000 : 0;
		lucasSetup (p, fftlen);
		MAXDIFF *= 2.0;
		if (LLDATA == NULL) {		/* Shouldn't ever happen... */
			OutputStr ("Out of memory\n");
			lucasDone ();
			return (FALSE);
		}

/* Determine how many gwnums we can allocate in the memory we are given */

		if (memory <= 8 || (iter & 1) == 0)
			num_gwnums = 1;
		else {
			num_gwnums = (unsigned int)
				(((double) memory * 1000000.0 -
				  (double) map_fftlen_to_memused (FFTLEN, 0)) /
				 (double) (gwnum_size (FFTLEN) + 32 + GW_ALIGNMENT));
			if (num_gwnums < 1) num_gwnums = 1;
			if (num_gwnums > ll_iters) num_gwnums = ll_iters;
		}

/* Allocate gwnums to eat up the available memory */

		gwarray = (gwnum *) malloc (num_gwnums * sizeof (gwnum));
		gwarray[0] = LLDATA;
		for (k = 1; k < num_gwnums; k++) {
			gwarray[k] = gwalloc ();
			if (gwarray[k] == NULL) {
				num_gwnums = k;
				break;
			}
		}

/* Init data area with a pre-determined value */

restart_test:	units_bit = 0;
		dbltogw (4.0, LLDATA);
		g = LLDATA;

/* Do Lucas-Lehmer iterations */

		for (k = 0; k < ll_iters; k++) {
			int	fd;
			short	type;
			unsigned long trash;

/* Copy previous squared value (so we plow through memory) */

			if (k && num_gwnums > 1) {
				gwnum	prev;
				prev = g;
				g = gwarray[k % num_gwnums];
				gwcopy (prev, g);
			}

/* One Lucas-Lehmer test with error checking */

			gwsetnormroutine (0, 1, 0);
			gwstartnextfft (k != 100 && k != ll_iters - 1);
			lucas_fixup (&units_bit);
			gwsquare (g);

/* If the sum of the output values is an error (such as infinity) */
/* then raise an error. */

			if (gw_test_illegal_sumout ()) {
				OutputBoth (SELFFAIL1);
				SELF_TEST_WARNINGS++;
				if (SELF_TEST_WARNINGS < 100) {
					OutputBoth (SELFFAIL4);
					goto restart_test;
				} else {
					OutputBoth (SELFFAIL6);
					lucasDone ();
					free (gwarray);
					return (FALSE);
				}
			}

/* Check that the sum of the input numbers squared is approximately */
/* equal to the sum of unfft results. */

			if (gw_test_mismatched_sums ()) {
				sprintf (buf, SELFFAIL2,
					 gwsumout (g),
					 gwsuminp (g));
				OutputBoth (buf);
				OutputBoth (SELFFAIL5);
				SELF_TEST_ERRORS++;
				lucasDone ();
				free (gwarray);
				return (FALSE);
			}

/* Make sure round off error is tolerable */

			if (MAXERR > 0.45) {
				sprintf (buf, SELFFAIL3, MAXERR);
				OutputBoth (buf);
				OutputBoth (SELFFAIL5);
				SELF_TEST_ERRORS++;
				lucasDone ();
				free (gwarray);
				return (FALSE);
			}

/* Abort if user demands it */

			if (escapeCheck ()) {
				lucasDone ();
				free (gwarray);
				return (FALSE);
			}

/* Test our ability to read and write files too. */

			if (k != 100) continue;
			if (--countdown) continue;
			if (g != LLDATA) gwcopy (g, LLDATA);
			if (! writeToFile (filename, 0, 0, 0)) {
				OutputBoth (SELFFAILW);
				SELF_TEST_ERRORS++;
				lucasDone ();
				free (gwarray);
				return (FALSE);
			}
			countdown = 9;
			dbltogw (0.0, LLDATA); /* clear memory */
			if (! readFileHeader (filename, &fd, &type, &trash) ||
			    ! readFileData (fd, &trash, &trash)) {
				OutputBoth (SELFFAILR);
				SELF_TEST_ERRORS++;
				lucasDone ();
				free (gwarray);
				return (FALSE);
			}
			_unlink (filename);
		}

/* Compare final 32 bits with the pre-computed array of correct residues */

		if (g != LLDATA) gwcopy (g, LLDATA);
		generateResidue64 (units_bit, &reshi, &reslo);
		lucasDone ();
		free (gwarray);
		if (reshi != test_data[i].reshi) {
			sprintf (buf, SELFFAIL, reshi, test_data[i].reshi);
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
/* An Advanced/Time 9999 corresponds to type 0, Advanced/Time 9998 */
/* corresponds to type 1, etc. */
/* Type 0 executes much like an LL test, error checking and doing a */
/* careful iteration occasionally */
/* Type 1 does roundoff checking every iteration and accumulates */
/* statistics on the round off data. */
/* Type 2 and higher have not been used much and may not work */

int lucas_QA (
	int	type)
{
	FILE	*fd;

/* Set the title, init random generator */

	title ("QA");
	srand ((unsigned) time (NULL));

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
		unsigned long i, word, bit_in_word, maxerrcnt, loops;
		double	maxsumdiff, maxerr, toterr, M, S;
		unsigned long ge_300, ge_325, ge_350, ge_375, ge_400;
		gwnum	t1, t2;
		unsigned int iters_unchecked;

/* Read a line from the file */

		p = 0;
		fscanf (fd, "%lu,%lu,%lu,%lu,%s\n",
			&p, &fftlen, &iters, &units_bit, &res);
		if (p == 0) break;

/* In a type 4 run, we decrement through exponents to find any with */
/* anamolously high average errors.  After selecting a tentative FFT */
/* crossover, we do a type 4 run looking for a higher average error */
/* below the crossover we selected. */

		for (loops = (type != 4 ? 1 : units_bit % 100); loops--; p-=2){

/* Now run Lucas setup */

		lucasSetup (p, fftlen);
		maxsumdiff = 0.0;
		ge_300 = ge_325 = ge_350 = ge_375 = ge_400 = 0;
		maxerr = 0.0; maxerrcnt = 0; toterr = 0.0;
		iters_unchecked = (type > 3) ? 2 : 40;

/* Check for a randomized units bit */

		if (units_bit >= p) {
			unsigned long hi, lo;
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
				       (type == 3 || type == 4) ?
					 (rand () & 1) ? rand () : -rand () :
				       (i == word) ? (1L << bit_in_word) : 0);
		}

/* The thorough, P-1, and ECM tests use more than one number */

		if (type == 2 || type == 3) {
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

			if (type == 0) {		/* Typical LL test */
				gwsetnormroutine (0, (i & 63) == 37, 0);
				gwstartnextfft (i < iters / 2);
				if (i > iters / 2 && (i & 63) == 44)
					careful_iteration (LLDATA, &units_bit);
				else {
					lucas_fixup (&units_bit);
					gwsquare (LLDATA);
				}
			} else if (type == 1 || type == 4) { /* Gather stats */
				gwsetnormroutine (0, 1, 0);
				gwstartnextfft (i < iters / 2);
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

			if (i > iters_unchecked) {
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

/* Maintain range info */

			if (MAXERR >= 0.300) ge_300++;
			if (MAXERR >= 0.325) ge_325++;
			if (MAXERR >= 0.350) ge_350++;
			if (MAXERR >= 0.375) ge_375++;
			if (MAXERR >= 0.400) ge_400++;

/* Maintain maximum error info */

			if (MAXERR > maxerr) maxerr = MAXERR, maxerrcnt = 1;
			else if (MAXERR == maxerr) maxerrcnt++;
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

		if (type == 1 || type == 3 || type == 4) {
			S = sqrt (S / (iters - iters_unchecked - 1));
			toterr /= iters - iters_unchecked;
			sprintf (buf, "avg: %6.6f, stddev: %6.6f, #stdev to 0.5: %6.6f\n",
				 toterr, S, (0.50 - toterr) / S);
			OutputBoth (buf);
		}

/* Compare residue with correct residue from the input file */

		sprintf (buf, "%08X%08X", reshi, reslo);
		if (type <= 2 && stricmp (res, buf)) {
			sprintf (buf, "Warning: Residue mismatch. Expected %s\n", res);
			OutputBoth (buf);
		}

/* Output message */

		sprintf (buf, "Exp/iters: %lu/%lu, res: %08X%08X, maxerr: %6.6f/%lu, %lu/%lu/%lu/%lu/%lu, maxdiff: %9.9f/%9.9f\n",
			 p, iters, reshi, reslo, maxerr, maxerrcnt,
			 ge_300, ge_325, ge_350, ge_375, ge_400,
			 maxsumdiff, MAXDIFF);
		OutputBoth (buf);
		}
	}
	fclose (fd);

	return (TRUE);
}

/* Generate random FFT data for timing the Lucas-Lehmer code */

void generateRandomData (void)
{
	unsigned long i;

/* Fill data space with random values. */

	srand ((unsigned) time (NULL));
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
	unsigned long i, j, saved, save_limit;
	char	buf[80];
	double	time, saved_times[SAVED_LIMIT];
	int	days, hours, minutes;
	unsigned long best_asm_timers[32]= {0};
	unsigned long *asm_timer;

/* Set the process/thread priority */

	SetPriority ();

/* Clear all timers */

	clear_timers ();

/* Look for special values to run QA suites */

	if (p >= 9994 && p <= 9999) {
		lucas_QA (9999 - p);
		return;
	}
	if (p == 9991) {
		ecm_QA ();
		return;
	}
	if (p == 9992) {
		pminus1_QA ();
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

/* Time a single squaring */

		start_timer (0);
		gwsquare (LLDATA);
		end_timer (0);
		timers[1] += timers[0];
		saved_times[saved++] = timers[0];
		timers[0] = 0;

/* Remember the best asm timers (used when I'm optimizing assembly code) */

		asm_timer = (unsigned long *) COPYZERO[0];
		if (asm_timer != NULL)
		for (j = 0; j < 32; j++)
			if (i == 0 || asm_timer[j] < best_asm_timers[j])
				best_asm_timers[j] = asm_timer[j];

/* Output timer squaring times */

		if (saved == save_limit || i == iterations - 1) {
			for (j = 0; j < saved; j++) {
				OutputStr ("p: ");
				OutputNum (p);
				OutputStr (".  Time: ");
				timers[0] = saved_times[j];
				print_timer (0, TIMER_MS | TIMER_NL | TIMER_CLR);
			}
			saved = 0;
		}

/* Abort early if so requested */

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
	sprintf (buf, minutes == 1 ? "%d minute.\n" : "%d minutes.\n", minutes);
	OutputStr (buf);
//for (i = 0; i < 32; i++) {
//sprintf (buf, "timer %d: %d\n", i, best_asm_timers[i]);
//if(best_asm_timers[i]) OutputBoth (buf);}
}

#define BENCH1 "\nYour timings will be written to the results.txt file.\n"
#define BENCH2 "Compare your results to other computers at http://www.mersenne.org/bench.htm\n"
#define BENCH3 "That web page also contains instructions on how your results can be included.\n\n"

/* Time a few iterations of many FFT lengths */

void primeBench (void)
{
	unsigned long num_lengths, i, j, iterations;
	double	best_time;
	char	buf[512];
	int	fft_lengths[] = {384, 448, 512, 640, 768, 896, 1024, 1280, 1536, 1792, 2048, 2560, 3072, 3584, 4096};

/* Set the process/thread priority */

	SetPriority ();

/* Output startup message */

	OutputStr (BENCH1);
	OutputBoth (BENCH2);
	OutputBoth (BENCH3);

/* Output to the results file a full CPU description */

	getCpuDescription (buf, 1);
	writeResults (buf);
	sprintf (buf, "Prime95 version %s, RdtscTiming=%d\n",
		 VERSION, RDTSC_TIMING);
	writeResults (buf);

/* Loop over all FFT lengths */

	num_lengths = sizeof (fft_lengths) / sizeof (int);
	if (!IniGetInt (INI_FILE, "FullBench", 0)) num_lengths -= 4;
	for (i = 0; i < num_lengths; i++) {

/* Initialize for this FFT length.  Compute the number of iterations to */
/* time.  This is based on the fact that it doesn't take too long for */
/* my 1400 MHz P4 to run 10 iterations of a 1792K FFT. */

		lucasSetup (5000000, fft_lengths[i] * 1024);
		iterations = (unsigned long) (10 * 1792 * CPU_SPEED / 1400 / fft_lengths[i]);
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
				OutputStr ("\nExecution halted.\n");
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
	OutputStr ("Benchmark complete.\n");
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
	FACHSW = hsw;
	FACMSW = msw;
	factorSetup (p);
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

	test_bits = factorLimit (p, work_type);
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

/* Loop testing larger and larger factors until we've tested to the */
/* appropriate number of bits.  Advance one bit at a time to minimize wasted */
/* time looking for a second factor after a first factor is found. */

	while (test_bits > bits) {
	    unsigned int end_bits;
	    unsigned long iters, iters_r;
	    int	stopping, saving;

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
		factorSetup (p);

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
/* Set flag if we are saving or stopping */

			stopping = stopCheck ();
			if ((iters_r & 0x7F) == 0 && !stopping) {
				if (!communicateWithServer ()) stopping = 1;
				if (!pauseWhileRunning ()) stopping = 1;
				time (&current_time);
				saving = (current_time-start_time > write_time);
			} else
				saving = 0;

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

			if (stopping || saving) {
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

/* Do next of the 16 passes */

nextpass:	;
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

		sprintf (buf, "M%ld no factor to 2^%d, WZ%d: %08lX\n",
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

char NOFAC[] = "M%ld no factor from 2^%d to 2^%d, WZ%d: %08lX\n";

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
	unsigned short bits;
	short	pass;
	int	stopping, saving;
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
		_read (fd, &p, sizeof (long));
		_read (fd, &pass, sizeof (short));
		_read (fd, &bits, sizeof (short));
		_read (fd, &FACHSW, sizeof (long));
		_read (fd, &FACMSW, sizeof (long));
		_close (fd);
		continuation = TRUE;
	} else
		continuation = FALSE;

/* Init filename */

	tempFileName (filename, 0);

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

/* Loop through all the bit levels */

		if (!continuation) bits = minbits;
		for ( ; bits <= maxbits; bits++) {

/* Determine how much we should factor */

		if (bits < 64) {
			endpthi = 0;
			endptlo = 1L << (bits-32);
		} else {
			endpthi = 1L << (bits-64);
			endptlo = 0;
		}

/* Sixteen passes! two for the 1 or 7 mod 8 factors times two for the */
/* 1 or 2 mod 3 factors times four for the 1, 2, 3, or 4 mod 5 factors. */

		if (!continuation) pass = 0;
		for ( ; pass < 16; pass++) {

/* Setup the factoring program */

		if (!continuation) {
			if (bits <= 32) {
				FACHSW = 0;
				FACMSW = 0;
			} else if (bits <= 64) {
				FACHSW = 0;
				FACMSW = 1L << (bits-33);
			} else {
				FACHSW = 1L << (bits-65);
				FACMSW = 0;
			}
			if (pass > 0 && FACMSW == 0) FACMSW = 1;
		}
		FACPASS = pass;
		factorSetup (p);
		continuation = FALSE;

/* Loop until all factors tested or factor found */

		iters = 0;
		iters_r = 0;
		for ( ; ; ) {
			int	res;

/* Test for completion */

			if (FACHSW > endpthi ||
			    (FACHSW == endpthi && FACMSW >= endptlo))
				break;

/* Factor some more */

			start_timer (0);
			res = factorAndVerify (p);
			end_timer (0);
			if (res != 2) goto bingo;

/* Send queued messages to the server every so often */
/* Set flag if we are saving or stopping */

			stopping = stopCheck ();
			if ((iters_r & 0x7F) == 0 && !stopping) {
				if (!communicateWithServer ()) stopping = 1;
				if (!pauseWhileRunning ()) stopping = 1;
				time (&current_time);
				saving = (current_time-start_time > write_time);
			} else
				saving = 0;

/* Output informative message */

			if (++iters >= ITER_OUTPUT) {
				char	fmt_mask[80];
				double	percent;
				percent = facpct (pass, bits-1, endpthi, endptlo);
				percent = trunc_percent (percent);
				sprintf (fmt_mask, FACMSG, PRECISION);
				sprintf (buf, fmt_mask, p, bits, percent);
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
				percent = facpct (pass, bits-1, endpthi, endptlo);
				percent = trunc_percent (percent);
				sprintf (fmt_mask, FACMSG, PRECISION);
				sprintf (buf, fmt_mask, p, bits, percent);
				strcat (buf, "\n");
				writeResults (buf);
				iters_r = 0;
			}

/* If an escape key was hit, write out the results and return */

			if (stopping || saving) {
				short	four = 4;
				fd = _open (filename, _O_BINARY | _O_WRONLY | _O_TRUNC | _O_CREAT, 0666);
				_write (fd, &four, sizeof (short));
				_write (fd, &p, sizeof (long)); /* dummy */
				_write (fd, &p, sizeof (long));
				_write (fd, &pass, sizeof (short));
				_write (fd, &bits, sizeof (short));
				_write (fd, &FACHSW, sizeof (long));
				_write (fd, &FACMSW, sizeof (long));
				_commit (fd);
				_close (fd);
				if (stopping) {
					factorDone ();
					return (FALSE);
				}
				start_time = current_time;
			}
		}

/* Next pass */

		if (FACMSW != 0xFFFFFFFF) {
			endpthi = FACHSW;
			endptlo = FACMSW+1;
		} else {
			endpthi = FACHSW+1;
			endptlo = 0;
		}
		}
		}

/* Output message if no factor found */

		sprintf (buf, NOFAC, p, minbits-1, maxbits, PORT, SEC4 (p));
		OutputBoth (buf);
		spoolMessage (PRIMENET_RESULT_MESSAGE, buf);
		goto nextp;

/* Format and output the factor found message */

bingo:		makestr (FACHSW, FACMSW, FACLSW, str);
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

/* Factor next prime */

nextp:		;
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
			if (fachi >= 4194304 ||
			    (fachi >= 4096 && !(CPU_FLAGS & CPU_SSE2))) {
				sprintf (buf, "%ld%s factor too big.\n", p, fac);
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
		factorSetup (p);

/* Factor found, is it a match? */

		do {
			if (factor64 () != 2 &&
			    FACHSW == fachi &&
			    FACMSW == facmid &&
			    FACLSW == faclo) {
				sprintf (buf, "%ld%s factored OK.\n", p, fac);
				OutputSomewhere (buf);
				goto nextp;
			}
		} while (FACMSW == facmid);

/* Uh oh. */

bad:		sprintf (buf, "%ld%s factor not found.\n", p, fac);
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
