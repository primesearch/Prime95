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

int STARTUP_IN_PROGRESS = 0;

/* Routine to eliminate odd puctuation characters from user ID */
/* and computer ID */

void sanitizeString (
	char	*p)
{
	int	i;
	for (i = strlen (p); i > 0 && isspace (p[i-1]); i--) p[i-1] = 0;
	while (*p) {
		if (!IsCharAlphaNumeric (*p) &&
		    *p != '.' && *p != '-' && *p != '_')
			*p = '_';
		p++;
	}
}

/* Create a status report message from the work-to-do file */

#define STAT0 "Below is a report on the work you have queued and any expected completion dates.\n"
#define STAT1 "The chance that one of the %d exponents you are testing will yield a Mersenne prime is about 1 in %ld. "
#define STAT1a "The chance that the exponent you are testing will yield a Mersenne prime is about 1 in %ld. "
#define STAT3 "You have no work queued up."

void rangeStatusMessage (
	char	*buf)		/* 2000 character buffer */
{
	unsigned int i, ll_cnt;
	double	prob, est;
	int	est_is_accurate;

	ll_cnt = 0;
	prob = 0.0;
	est = 0.0; est_is_accurate = TRUE;
	strcpy (buf, STAT0);
	if (!WELL_BEHAVED_WORK) IniFileOpen (WORKTODO_FILE, 0);
	for (i = 1; ; i++) {
		struct work_unit w;
		unsigned long iteration;
		double	pct, work;
		char	timebuf[30];
		char	*work_type;
		unsigned int len;

/* Stop processing worktodo file if buffer is full */

		len = strlen (buf);
		if (len >= 1800) break;

/* Read the line of the work file */

		if (!parseWorkToDoLine (i, &w)) break;

/* If testing then adjust our probabilities */
/* This assumes our error rate is roughly 1.8% */

		if (w.bits < 32) w.bits = 32;
		if (w.work_type == WORK_TEST) {
			ll_cnt++;
			if (w.pminus1ed)
				prob += (double) ((w.bits - 1) * 1.803) / w.p;
			else
				prob += (double) ((w.bits - 1) * 1.733) / w.p;
		}
		if (w.work_type == WORK_DBLCHK) {
			ll_cnt++;
			if (w.pminus1ed)
				prob += (double) ((w.bits - 1) * 1.803 * ERROR_RATE) / w.p;
			else
				prob += (double) ((w.bits - 1) * 1.733 * ERROR_RATE) / w.p;
		}

/* Adjust our time estimate */

		work = work_estimate (&w);
		pct = pct_complete (w.work_type, w.p, &iteration);
		if (pct == 999.0) est_is_accurate = FALSE;
		work *= (1.0 - pct);
		est += work;

/* Add the exponent to the output message */

		buf[len] = 'M';

		if (w.work_type == WORK_ECM)
			work_type = "ECM";
		if (w.work_type == WORK_PMINUS1)
			work_type = "P-1";
		if (w.work_type == WORK_TEST ||
		    w.work_type == WORK_ADVANCEDTEST)
			work_type = "Lucas-Lehmer test";
		if (w.work_type == WORK_FACTOR)
			work_type = "Factoring";
		if (w.work_type == WORK_PFACTOR)
			work_type = "P-1 factoring";
		if (w.work_type == WORK_DBLCHK)
			work_type = "Double-checking";

		if (w.work_type == WORK_ADVANCEDFACTOR) {
			sprintf (buf+len, "Advanced/Factor, %d to %d\n",
				 w.bits, w.B1);
			est_is_accurate = FALSE;
			continue;
		}

		if (w.work_type == WORK_ECM) {
			sprintf (timebuf, "%d curves with B1=%lu\n",
				 w.curves_to_do, w.B1);
			if (w.plus1) buf[len] = 'P';
			est_is_accurate = FALSE;
		} else if (w.work_type == WORK_PMINUS1) {
			sprintf (timebuf, "B1=%lu\n", w.B1);
			if (w.plus1) buf[len] = 'P';
			est_is_accurate = FALSE;
		} else if (est_is_accurate) {
			time_t	this_time;
			time (&this_time);
			if (est + this_time < 2147483640) {
				this_time += (long) est;
				strcpy (timebuf, ctime (&this_time));
				strcpy (timebuf+16, timebuf+19);
			} else
				strcpy (timebuf, "after Jan 1 2038\n");
		} else
			strcpy (timebuf, "unknown completion date\n");

		sprintf (buf+len+1, "%ld, %s, %s", w.p, work_type, timebuf);
	}
	if (!WELL_BEHAVED_WORK) IniFileClose (WORKTODO_FILE);

/* Format more of the message */

	if (est == 0.0) {
		strcpy (buf, STAT3);
		return;
	}

/* Make a pretty message */

	if (ll_cnt == 1)
		sprintf (buf+strlen(buf), STAT1a, (long) (1.0 / prob));
	if (ll_cnt > 1)
		sprintf (buf+strlen(buf), STAT1, ll_cnt, (long) (1.0 / prob));
}
