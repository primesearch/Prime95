/* Copyright 1995-2000 Just For Fun Software, Inc. */
/* Author:  George Woltman */
/* Email: woltman@alum.mit.edu */

/* Include files */

#include "prime.h"
#include <string.h>

/* Routine definitions */

void rangeStatus ();
void options_cpu ();

/* Get line from the user (stdin) */

void get_line (
	char	*buf)
{
	char	c;
	int	i;
	for (i = 0; ; i++) {
		if (_read (0, &c, 1) != 1) break;
		if (c == '\n' || c == 0) break;
		if (i < 80) *buf++ = c;
	}
	*buf++ = 0;
}

/* Get a number from the user */

unsigned long get_number (
	unsigned long dflt)
{
	char	line[80];
	unsigned long i;
	get_line (line);
	if (line[0] == 0) return (dflt);
	return (atol (line));
}

/* Ask a Yes/No question */

void askYN (
	char	*str,
	int	*val)
{
	char	buf[80];
	printf ("%s (%s): ", str, *val ? "Y" : "N");
	get_line (buf);
	if (buf[0] == 0) return;
	*val = (buf[0] == 'Y' || buf[0] == 'y');
}

/* Ask a number question */

void askNum (
	char	*str,
	unsigned long *val,
	unsigned long min,
	unsigned long max)
{
	char	buf[80];
	unsigned long newval;
	printf ("%s (%ld): ", str, *val);
loop:	get_line (buf);
	if (buf[0] == 0) return;
	newval = atol (buf);
	if (min || max) {
		if (newval < min || newval > max) {
			printf ("Please enter a value between %ld and %ld. ",
				min, max);
			goto loop;
		}
	}
	*val = newval;
}

/* Ask a number question */

void askNumNoDflt (
	char	*str,
	unsigned long *val,
	unsigned long min,
	unsigned long max)
{
	char	buf[80];
	unsigned long newval;
	printf ("%s: ", str);
loop:	get_line (buf);
	if (buf[0] == 0) goto loop;
	newval = atol (buf);
	if (min || max) {
		if (newval < min || newval > max) {
			printf ("Please enter a value between %ld and %ld. ",
				min, max);
			goto loop;
		}
	}
	*val = newval;
}

/* Ask a string question */

void askStr (
	char	*str,
	char	*val,
	unsigned long maxlen)
{
	char	buf[80];
	if (val[0])
		printf ("%s (%s): ", str, val);
	else
		printf ("%s: ", str);
loop:	get_line (buf);
	if (buf[0] == 0) return;
	if (strlen (buf) > maxlen) {
		printf ("Maximum string length is %ld characters. ", maxlen);
		goto loop;
	}
	strcpy (val, buf);
}

/* Wait for user input - gives the user time to read the screen */

void askOK ()
{
	char	str[80];
	if (THREAD_KILL) return;
	printf ("\nHit enter to continue: ");
	get_line (str);
}

/* Ask user if he is satisfied with his dialog responses */

int askOkCancel ()
{
	char	buf[80];
	if (THREAD_KILL) return (FALSE);
	printf ("\nAccept the answers above? (Y): ");
	get_line (buf);
	return (buf[0] == 0 || buf[0] == 'Y' || buf[0] == 'y');
}

/* Ask user if he is satisfied with his dialog responses */

int askYesNo (
	char	dflt)
{
	char	buf[80];
	if (THREAD_KILL) return (FALSE);
	printf (" (%c): ", dflt);
	get_line (buf);
	if (buf[0] == 0) buf[0] = dflt;
	return (buf[0] == 'Y' || buf[0] == 'y');
}

/* Output a long string with a max of 75 characters to a line */

void outputLongLine (
	char	*buf)
{
	char	line[80];
	char	*p;
	int	i, j;

	for (p = buf; ; ) {
		for (i = 0; i < 75; i++) {
			line[i] = p[i];
			if (p[i] == 0 || p[i] == '\n') { j = i; break; }
			if (p[i] == ' ' || p[i] == '.' || p[i] == ',') j = i;
		}
		line[j+1] = 0;
		printf ("%s", line);
		if (p[j] == 0) break;
		if (p[j] != '\n') printf ("\n");
		p += j + 1;
		while (*p == ' ') p++;
	}
}

/* Test/Primenet dialog */

#define MSG_BIGONES1 "\nA 500 MHz Pentium-III computer will take a full year to test just\n"
#define MSG_BIGONES2 "one 10,000,000 digit number.  Your chance of finding a new prime\n"
#define MSG_BIGONES3 "is roughly 1 in 250,000.  Before proceeding, read the prize rules\n"
#define MSG_BIGONES4 "at http://www.mersenne.org/prize.htm.  Do you accept these rules\n"
#define MSG_BIGONES5 "and still want to search for 10,000,000 digit primes"

void test_primenet ()
{
	int	m_primenet, m_dialup, m_work_dflt;
	int	m_bigones, m_lucas, m_factor, m_dblchk;
	unsigned long m_work;
	short	work_pref;

	m_primenet = USE_PRIMENET;
	m_dialup = DIAL_UP;
	m_work = DAYS_OF_WORK;
	if (WORK_PREFERENCE == 0) {
		m_work_dflt = 1;
		work_pref = default_work_type ();
	} else {
		m_work_dflt = 0;
		work_pref = WORK_PREFERENCE;
	}
	m_bigones = !! (work_pref & PRIMENET_ASSIGN_BIGONES);
	m_lucas = !! (work_pref & PRIMENET_ASSIGN_TEST);
	m_factor = !! (work_pref & PRIMENET_ASSIGN_FACTOR);
	m_dblchk = !! (work_pref & PRIMENET_ASSIGN_DBLCHK);

	askYN ("Use PrimeNet to get work and report results", &m_primenet);
	if (!m_primenet) goto done;
	askYN ("Use a dial-up connection to the Internet", &m_dialup);
	askNum ("Always have at least this many days of work queued up",
		&m_work, 1, 90);
	askYN ("Request whatever type of work makes the most sense",
		&m_work_dflt);
	if (!m_work_dflt) {
		askYN ("Request 10,000,000 digit numbers to test", &m_bigones);
		askYN ("Request Mersenne numbers to run primality tests",
			&m_lucas);
		askYN ("Request Mersenne numbers to double-check", &m_dblchk);
		askYN ("Request Mersenne numbers to factor", &m_factor);
	}

done:	if (askOkCancel ()) {
		if (!USE_PRIMENET && m_primenet) {
			USE_PRIMENET = 1;
			spoolMessage (PRIMENET_MAINTAIN_USER_INFO, NULL);
			spoolMessage (PRIMENET_SET_COMPUTER_INFO, NULL);
			spoolExistingResultsFile ();
		}
		USE_PRIMENET = m_primenet;
		DIAL_UP = m_dialup;
		DAYS_OF_WORK = m_work;
		if (m_work_dflt)
			WORK_PREFERENCE = 0;
		else {
			work_pref =
				(m_bigones ? PRIMENET_ASSIGN_BIGONES : 0) +
				(m_lucas ? PRIMENET_ASSIGN_TEST : 0) +
				(m_factor ? PRIMENET_ASSIGN_FACTOR : 0) +
				(m_dblchk ? PRIMENET_ASSIGN_DBLCHK : 0);
			if (! (WORK_PREFERENCE & PRIMENET_ASSIGN_BIGONES) &&
			    work_pref & PRIMENET_ASSIGN_BIGONES) {
				int	m_ask;
				m_ask = 0;
				printf (MSG_BIGONES1);
				printf (MSG_BIGONES2);
				printf (MSG_BIGONES3);
				printf (MSG_BIGONES4);
				askYN (MSG_BIGONES5, &m_ask);
				if (m_ask) WORK_PREFERENCE = work_pref;
			} else
				WORK_PREFERENCE = work_pref;
		}

		IniWriteInt (INI_FILE, "UsePrimenet", USE_PRIMENET);
		IniWriteInt (INI_FILE, "DialUp", DIAL_UP);
		IniWriteInt (INI_FILE, "DaysOfWork", DAYS_OF_WORK);
		IniWriteInt (INI_FILE, "WorkPreference", WORK_PREFERENCE);
		CHECK_WORK_QUEUE = 1;
		if (STARTUP_IN_PROGRESS) {
			STARTUP_IN_PROGRESS = 0;
			if (USE_PRIMENET) linuxContinue (NULL);
		}
	} else
		STARTUP_IN_PROGRESS = 0;
}

/* Test/User Information dialog */

void test_user ()
{
	char	m_name[80], m_email[80], m_userid[20];
	char	m_password[20], m_compid[20];
	int	m_sendemail, m_team;

	if (IniGetInt (INI_FILE, "LockUserInfo", 0)) {
		OutputStr ("The user information cannot be changed.\n");
		askOK ();
		return;
	}

	strcpy (m_name, USER_NAME);
	strcpy (m_email, USER_ADDR);
	strcpy (m_userid, USERID);
	strcpy (m_password, USER_PWD);
	strcpy (m_compid, COMPID);
	m_sendemail = NEWSLETTERS;
	m_team = 0;

	askStr ("Your name", m_name, 76);
	askStr ("Your email address", m_email, 76);
	askYN ("Receive occasional newsletters by email from the PrimeNet server", &m_sendemail);
	outputLongLine ("\nYou can either pick your own user ID, password, and computer ID or let the PrimeNet server assign one to you.  See the readme.txt file for details.\n");
	askStr ("Your user ID", m_userid, 14);
	askStr ("Your password", m_password, 8);
	askStr ("Your computer ID", m_compid, 12);
	if (!STARTUP_IN_PROGRESS) {
		printf ("Create a team using the information above.  Team accounts are protected by\n");
		askYN ("the PrimeNet server against unintentional or unauthorized changes", &m_team);
	}

	if (askOkCancel ()) {
		if (strcmp (COMPID, m_compid) != 0) {
			strcpy (COMPID, m_compid);
			sanitizeString (COMPID);
			IniWriteString (LOCALINI_FILE, "ComputerID", COMPID);
			spoolMessage (PRIMENET_SET_COMPUTER_INFO, NULL);
		}
		if (OLD_USERID[0] == 0 &&
		    strcmp (USERID, m_userid) != 0) {
			strcpy (OLD_USERID, USERID);
			strcpy (OLD_USER_PWD, USER_PWD);
			IniWriteString (INI_FILE, "OldUserID", OLD_USERID);
			IniWriteString (INI_FILE, "OldUserPWD", OLD_USER_PWD);
		}
		strcpy (USER_PWD, m_password);
		sanitizeString (USER_PWD);
		IniWriteString (INI_FILE, "UserPWD", USER_PWD);
		if (strcmp (USER_NAME, m_name) != 0 ||
		    strcmp (USER_ADDR, m_email) != 0 ||
		    NEWSLETTERS != m_sendemail ||
		    m_team ||
		    strcmp (USERID, m_userid) != 0) {
			strcpy (USER_NAME, m_name);
			strcpy (USER_ADDR, m_email);
			strcpy (USERID, m_userid);
			sanitizeString (USERID);
			NEWSLETTERS = m_sendemail;
			IniWriteString (INI_FILE, "UserName", USER_NAME);
			IniWriteString (INI_FILE, "UserEmailAddr", USER_ADDR);
			IniWriteInt (INI_FILE, "Newsletters", NEWSLETTERS);
			IniWriteString (INI_FILE, "UserID", USERID);
			spoolMessage (m_team ?
				PRIMENET_MAINTAIN_USER_INFO + 0x80 :
				PRIMENET_MAINTAIN_USER_INFO, NULL);
		}
		if (STARTUP_IN_PROGRESS) options_cpu ();
	} else
		STARTUP_IN_PROGRESS = 0;
}

/* Test/Vacation or Holiday dialog */

void test_vacation ()
{
	unsigned long m_vacation_days;
	int	m_computer_on;

	m_vacation_days = (secondsUntilVacationEnds () + 43200) / 86400;
	m_computer_on = ON_DURING_VACATION;

	outputLongLine ("If you are going on vacation and your computer will not be running or your computer will not be connected to the Internet, then the PrimeNet server must be informed of new expected completion dates.\n");
	outputLongLine ("\nPlease make sure you are connected to the Internet long enough to send the new completion dates.\n");
	askNum ("Days of vacation", &m_vacation_days, 0, 120);
	askYN ("Computer will be on during your absence", &m_computer_on);

	if (askOkCancel ()) {
		if (m_vacation_days) {
			time (&VACATION_END);
			VACATION_END += m_vacation_days * 86400;
		} else
			VACATION_END = 0;
		ON_DURING_VACATION = m_computer_on;
                IniWriteInt (LOCALINI_FILE, "VacationEnd", VACATION_END);
                IniWriteInt (LOCALINI_FILE, "VacationOn", ON_DURING_VACATION);
		if (VACATION_END && !ON_DURING_VACATION)
			IniWriteInt (LOCALINI_FILE, "RollingStartTime", 0);
		UpdateEndDates ();
		MANUAL_COMM |= 0x2;
		CHECK_WORK_QUEUE = 1;
		communicateWithServer ();
	}
}

/* Output a status report for the range */

void rangeStatus ()
{
	char	buf[2000];

	rangeStatusMessage (buf);
	outputLongLine (buf);
	askOK ();
}

/* Advanced/Test dialog */

void advanced_test ()
{
	unsigned long m_p;
#define NOTPRIMEERR "This number is not prime, there is no need to test it.\n"

loop:	m_p = 0;
	askNumNoDflt ("Exponent to test", &m_p, MIN_PRIME, MAX_PRIME);

	if (askOkCancel ()) {
		if (! isPrime (m_p)) {
			printf (NOTPRIMEERR);
			goto loop;
		}
		IniFileOpen (WORKTODO_FILE, 1);
		IniInsertLineAsInt (WORKTODO_FILE, 1, "AdvancedTest", m_p);
		linuxContinue ("\nWork added to worktodo.ini file.  Another mprime is running.\n");
		askOK ();
	}
}

/* Advanced/Time dialog */

void advanced_time ()
{
	unsigned long m_p, m_iter;

	m_p = 10000000;
	m_iter = 10;

	askNum ("Exponent to time", &m_p, MIN_PRIME, MAX_PRIME);
	askNum ("Number of Iterations", &m_iter, 1, 1000);
	if (askOkCancel ()) {
		primeTime (m_p, m_iter);
		askOK ();
	}
}

/* Advanced/P-1 dialog */

void advanced_pminus1 ()
{
	unsigned long m_p, m_bound1, m_bound2;
	int	m_plus1;

	m_p = 0;
	m_bound1 = 1000000;
	m_bound2 = 0;
	m_plus1 = 0;

	askNumNoDflt ("Exponent", &m_p, 100, 20500000);
	askNum ("Bound #1", &m_bound1, 100, 1000000000);
	askNum ("Bound #2", &m_bound2, 0, 4000000000UL);
	askYN ("Factor 2^N+1", &m_plus1);

	if (askOkCancel ()) {
		struct work_unit w;
		w.work_type = WORK_PMINUS1;
		w.p = m_p;
		w.B1 = m_bound1;
		w.B2_start = 0;
		w.B2_end = m_bound2;
		w.plus1 = m_plus1;
		addWorkToDoLine (&w);
		linuxContinue ("\nWork added to worktodo.ini file.  Another mprime is running.\n");
		askOK ();
	}
}

/* Advanced/ECM dialog */

void advanced_ecm ()
{
	unsigned long m_p, m_bound1, m_bound2, m_num_curves;
	int	m_plus1;
	double	m_curve;

	m_p = 0;
	m_bound1 = 1000000;
	m_bound2 = 0;
	m_curve = 0.0;
	m_num_curves = 100;
	m_plus1 = 0;

	askNumNoDflt ("Exponent", &m_p, 100, 20500000);
	askNum ("Bound #1", &m_bound1, 100, 1000000000);
	askNum ("Bound #2", &m_bound2, 0, 4000000000UL);
	askNum ("Curves to test", &m_num_curves, 1, 100000);
	askYN ("Factor 2^N+1", &m_plus1);

	if (askOkCancel ()) {
		struct work_unit w;
		w.work_type = WORK_ECM; 
		w.p = m_p;
		w.B1 = m_bound1;
		w.B2_start = 0;
		w.B2_end = m_bound2;
		w.curves_to_do = m_num_curves;
		w.curves_completed = 0;
		w.curve = m_curve;
		w.plus1 = m_plus1;
		addWorkToDoLine (&w);
		linuxContinue ("\nWork added to worktodo.ini file.  Another mprime is running.\n");
		askOK ();
	}
}

/* Advanced/Priority dialog */

void advanced_priority ()
{
	unsigned long m_priority;

	m_priority = PRIORITY;

	outputLongLine ("Pick a priority between 1 and 10 where 1 is the lowest priority and 10 is the highest.\n");
	outputLongLine ("It is strongly recommended that you use the default priority of 1.  Your throughput will probably not improve by using a higher priority.  The only time you should raise the priority is when another process, such as a screen saver, is stealing CPU cycles from this program.\n");
	askNum ("Priority", &m_priority, 1, 10);

	if (askOkCancel ()) {
		PRIORITY = m_priority;
		IniWriteInt (INI_FILE, "Priority", PRIORITY);
	}
}

/* Advanced/Manual Communication dialog */

void advanced_manualcomm ()
{
	int	m_manual_comm, m_comm_now, m_new_dates;

	m_manual_comm = MANUAL_COMM & 0x1;
	m_comm_now = 1;
	m_new_dates = 0;

	m_manual_comm = !m_manual_comm;
	askYN ("Contact PrimeNet server automatically", &m_manual_comm);
	m_manual_comm = !m_manual_comm;
	askYN ("Contact PrimeNet server now", &m_comm_now);
	askYN ("Send new expected completion dates to server", &m_new_dates);

	if (askOkCancel ()) {
		MANUAL_COMM = m_manual_comm;
		IniWriteInt (INI_FILE, "ManualComm", MANUAL_COMM);
		if (m_new_dates) UpdateEndDates ();
		if (m_comm_now) {
			MANUAL_COMM |= 0x2;
			CHECK_WORK_QUEUE = 1;
			communicateWithServer ();
		}
	}
}

/* Advanced/Quit Gimps dialog */

void advanced_quit ()
{

	if (!USE_PRIMENET) {
		outputLongLine (MANUAL_QUIT);
		if (askYesNo ('N')) {
			writeResults ("Quitting GIMPS.\n");
			IniDeleteAllLines (WORKTODO_FILE);
		}
	} else {
		outputLongLine (PRIMENET_QUIT);
		if (askYesNo ('N')) {
			writeResults ("Quitting GIMPS.\n");
			spoolMessage (999, NULL);
			MANUAL_COMM |= 0x2;
			CHECK_WORK_QUEUE = 1;
			communicateWithServer ();
			askOK ();
		}
	}
}

/* Options/CPU dialog */

void options_cpu ()
{
	unsigned long m_cpu_type, m_speed, m_hours;
	unsigned long m_day_memory, m_night_memory, max_mem;
	char m_start_time[13];
	char m_end_time[13];

again:	m_cpu_type =
		(CPU_TYPE == 10) ? 0 : (CPU_TYPE == 9) ? 1 :
		(CPU_TYPE == 8) ? 2 : (CPU_TYPE == 6) ? 3 :
		(CPU_TYPE == 5) ? 4 : (CPU_TYPE == 4) ? 5 :
		(CPU_TYPE == 11) ? 6 : (CPU_TYPE == 7) ? 7 : 8;
	m_speed = CPU_SPEED;
	m_hours = CPU_HOURS;
	m_day_memory = DAY_MEMORY;
	m_night_memory = NIGHT_MEMORY;
	minutesToStr (DAY_START_TIME, m_start_time);
	minutesToStr (DAY_END_TIME, m_end_time);

	m_cpu_type++;
	/*if (isPentiumPro ()) limit = 9;
	else if (isPentium ()) limit = 5;
	else limit = 2;*/
	printf ("CPU Type, 1=Pentium III, 2=Pentium II, 3=Celeron,\n");
	printf ("          4=Pentium Pro, 5=Pentium, 6=486,\n");
	askNum ("          7=AMD Athlon, 8=AMD K6, 9 = Cyrix", &m_cpu_type, 1, 9);
	m_cpu_type--;
	askNum ("CPU speed in MHz", &m_speed, 25, 1000);
	askNum ("Hours per day this program will run", &m_hours, 1, 24);

	printf ("\nPlease see the readme.txt file for very important\n");
	printf ("information on the available memory settings.\n\n");

	max_mem = physical_memory () - 8;
	if (max_mem < 8) max_mem = 8;
	askNum ("Daytime available memory in MB", &m_day_memory, 8, max_mem);
	askNum ("Nighttime available memory in MB", &m_night_memory, 8, max_mem);
	askStr ("Daytime begins at", (char *) &m_start_time, 12);
	askStr ("Daytime ends at", (char *) &m_end_time, 12);

	if (askOkCancel ()) {
		unsigned int new_cpu_type, new_day_start_time, new_day_end_time;

		new_cpu_type = (m_cpu_type == 0) ? 10 :
			       (m_cpu_type == 1) ? 9 :
			       (m_cpu_type == 2) ? 8 :
			       (m_cpu_type == 3) ? 6 :
			       (m_cpu_type == 4) ? 5 :
			       (m_cpu_type == 5) ? 4 :
			       (m_cpu_type == 6) ? 11 :
			       (m_cpu_type == 7) ? 7 : 3;
		if (CPU_SPEED != m_speed ||
		    CPU_TYPE != new_cpu_type ||
		    CPU_HOURS != m_hours) {
			ROLLING_AVERAGE = 1000;
			IniWriteInt (LOCALINI_FILE, "RollingAverage", 1000);
			IniWriteInt (LOCALINI_FILE, "RollingStartTime", 0);
			spoolMessage (PRIMENET_SET_COMPUTER_INFO, NULL);
			UpdateEndDates ();
		}
		new_day_start_time = strToMinutes ((char *) &m_start_time);
		new_day_end_time = strToMinutes ((char *) &m_end_time);
		if (DAY_MEMORY != m_day_memory ||
		    NIGHT_MEMORY != m_night_memory ||
		    DAY_START_TIME != new_day_start_time ||
		    DAY_END_TIME != new_day_end_time)
			memSettingsChanged ();
		CPU_SPEED = m_speed;
		CPU_TYPE = new_cpu_type;
		CPU_HOURS = m_hours;
		DAY_MEMORY = m_day_memory;
		NIGHT_MEMORY = m_night_memory;
		DAY_START_TIME = new_day_start_time;
		DAY_END_TIME = new_day_end_time;
		IniWriteInt (LOCALINI_FILE, "CPUType", CPU_TYPE);
		IniWriteInt (LOCALINI_FILE, "CPUSpeed", CPU_SPEED);
		IniWriteInt (LOCALINI_FILE, "CPUHours", CPU_HOURS);
		IniWriteInt (LOCALINI_FILE, "DayMemory", DAY_MEMORY);
		IniWriteInt (LOCALINI_FILE, "NightMemory", NIGHT_MEMORY);
		IniWriteInt (LOCALINI_FILE, "DayStartTime", DAY_START_TIME);
		IniWriteInt (LOCALINI_FILE, "DayEndTime", DAY_END_TIME);

		if (!IniGetInt (INI_FILE, "AskedAboutMemory", 0)) {
			IniWriteInt (INI_FILE, "AskedAboutMemory", 1);
			if (DAY_MEMORY == 8 && NIGHT_MEMORY == 8) {
				outputLongLine (MSG_MEMORY);
				if (askYesNo ('Y')) goto again;
			}
		}
		 	
		if (STARTUP_IN_PROGRESS) test_primenet ();
	} else
		STARTUP_IN_PROGRESS = 0;
}

/* Options/Preferences dialog */

void options_preferences ()
{
	unsigned long m_iter, m_r_iter, m_disk_write_time;
	unsigned long m_modem, m_retry, m_end_dates;
	int	m_backup, m_noise;

	m_iter = ITER_OUTPUT;
	m_r_iter = ITER_OUTPUT_RES;
	m_disk_write_time = DISK_WRITE_TIME;
	m_modem = MODEM_RETRY_TIME;
	m_retry = NETWORK_RETRY_TIME;
	m_end_dates = DAYS_BETWEEN_CHECKINS;
	m_backup = TWO_BACKUP_FILES;
	m_noise = !SILENT_VICTORY;

	askNum ("Iterations between screen outputs", &m_iter, 1, 999999999);
	askNum ("Iterations between results file outputs",
		&m_r_iter, 10000, 999999999);
	askNum ("Minutes between disk writes", &m_disk_write_time, 10, 999999);
	if (PRIMENET)
		askNum ("Minutes between modem retries", &m_modem, 1, 300);
	if (PRIMENET)
		askNum ("Minutes between network retries", &m_retry, 1, 300);
	if (PRIMENET)
		askNum ("Days between sending end dates", &m_end_dates, 1, 60);
	askYN ("Create Two Backup Files", &m_backup);
	askYN ("Make noise if new Mersenne prime is found", &m_noise);

	if (askOkCancel ()) {
		ITER_OUTPUT = m_iter;
		ITER_OUTPUT_RES = m_r_iter;
		DISK_WRITE_TIME = m_disk_write_time;
		MODEM_RETRY_TIME = m_modem;
		NETWORK_RETRY_TIME = m_retry;
		DAYS_BETWEEN_CHECKINS = m_end_dates;
		TWO_BACKUP_FILES = m_backup;
		SILENT_VICTORY = !m_noise;
		IniWriteInt (INI_FILE, "OutputIterations", ITER_OUTPUT);
		IniWriteInt (INI_FILE, "ResultsFileIterations", ITER_OUTPUT_RES);
		IniWriteInt (INI_FILE, "DiskWriteTime", DISK_WRITE_TIME);
		IniWriteInt (INI_FILE, "NetworkRetryTime", MODEM_RETRY_TIME);
		IniWriteInt (INI_FILE, "NetworkRetryTime2", NETWORK_RETRY_TIME);
		IniWriteInt (INI_FILE, "DaysBetweenCheckins", DAYS_BETWEEN_CHECKINS);
		IniWriteInt (INI_FILE, "TwoBackupFiles", TWO_BACKUP_FILES);
		IniWriteInt (INI_FILE, "SilentVictory", SILENT_VICTORY);
	}
}

/* Help/About */

void help_about ()
{
	printf ("Mersenne Prime Finder - Version %s.2\n", VERSION);
	printf ("Copyright 1996-2000 Just For Fun Software, Inc.\n");
	printf ("Author: George Woltman\n");
	printf ("Email:  woltman@alum.mit.edu\n");
	askOK ();
}

/* Help/About PrimeNet Server */

void help_about_server ()
{
	struct primenetPingServerInfo pkt;

	memset (&pkt, 0, sizeof (pkt));
	pkt.u.serverInfo.versionNumber = PRIMENET_VERSION;
	if (sendMessage (PRIMENET_PING_SERVER_INFO, &pkt)) {
		printf (PING_ERROR);
	} else {
		printf ("Version: %s\n", pkt.u.serverInfo.buildID);
		printf ("Server Name: %s\n", pkt.u.serverInfo.primenetServerName);
		printf ("Email Addr: %s\n", pkt.u.serverInfo.adminEmailAddr);
	}
	askOK ();
}

/* Display the main menu */

void main_menu ()
{
	unsigned long choice;

mloop:	if (THREAD_KILL) return;
	printf ("\t     Main Menu\n");
loop:	printf ("\n");
	printf ("\t 1.  Test/Primenet\n");
	printf ("\t 2.  Test/User Information\n");
	printf ("\t 3.  Test/Vacation or Holiday\n");
	printf ("\t 4.  Test/Status\n");
	printf ("\t 5.  Test/Continue\n");
	printf ("\t 6.  Test/Exit\n");
	printf ("\n");
	printf ("\t 7.  Advanced/Test\n");
	printf ("\t 8.  Advanced/Time\n");
	printf ("\t 9.  Advanced/P-1\n");
	printf ("\t10.  Advanced/ECM\n");
	printf ("\t11.  Advanced/Priority\n");
	printf ("\t12.  Advanced/Manual Communication\n");
	printf ("\t13.  Advanced/Quit Gimps\n");
	printf ("\n");
	printf ("\t14.  Options/CPU\n");
	printf ("\t15.  Options/Preferences\n");
	printf ("\t16.  Options/Self Test\n");
	printf ("\t17.  Options/Torture Test\n");
	printf ("\n");
	printf ("\t18.  Help/About\n");
	printf ("\t19.  Help/About PrimeNet Server\n");
	printf ("Your choice: ");
	choice = get_number (0);
	if (choice <= 0 || choice >= 20) {
		printf ("\t     Invalid choice\n");
		goto loop;
	}

/* Display the main menu and switch off the users choice */

	printf ("\n");
	switch (choice) {

/* Test/Primenet dialog */

	case 1:
		test_primenet ();
		break;

/* Test/User Information dialog */

	case 2:
		test_user ();
		break;

/* Test/Vacation or Holiday dialog */

	case 3:
		test_vacation ();
		break;

/* Test/Status message */

	case 4:
		rangeStatus ();
		break;

/* Test/Continue */

	case 5:
		linuxContinue ("Another mprime is running.\n");
		askOK ();
		break;

/* Test/Exit */

	case 6:
		return;

/* Advanced/Test dialog */

	case 7:
		advanced_test ();
		break;

/* Advanced/Time dialog */

	case 8:
		advanced_time ();
		break;

/* Advanced/ECM dialog */

	case 9:
		advanced_pminus1 ();
		break;

/* Advanced/ECM dialog */

	case 10:
		advanced_ecm ();
		break;

/* Advanced/Priority dialog */

	case 11:
		advanced_priority ();
		break;

/* Advanced/Manual Communication dialog */

	case 12:
		advanced_manualcomm ();
		break;

/* Advanced/Quit Gimps dialog */

	case 13:
		advanced_quit ();
		break;

/* Options/CPU dialog */

	case 14:
		options_cpu ();
		break;

/* Options/Preferences dialog */

	case 15:
		options_preferences ();
		break;

/* Options/Self Test */

	case 16:
		selfTest (0);
		askOK ();
		break;

/* Options/Torture Test */

	case 17:
		selfTest (1);
		askOK ();
		break;

/* Help/About */

	case 18:
		help_about ();
		break;

/* Help/About PrimeNet Server */

	case 19:
		help_about_server ();
		break;
	}
	goto mloop;
}
