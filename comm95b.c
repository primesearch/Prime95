/*
 * Common routines and variables used by Prime95, Saver95, and NTPrime
 *
 * Comm95a contains information used only during setup
 * Comm95b contains information used only during execution
 * Comm95c contains information used during setup and execution
 */ 

/* Common global variables */

HANDLE	CURRENT_THREAD = 0;

/* Common routines */

/* Is this Windows 95/98 or Windows NT? */

int isWindows95 ()
{
	OSVERSIONINFO info;

	info.dwOSVersionInfoSize = sizeof(OSVERSIONINFO);
	if (GetVersionEx (&info)) {
		return (info.dwPlatformId == VER_PLATFORM_WIN32_WINDOWS);
	}
	return (FALSE);
}

/* Set the thread priority correctly */

void SetPriority ()
{
	SetPriorityClass (GetCurrentProcess (),
		(PRIORITY < 2 || PRIORITY > 6) ?
			NORMAL_PRIORITY_CLASS :
			IDLE_PRIORITY_CLASS);
	CURRENT_THREAD = GetCurrentThread ();
	SetThreadPriority (CURRENT_THREAD,
		(PRIORITY == 1) ? THREAD_PRIORITY_IDLE :
		(PRIORITY == 2 || PRIORITY == 7) ? THREAD_PRIORITY_LOWEST :
		(PRIORITY == 3 || PRIORITY == 8) ? THREAD_PRIORITY_BELOW_NORMAL :
		(PRIORITY == 4 || PRIORITY == 9) ? THREAD_PRIORITY_NORMAL :
		(PRIORITY == 5 || PRIORITY == 10) ? THREAD_PRIORITY_ABOVE_NORMAL :
		THREAD_PRIORITY_HIGHEST);
	if (CPU_AFFINITY != 99 && !isWindows95 ())
		SetThreadAffinityMask (CURRENT_THREAD, 1 << CPU_AFFINITY);
}
