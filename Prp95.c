/*
 * Common routines and variables used by Prime95, Saver95, and NTPrime
 *
 * Comm95a contains information used only during setup
 * Comm95b contains information used only during execution
 * Comm95c contains information used during setup and execution
 */ 

/* Common global variables */

HANDLE	WORKER_THREAD = 0;

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
		(PRIORITY > 6) ? NORMAL_PRIORITY_CLASS : IDLE_PRIORITY_CLASS);
	WORKER_THREAD = GetCurrentThread ();
	SetThreadPriority (WORKER_THREAD,
		(PRIORITY == 1) ? THREAD_PRIORITY_IDLE :
		(PRIORITY == 2 || PRIORITY == 7) ? THREAD_PRIORITY_LOWEST :
		(PRIORITY == 3 || PRIORITY == 8) ? THREAD_PRIORITY_BELOW_NORMAL :
		(PRIORITY == 4 || PRIORITY == 9) ? THREAD_PRIORITY_NORMAL :
		(PRIORITY == 5 || PRIORITY == 10) ? THREAD_PRIORITY_ABOVE_NORMAL :
		THREAD_PRIORITY_HIGHEST);
	if (CPU_AFFINITY != 99 && !isWindows95 ())
		SetThreadAffinityMask (WORKER_THREAD, 1 << CPU_AFFINITY);
}

/* Call routines provided by Intel to guess the cpu type and speed */

void guessCpuType ()
{
	WORD cpuid;
	struct FREQ_INFO x;

	cpuid = wincpuid ();
	if (cpuid & CLONE_MASK) {
		CPU_TYPE = 3;
		CPU_SPEED = 166;
	} else {
		CPU_TYPE = cpuid;
		x = cpuspeed (0);
		CPU_SPEED = x.norm_freq;
		if (CPU_SPEED == 0) CPU_SPEED = 100;
	}
}

/* Return the number of CPUs in the system */

unsigned long num_cpus ()
{
	SYSTEM_INFO sys;

	GetSystemInfo (&sys);
	return (sys.dwNumberOfProcessors);
}
