/*
 * Windoows specific routines used by PRP program
 */ 

/* Common global variables */

HANDLE	WORKER_THREAD = 0;

/* Common routines */

/* Is this Windows 95/98 or Windows NT? */

int isWindows95 (void)
{
	OSVERSIONINFO info;

	info.dwOSVersionInfoSize = sizeof(OSVERSIONINFO);
	if (GetVersionEx (&info)) {
		return (info.dwPlatformId == VER_PLATFORM_WIN32_WINDOWS);
	}
	return (FALSE);
}

/* Is this Windows 2000 or a later Windows NT version? */

int isWindows2000 ()
{
	OSVERSIONINFO info;

	info.dwOSVersionInfoSize = sizeof(OSVERSIONINFO);
	if (GetVersionEx (&info)) {
		return (info.dwPlatformId == VER_PLATFORM_WIN32_NT &&
			info.dwMajorVersion >= 5);
	}
	return (FALSE);
}

/* Set the thread priority correctly */

void SetPriority (void)
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

/* Return the number of CPUs in the system */

unsigned long num_cpus (void)
{
	SYSTEM_INFO sys;

	GetSystemInfo (&sys);
	return (sys.dwNumberOfProcessors);
}

/* Routines to access the high resolution performance counter */

int isHighResTimerAvailable (void)
{
	LARGE_INTEGER large;
	return (QueryPerformanceCounter (&large));
}

double getHighResTimer (void)
{
	LARGE_INTEGER large;

	QueryPerformanceCounter (&large);
	return ((double) large.HighPart * 4294967296.0 +
		(double) large.LowPart);
}

double getHighResTimerFrequency (void)
{
	LARGE_INTEGER large;

	QueryPerformanceFrequency (&large);
	return ((double) large.HighPart * 4294967296.0 +
		(double) large.LowPart);
}
