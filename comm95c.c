/*
 * Common routines and variables used by Prime95, Saver95, and NTPrime
 *
 * Comm95a contains information used only during setup
 * Comm95b contains information used only during execution
 * Comm95c contains information used during setup and execution
 */ 

#include <ras.h>

/* Global variables */

HMODULE HLIB;
int (__stdcall *PRIMENET)(short, void *);


/* Common routines */

/* Load the PrimeNet DLL, make sure an internet connection is active */

int LoadPrimeNet ()
{
static	int	RAS_NOT_AVAILABLE = 0;
static	HMODULE	HRAS = 0;
static	DWORD (APIENTRY *RAS_ENUM)(LPRASCONNA, LPDWORD, LPDWORD);
static	DWORD (APIENTRY *RAS_STAT)(HRASCONN, LPRASCONNSTATUSA);
	RASCONN connections[10];
	DWORD	bufsize;
	DWORD	i, num_connections;
	DWORD	ret;

/* Load the DLL if necessary. */

	if (!HLIB) {
		char	*DLLname;
		DLLname = USE_HTTP ? "HTTPNET.DLL" : "RPCNET.DLL";
		HLIB = LoadLibrary (DLLname);
		if (HLIB) {
			PRIMENET = (int (__stdcall *)(short, void *))
				GetProcAddress (HLIB, "PrimeNet");
			if (PRIMENET == NULL) {
				FreeLibrary (HLIB);
				HLIB = 0;
			}
		}
		if (!HLIB) {
			char	buf[80];
			sprintf (buf, "Unable to load %s.\n", DLLname);
			OutputStr (buf);
			return (FALSE);
		}
	}

/* Special handling prior to first primenet call. */
/* Since Windows 95 can bring up a "Connect To" dialog box */
/* on any call to primenet, we try to make sure we are */
/* already connected before we call primenet.  Otherwise, if */
/* no one is at the computer to respond to the "Connect To" */
/* dialog, the thread hangs until some one does respond. */

// If we're not using a dial-up connection, let primenet try
// to contact the server.

	if (!DIAL_UP) return (TRUE);

// Unfortunately, the RASAPI32.DLL is not installed on every
// system.  We must load it dynamically.  If the RAS library
// is not found, let primenet.dll try to contact the server.

	if (RAS_NOT_AVAILABLE) return (TRUE);
	if (HRAS == 0) {
		RAS_NOT_AVAILABLE = 1;
		HRAS = LoadLibrary ("rasapi32.dll");
		if (HRAS == 0) return (TRUE);
		RAS_ENUM = (DWORD (APIENTRY *)(LPRASCONNA, LPDWORD, LPDWORD))
			GetProcAddress (HRAS, "RasEnumConnectionsA");
		if (RAS_ENUM == NULL) return (TRUE);
		RAS_STAT = (DWORD (APIENTRY *)(HRASCONN, LPRASCONNSTATUSA))
			GetProcAddress (HRAS, "RasGetConnectStatusA");
		if (RAS_STAT == NULL) return (TRUE);
		RAS_NOT_AVAILABLE = 0;
	}

// Call RAS to see if there are any active connections to the Internet

	connections[0].dwSize = sizeof (RASCONN);
	bufsize = sizeof (connections);
        ret = (*RAS_ENUM) ((RASCONN *) &connections, &bufsize, &num_connections);

// If RAS returns an error who knows what went wrong. 
// Let primenet try to connect anyway.

	if (ret) return (TRUE);

// See if any of these connections are really connected

	for (i = 0; i < num_connections; i++) {
		RASCONNSTATUS status;
		status.dwSize = sizeof (RASCONNSTATUS);
		ret = (*RAS_STAT) (connections[i].hrasconn, &status);
		if (ret) continue;
		if (status.rasconnstate == RASCS_Connected) return (TRUE);
	}

// Print error message if no there are no connections

	OutputStr ("Dial-up connection not active.\n");
	return (FALSE);
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

/* Return the number of MB of physical memory */

unsigned long physical_memory ()
{
	MEMORYSTATUS mem;

	GlobalMemoryStatus (&mem);
	return (mem.dwTotalPhys >> 20);
}

/* Return the number of CPUs in the system */

unsigned long num_cpus ()
{
	SYSTEM_INFO sys;

	GetSystemInfo (&sys);
	return (sys.dwNumberOfProcessors);
}

/* Return 1 to print time in AM/PM format.  Return 2 to print */
/* times using a 24-hour clock. */

int getDefaultTimeFormat ()
{
	char	buf[10];

	GetLocaleInfo (
		LOCALE_USER_DEFAULT, LOCALE_ITIME, (LPTSTR) buf, sizeof (buf));
	return (buf[0] == '0' ? 1 : 2);
}
