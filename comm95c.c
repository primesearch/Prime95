/*
 * Common routines and variables used by Prime95 and NTPrime
 *
 * Comm95a contains information used only during setup
 * Comm95b contains information used only during execution
 * Comm95c contains information used during setup and execution
 */ 

#include <ras.h>
#include <winsock.h>
#include <wininet.h>

int	SOCKETS_INITIALIZED = 0;

/* Common routines */

/* Load the PrimeNet DLL, make sure an internet connection is active */

int LoadPrimeNet (void)
{
static	int	RAS_NOT_AVAILABLE = 0;
static	HMODULE	HRAS = 0;
static	DWORD (APIENTRY *RAS_ENUM)(LPRASCONNA, LPDWORD, LPDWORD);
static	DWORD (APIENTRY *RAS_STAT)(HRASCONN, LPRASCONNSTATUSA);
	RASCONN connections[10];
	DWORD	bufsize;
	DWORD	i, num_connections;
	DWORD	ret;

/* Special handling prior to first primenet call. */
/* Init Winsock, requesting version 1.1 */

	if (! SOCKETS_INITIALIZED) {
		static WSADATA zz;
		int	res;
		res = WSAStartup (MAKEWORD (1, 1), &zz);
		if (res != 0) {
			char buf[80];
			sprintf (buf, "ERROR: Winsock initialization returned %d.\n", res);
			OutputStr (buf);
			return (FALSE);
		}
		SOCKETS_INITIALIZED = 1;
	}

/* If we're not using a dial-up connection, let primenet try */
/* to contact the server. */

	if (!DIAL_UP) return (TRUE);

/* Since Windows 95 can bring up a "Connect To" dialog box */
/* on any call to primenet, we try to make sure we are */
/* already connected before we call primenet.  Otherwise, if */
/* no one is at the computer to respond to the "Connect To" */
/* dialog, the thread hangs until some one does respond. */

/* RAS calls, see below, is no longer the MS-prefered method of detecting */
/* an Internet connection.  Starting in version 22.10 we offer a way for */
/* for users to use the prefered wininet.dll method. */
/* InternetGetConnectedState should return FALSE if the modem is not */
/* connected to the Internet. */

	if (IniGetInt (INI_FILE, "AlternateModemDetection", 0)) {
		DWORD	flags;
		if (InternetGetConnectedState (&flags, 0)) return (TRUE);
		goto no_modem_connection;
	}

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

no_modem_connection:
	OutputStr ("Dial-up connection not active.\n");
	return (FALSE);
}

/* Unload the PrimeNet DLL */

void UnloadPrimeNet (void)
{

/* Tell winsock we are done. */

	if (SOCKETS_INITIALIZED) {
		// Should we call WSACancelBlockingCall first??
		// Should we check error code from WSACleanup?
		// Should we call WSACleanup after each communication session
		// with the server?  That is, are we tying up any resources?
		WSACleanup ();
		SOCKETS_INITIALIZED = 0;
	}
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

/* Return the number of MB of physical memory */

unsigned long physical_memory (void)
{
	MEMORYSTATUS mem;

	GlobalMemoryStatus (&mem);
	return (mem.dwTotalPhys >> 20);
}

/* Return the number of CPUs in the system */

unsigned long num_cpus (void)
{
	SYSTEM_INFO sys;

	GetSystemInfo (&sys);
	return (sys.dwNumberOfProcessors);
}

/* Return 1 to print time in AM/PM format.  Return 2 to print */
/* times using a 24-hour clock. */

int getDefaultTimeFormat (void)
{
	char	buf[10];

	GetLocaleInfo (
		LOCALE_USER_DEFAULT, LOCALE_ITIME, (LPTSTR) buf, sizeof (buf));
	return (buf[0] == '0' ? 1 : 2);
}
