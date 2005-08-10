// Prime95.cpp : Defines the class behaviors for the application.
//

#include "stdafx.h"
#include "Prime95.h"

#include "MainFrm.h"
#include "Prime95Doc.h"
#include "Prime95View.h"

#include <aclapi.h>
#include <direct.h>
#include <fcntl.h>
#include <io.h>
#include <memory.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#ifdef _DEBUG
#define new DEBUG_NEW
#undef THIS_FILE
static char THIS_FILE[] = __FILE__;
#endif

// We need a global mutex handle (could be a member of WinApp)
HANDLE	g_hMutexInst = NULL;
LONG	g_MutexNum = 0;
BOOL CALLBACK MyEnumProc (
	HWND	hwnd,	// handle to parent window
	LPARAM	lParam)	// application-defined value
{

// Return true to continue enumerating Windows

	if (GetWindowLongPtr (hwnd, GWLP_USERDATA) != g_MutexNum)
		return (TRUE);

// We've found the instance, return the window handle

	* (HWND *) lParam = hwnd;
	return (FALSE);
}

/////////////////////////////////////////////////////////////////////////////
// CPrime95App

BEGIN_MESSAGE_MAP(CPrime95App, CWinApp)
	//{{AFX_MSG_MAP(CPrime95App)
	ON_COMMAND(ID_APP_ABOUT, OnAppAbout)
		// NOTE - the ClassWizard will add and remove mapping macros here.
		//    DO NOT EDIT what you see in these blocks of generated code!
	//}}AFX_MSG_MAP
	// Standard file based document commands
	ON_COMMAND(ID_FILE_NEW, CWinApp::OnFileNew)
	ON_COMMAND(ID_FILE_OPEN, CWinApp::OnFileOpen)
END_MESSAGE_MAP()

/////////////////////////////////////////////////////////////////////////////
// CPrime95App construction

CPrime95App::CPrime95App()
{
	EnableHtmlHelp ();

	// TODO: add construction code here,
	// Place all significant initialization in InitInstance
}

/////////////////////////////////////////////////////////////////////////////
// The one and only CPrime95App object

CPrime95App theApp;

/////////////////////////////////////////////////////////////////////////////
// CPrime95App initialization

BOOL CPrime95App::InitInstance()
{
	int	orig_cmdShow;
	int	named_ini_files = -1;
	int	torture_test = 0;
	char	*p;

	// Standard initialization
	// If you are not using these features and wish to reduce the size
	//  of your final executable, you should remove from the following
	//  the specific initialization routines you do not need.

	//LoadStdProfileSettings(0);  // Load standard INI file options (including MRU)

	// Register the application's document templates.  Document templates
	//  serve as the connection between documents, frame windows and views.

	CSingleDocTemplate* pDocTemplate;
	pDocTemplate = new CSingleDocTemplate(
		IDR_MAINFRAME,
		RUNTIME_CLASS(CPrime95Doc),
		RUNTIME_CLASS(CMainFrame),       // main SDI frame window
		RUNTIME_CLASS(CPrime95View));
	AddDocTemplate(pDocTemplate);

	// Parse command line for standard shell commands, DDE, file open
	CCommandLineInfo cmdInfo;
	ParseCommandLine(cmdInfo);

	// Dispatch commands specified on the command line
	orig_cmdShow = m_nCmdShow;
	if (!ProcessShellCommand(cmdInfo))
		return FALSE;

/* Change the working directory to the same directory that */
/* the executable is located.  This is especially important */
/* for running prime95 as a Windows 95 service */

	{
		char	buf[256];
		GetModuleFileName (NULL, buf, sizeof (buf));
		strrchr (buf, '\\')[1] = 0;
		_chdir (buf);
	}

/* Initialize gwnum call back routines.  Using callback routines lets the */
/* gwnum library have a nice clean interface for users that do not need */
/* additional functionality that only prime95 uses. */

	StopCheckRoutine = stopCheck;
	OutputBothRoutine = OutputBoth;

/* NT services are not passed command line arguments.  In this case we */
/* encode the -An information in the NT service name. */

	if (NTSERVICENAME[0] && NTSERVICENAME[15] == '-')
		named_ini_files = atoi (&NTSERVICENAME[16]);

// Process command line switches

	for (p = m_lpCmdLine; *p == '//' || *p == '-'; ) {
		p++;
		switch (*p++) {

// Accept a -A switch indicating an alternate set of INI files
// are to be used.

		case 'A':
		case 'a':
			named_ini_files = 0;
			while (isspace (*p)) p++;
			while (isdigit (*p)) {
				named_ini_files = named_ini_files * 10 + (*p - '0');
				p++;
			}
			break;

// Accept a -T switch to run the torture test.

		case 'T':
		case 't':
			torture_test = 1;
			break;

// Accept a -W switch indicating an alternate working directory.

		case 'W':
		case 'w':
			{
			char	buf[256];
			char	*bufp = buf;
			while (isspace (*p)) p++;
			while (*p && !isspace (*p)) *bufp++ = *p++;
			*bufp = 0;
			_chdir (buf);
			}
			break;
		}

// Skip whitespace between switches

		while (isspace (*p)) p++;
	}

// Make sure only one copy of prime95 is running at a time.
// This code is courtesy of Jeroen C. van Gelderen
// I enhanced it to allow multiple copies if they are running
// from different directories or use different -A switches.

	{
		char	buf[256];
		char	*p;
		DWORD mutex_error_code;
		PSID pEveryoneSID = NULL, pAdminSID = NULL;
		PACL pACL = NULL;
		PSECURITY_DESCRIPTOR pSD = NULL;
		EXPLICIT_ACCESS ea[1];
		SID_IDENTIFIER_AUTHORITY SIDAuthWorld = SECURITY_WORLD_SID_AUTHORITY;
		SECURITY_ATTRIBUTES sa;

// Turn directory name into a (likely) unique integer
// Add in the -A value.  Use this integer to create a mutex name.

		_getcwd (buf, 255);
		for (p = buf; *p; p++)
			g_MutexNum = g_MutexNum * 17 + *p;
		g_MutexNum += named_ini_files;
		sprintf (buf, "Global\\GIMPS%ld", g_MutexNum);

/* Create a world access security descriptor to share the Mutex we */
/* are about to create.  If we run into any troubles, assume this is */
/* a Windows 95/98/Me system and create a simple Mutex. */
		
// Create a well-known SID for the Everyone group.

		if (! AllocateAndInitializeSid (
				&SIDAuthWorld, 1, SECURITY_WORLD_RID,
				0, 0, 0, 0, 0, 0, 0, &pEveryoneSID))
			goto simple_mutex;

// Initialize an EXPLICIT_ACCESS structure for an ACE.
// The ACE will allow the Administrators group full access to the key.

		ZeroMemory (&ea, sizeof (EXPLICIT_ACCESS));
		ea[0].grfAccessPermissions = EVENT_ALL_ACCESS;
		ea[0].grfAccessMode = SET_ACCESS;
		ea[0].grfInheritance= NO_INHERITANCE;
		ea[0].Trustee.TrusteeForm = TRUSTEE_IS_SID;
		ea[0].Trustee.TrusteeType = TRUSTEE_IS_WELL_KNOWN_GROUP;
		ea[0].Trustee.ptstrName  = (LPTSTR) pEveryoneSID;

// Create a new ACL that contains the new ACEs.

		if (SetEntriesInAcl (1, ea, NULL, &pACL) != ERROR_SUCCESS)
			goto simple_mutex;

// Initialize a security descriptor.

		pSD = (PSECURITY_DESCRIPTOR)
			LocalAlloc (LPTR, SECURITY_DESCRIPTOR_MIN_LENGTH);
		if (pSD == NULL) goto simple_mutex;
		if (! InitializeSecurityDescriptor (pSD, SECURITY_DESCRIPTOR_REVISION))
			goto simple_mutex;

// Add the ACL to the security descriptor.

		if (! SetSecurityDescriptorDacl (pSD, TRUE, pACL, FALSE))
			goto simple_mutex;

// Initialize a security attributes structure.

		sa.nLength = sizeof (SECURITY_ATTRIBUTES);
		sa.lpSecurityDescriptor = pSD;
		sa.bInheritHandle = FALSE;

// Create our mutex.  Windows XP uses terminal services to support user
// switching.  The "Global\" prefix is required so that this mutex is
// is created in the global kernel objects namespace.  Unfortunately,
// the "\" character raises an error on Windows 95/98/Me systems.

		g_hMutexInst = CreateMutex (
			&sa,   // World access
			FALSE, // Not owned !!
			buf);  // Unique name
		if (g_hMutexInst == NULL)
simple_mutex:	 	g_hMutexInst = CreateMutex (
				NULL,  // No security stuff
				FALSE, // Not owned !!
				buf+7);  // Unique name
		mutex_error_code = GetLastError ();

// Cleanup all the security structures we initialized

		if (pEveryoneSID) FreeSid (pEveryoneSID);
		if (pACL) LocalFree (pACL);
	        if (pSD) LocalFree (pSD);

// Test for failure

		if (g_hMutexInst == NULL)
			return 0;

// If mutex already exists then another instance is already running

		if (mutex_error_code == ERROR_ALREADY_EXISTS) {
			HWND	hwndPrevInst = 0;

// Give other instance a little time to display it's main window

			Sleep (750);

// Find the window handle

			EnumWindows (&MyEnumProc, (LPARAM) &hwndPrevInst);

// Unhide the other instance's window

			if (hwndPrevInst) {
				ShowWindow (hwndPrevInst, SW_HIDE);
				ShowWindow (hwndPrevInst, SW_SHOWMINIMIZED);
				ShowWindow (hwndPrevInst, SW_SHOWNORMAL);
			}
			CloseHandle (g_hMutexInst);
			return 0;
		}

// Set the window user data so we can be identified by
// another instance of this program.

		SetWindowLongPtr (m_pMainWnd->m_hWnd, GWLP_USERDATA, (LONG_PTR) g_MutexNum);
	}

/* Read the INI files. */

	nameIniFiles (named_ini_files);
	readIniFiles ();

/* Before processing the rest of the INI file, hide and/or */
/* position the main window */

	m_pMainWnd->SetWindowText ("Prime95");

	WINDOWPLACEMENT wp;
	m_pMainWnd->GetWindowPlacement (&wp);
	int left = IniGetInt (INI_FILE, "Left", 0);
	int top = IniGetInt (INI_FILE, "Top", 0);
	int right = IniGetInt (INI_FILE, "Right", 0);
	int bottom = IniGetInt (INI_FILE, "Bottom", 0);
	if (right + left + top + bottom != 0) {
		wp.rcNormalPosition.left = left;
		wp.rcNormalPosition.top = top;
		wp.rcNormalPosition.right = right;
		wp.rcNormalPosition.bottom = bottom;
	}
	wp.showCmd = HIDE_ICON ? SW_HIDE : SW_SHOWMINIMIZED;
	m_pMainWnd->SetWindowPlacement (&wp);

	if (TRAY_ICON) TrayMessage (NIM_ADD, "Prime95", 0);

	// See if we are running as a Windows95 service
	WINDOWS95_SERVICE = IniGetInt (INI_FILE, "Windows95Service", 0);
	WINDOWS95_A_SWITCH = named_ini_files;
	Service95 ();

	// Set flag to read spool file.  We must see if there
	// are messages queued up for the server.
	if (USE_PRIMENET) SPOOL_FILE_CHANGED = 1;

	// To work around a bug where the "Connect To" dialog
	// box comes up on the first attempt to contact the 
	// server, bring that dialog box up now.  We found a way
	// to work around this problem, but maybe it will be
	// useful to do this in other cases.
	if (USE_PRIMENET && IniGetInt (INI_FILE, "PingAtStartup", 0)) {
		struct primenetPingServerInfo pkt;
		memset (&pkt, 0, sizeof (pkt));
		pkt.u.serverInfo.versionNumber = PRIMENET_VERSION;
		sendMessage (PRIMENET_PING_SERVER_INFO, &pkt);
	}

	// Run the torture test if asked to
	if (torture_test) {
		m_pMainWnd->ShowWindow (orig_cmdShow);
		m_pMainWnd->PostMessage (WM_COMMAND, USR_TORTURE, 0);
	}

	// On first run, get user name and email address
	// before contacting server for a work assignment
	else if (USE_PRIMENET && USERID[0] == 0 &&
		 !IniGetInt (INI_FILE, "StressTester", 0)) {
		m_pMainWnd->ShowWindow (orig_cmdShow);
		STARTUP_IN_PROGRESS = 1;
		m_pMainWnd->PostMessage (WM_COMMAND, USR_WELCOME, 0);
	}

	// Auto-continue if there is any work to do.
	else if (USE_PRIMENET || WELL_BEHAVED_WORK || IniGetNumLines (WORKTODO_FILE)) {
		m_pMainWnd->PostMessage (WM_COMMAND, IDM_CONTINUE, 0);
	}

	// Otherwise, show the window
	else if (!HIDE_ICON) {
		m_pMainWnd->ShowWindow (orig_cmdShow);
		ChangeIcon (IDLE_ICON);
	}

	// If a broadcast message from the server has been received but
	// never viewed by the user, then try to display it now.
	BroadcastMessage (NULL);

	// Initialization complete
	return TRUE;
}

void CALLBACK EXPORT TimerCallback (
	HWND	hWnd,		//handle of CWnd that called SetTimer
	UINT	nMsg,		//WM_TIMER
	UINT_PTR nIDEvent,	//timer identification
	DWORD	dwTime)		//system time
{
	BlinkIcon (-2);
}

void ChangeIcon (
	int	icon_id)
{
	CPrime95App *app;
static	UINT	icon = IDR_MAINFRAME;

	if (EXIT_IN_PROGRESS) return;

	app = (CPrime95App *) AfxGetApp();
	if (icon_id >= 0)
		icon = (icon_id == WORKING_ICON) ? IDR_MAINFRAME : IDI_YELLOW_ICON;
	app->m_pMainWnd->SetIcon (app->LoadIcon (icon), 1);
	if (TRAY_ICON) app->TrayMessage (NIM_MODIFY, NULL, icon);
}

void BlinkIcon (
	int	duration)		/* -2 = change icon (called from timer) */
					/* -1 = blinking off, 0 = indefinite */
{
static	int	state = -1;		/* Current duration state */
static	int	icon = WORKING_ICON;	/* Current icon */

/* If this is the first blinking call, start the timer */

	if (state == -1) {
		if (duration < 0) return;
		((CPrime95App *)AfxGetApp())->m_pMainWnd->SetTimer (363, 1000, &TimerCallback);
	}

/* Remember how long we are to blink */

	if (duration >= 0)
		state = duration;

/* If this is the last blinking call, kill the timer */

	else if (duration == -1 || (state && --state == 0)) {
		((CPrime95App *)AfxGetApp())->m_pMainWnd->KillTimer (363);
		state = -1;
		icon = IDLE_ICON;
	}

/* Toggle the icon */

	icon = (icon == IDLE_ICON) ? WORKING_ICON : IDLE_ICON;
	ChangeIcon (icon);
}

void BroadcastMessage (
	char	*message)
{
	char	filename[33];
	int	fd, len;

/* Generate broadcast message file name */

        strcpy (filename, "bcastmsg");
        strcat (filename, EXTENSION);

/* If this is a call to check if a broadcast message exists, then do so */

	if (message == NULL) {
		if (! fileExists (filename)) return;
	}

/* Otherwise, this is a new message - write it to the file */

	else {
		fd = _open (filename, _O_TEXT | _O_RDWR | _O_CREAT | _O_APPEND, 0666);
		if (fd < 0) return;
		_write (fd, message, strlen (message));
		_close (fd);
        }

/* Read in the message from the file */

	fd = _open (filename, _O_TEXT | _O_RDONLY, 0);
	if (fd < 0) return;
	BROADCAST_MESSAGE = (char *) malloc (1024);
	len = _read (fd, BROADCAST_MESSAGE, 1024);
	_close (fd);
	if (len < 0) return;
	BROADCAST_MESSAGE[len] = 0;

/* Blink the icon and display the message if prime95 is visible. */
/* Otherwise, wait until prime95 is activated to display the message. */

	if (!HIDE_ICON) BlinkIcon (0);
	if (AfxGetApp()->m_pMainWnd->IsWindowVisible ())
		 AfxGetApp()->m_pMainWnd->PostMessage (WM_COMMAND, USR_BROADCAST, 0);
}

void CPrime95App::TrayMessage (UINT message, LPCSTR prompt, UINT icon)
{
	NOTIFYICONDATA tnd;

	switch (message) {
	case NIM_ADD :
		WINDOWS95_TRAY_ADD = 0;
		tnd.uFlags = NIF_MESSAGE | NIF_ICON | NIF_TIP;
		break;
	case NIM_MODIFY :
		tnd.uFlags = 0;
		if (prompt != NULL) tnd.uFlags |= NIF_TIP;
		if (icon) tnd.uFlags |= NIF_ICON;
		break;
	case NIM_DELETE :
		tnd.uFlags = 0;
	}

	if (prompt)
		lstrcpy (tnd.szTip, prompt);
	else
		lstrcpy (tnd.szTip, "Prime95");

	tnd.uID = 352;
	tnd.cbSize = sizeof(tnd);
	tnd.hWnd = m_pMainWnd->m_hWnd;
	tnd.uCallbackMessage = MYWM_TRAYMESSAGE;
	tnd.hIcon = LoadIcon (icon ? icon : IDR_MAINFRAME);

	if (Shell_NotifyIcon(message, &tnd)) return;

// The shell command failed.  This could be due to running
// as a Windows 95 service.  Set a flag to try and add the
// icon at a later time

	if (message == NIM_ADD)
		WINDOWS95_TRAY_ADD = 1;
}

/////////////////////////////////////////////////////////////////////////////
// CAboutDlg dialog used for App About

class CAboutDlg : public CDialog
{
public:
	CAboutDlg();

// Dialog Data
	//{{AFX_DATA(CAboutDlg)
	enum { IDD = IDD_ABOUTBOX };
	//}}AFX_DATA

	// ClassWizard generated virtual function overrides
	//{{AFX_VIRTUAL(CAboutDlg)
	protected:
	virtual void DoDataExchange(CDataExchange* pDX);    // DDX/DDV support
	//}}AFX_VIRTUAL

// Implementation
protected:
	//{{AFX_MSG(CAboutDlg)
		// No message handlers
	//}}AFX_MSG
	DECLARE_MESSAGE_MAP()
};

CAboutDlg::CAboutDlg() : CDialog(CAboutDlg::IDD)
{
	//{{AFX_DATA_INIT(CAboutDlg)
	//}}AFX_DATA_INIT
}

void CAboutDlg::DoDataExchange(CDataExchange* pDX)
{
	CDialog::DoDataExchange(pDX);
	//{{AFX_DATA_MAP(CAboutDlg)
	//}}AFX_DATA_MAP
}

BEGIN_MESSAGE_MAP(CAboutDlg, CDialog)
	//{{AFX_MSG_MAP(CAboutDlg)
		// No message handlers
	//}}AFX_MSG_MAP
END_MESSAGE_MAP()

// App command to run the dialog
void CPrime95App::OnAppAbout()
{
	CAboutDlg aboutDlg;
	aboutDlg.DoModal();
}

// In ExitInstance clean up previously allocated mutex
int CPrime95App::ExitInstance() 
{
	if( NULL != g_hMutexInst ) {
		CloseHandle( g_hMutexInst );
	}
	
	return CWinApp::ExitInstance();
}


/////////////////////////////////////////////////////////////////////////////
// CPrime95App commands


/////////////////////////////////////////////////////////////////////////////
// My application-wide stuff went here
/////////////////////////////////////////////////////////////////////////////

int volatile THREAD_STOP = 0;
int	EXIT_IN_PROGRESS = 0;

char	*lines[NumLines] = {NULL};
int	charHeight = 0;

int	WINDOWS95_SERVICE = 0;
int	WINDOWS95_A_SWITCH = 0;
LONG	WM_ENDSESSION_LPARAM = 0;
int	WINDOWS95_TRAY_ADD = 0;
int	CHECK_BATTERY = 0;
int	STOPPED_ON_BATTERY = 0;
char	*BROADCAST_MESSAGE = NULL;

