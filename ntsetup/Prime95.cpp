// Prime95.cpp : Defines the class behaviors for the application.
//

#include "stdafx.h"
#include "Prime95.h"

#include "MainFrm.h"
#include "Prime95Doc.h"
#include "Prime95View.h"

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
	int	named_ini_files = -1;
	char	*p;

	// Standard initialization
	// If you are not using these features and wish to reduce the size
	//  of your final executable, you should remove from the following
	//  the specific initialization routines you do not need.

#ifdef _AFXDLL
	Enable3dControls();			// Call this when using MFC in a shared DLL
#else
	Enable3dControlsStatic();	// Call this when linking to MFC statically
#endif

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

/* Determine the names of the INI files and read them */

	nameIniFiles (named_ini_files);
	readIniFiles ();

// On first run, get user name and email address
// before contacting server for a work assignment

	if (USE_PRIMENET && USERID[0] == 0) {
		STARTUP_IN_PROGRESS = 1;
		m_pMainWnd->PostMessage (WM_COMMAND, ID_RANGE_USERINFORMATION, 0);
	}

// If a broadcast message from the server has been received but
// never viewed by the user, then try to display it now.

	BroadcastMessage (NULL);
	return TRUE;
}



void CALLBACK EXPORT TimerCallback (
	HWND	hWnd,		//handle of CWnd that called SetTimer
	UINT	nMsg,		//WM_TIMER
	UINT	nIDEvent,	//timer identification
	DWORD	dwTime)		//system time
{
	BlinkIcon (-2);
}

void ChangeIcon (
	int	icon_id)
{
	CPrime95App *app;
static	UINT	icon = IDR_MAINFRAME;

	app = (CPrime95App *) AfxGetApp();
	if (icon_id >= 0)
		icon = (icon_id == WORKING_ICON) ? IDR_MAINFRAME : IDI_YELLOW_ICON;
	app->m_pMainWnd->SetIcon (app->LoadIcon (icon), 1);
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

	BlinkIcon (0);
	if (AfxGetApp()->m_pMainWnd->IsWindowVisible ())
		 AfxGetApp()->m_pMainWnd->PostMessage (WM_COMMAND, IDM_BROADCAST, 0);
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

/////////////////////////////////////////////////////////////////////////////
// CPrime95App commands


/////////////////////////////////////////////////////////////////////////////
// My application-wide stuff went here
/////////////////////////////////////////////////////////////////////////////

char	*lines[NumLines] = {NULL};
int	charHeight = 0;

char	*BROADCAST_MESSAGE = NULL;
