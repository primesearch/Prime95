// Prime95.h : main header file for the PRIME95 application
//

#ifndef __AFXWIN_H__
	#error include 'stdafx.h' before including this file for PCH
#endif

#include "resource.h"       // main symbols
#include "EditDropFiles.h"
#define NO_GUI		0
#define EXTERNC		extern "C"
#include "cpuid.h"
#include "gwnum.h"
#include "prp.h"
#include "prp95.h"

/////////////////////////////////////////////////////////////////////////////
// CPrime95App:
// See Prime95.cpp for the implementation of this class
//

class CPrime95App : public CWinApp
{
public:
	CPrime95App();
	void TrayMessage (UINT, LPCSTR, UINT);

// Overrides
	// ClassWizard generated virtual function overrides
	//{{AFX_VIRTUAL(CPrime95App)
	public:
	virtual BOOL InitInstance();
	virtual int ExitInstance();
	//}}AFX_VIRTUAL

// Implementation

	//{{AFX_MSG(CPrime95App)
	afx_msg void OnAppAbout();
		// NOTE - the ClassWizard will add and remove member functions here.
		//    DO NOT EDIT what you see in these blocks of generated code !
	//}}AFX_MSG
	DECLARE_MESSAGE_MAP()
};


/////////////////////////////////////////////////////////////////////////////
// My non-MFC stuff went here
/////////////////////////////////////////////////////////////////////////////

// Global variables

extern HANDLE g_hMutexInst;

extern int THREAD_ACTIVE;		// True if a thread is active
extern int volatile THREAD_STOP;	// one if thread should stop with msg
					// two if thread should stop quietly
extern int EXIT_IN_PROGRESS;		// True if we are exiting

const int NumLines = 60;
extern char *lines[NumLines];
extern int charHeight;

extern int WINDOWS95_SERVICE;		// True if we're running as a Win95 service
extern int WINDOWS95_A_SWITCH;		// Value of the -A command line switch
extern LONG WM_ENDSESSION_LPARAM;	// LPARAM of WM_ENDSESSION message
extern int WINDOWS95_TRAY_ADD;		// True if we need to add the icon
					// to the shell once the user logs in
extern int CHECK_BATTERY;		// TRUE if we need to test battery status
extern int STOPPED_ON_BATTERY;		// TRUE if we need to autocontinue on a power
					// change event
// Variables used to communicate with the NT service code

extern "C" char NTSERVICENAME[32];	// name of the NT service
extern "C" HWND MAINFRAME_HWND;		// Handle of main frame window

// Internal routines

UINT threadDispatch (LPVOID);
void Service95 ();
