// Prime95.h : main header file for the PRIME95 application
//

#ifndef __AFXWIN_H__
	#error include 'stdafx.h' before including this file for PCH
#endif

#include "resource.h"       // main symbols
#include "cpuid.h"
#include "speed.h"
//#define SERVER_TESTING
#define NO_GUI		0
#define EXTERNC		extern "C"
#define PORT	5
#include "gwnum.h"
#include "commona.h"
#include "commonc.h"
#include "comm95a.h"
#include "comm95c.h"
#include "primenet.h"

/////////////////////////////////////////////////////////////////////////////
// CPrime95App:
// See Prime95.cpp for the implementation of this class
//

class CPrime95App : public CWinApp
{
public:
	CPrime95App();
	void TrayMessage (UINT, LPCSTR);

// Overrides
	// ClassWizard generated virtual function overrides
	//{{AFX_VIRTUAL(CPrime95App)
	public:
	virtual BOOL InitInstance();
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

const int NumLines = 60;
extern char *lines[NumLines];
extern int charHeight;

extern char *BROADCAST_MESSAGE;
