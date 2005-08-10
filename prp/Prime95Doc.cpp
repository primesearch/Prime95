// Prime95Doc.cpp : implementation of the CPrime95Doc class
//

#include "stdafx.h"
#include "MainFrm.h"
#include "Prime95.h"
#include "Prime95Doc.h"

#include <direct.h>
#include "math.h"

#include "AffinityDlg.h"
#include "ContentsDlg.h"
#include "CpuDlg.h"
#include "PreferencesDlg.h"
#include "Priority.h"
#include "TestDlg.h"

#include <ctype.h>
#include <dos.h>
#include <fcntl.h>
#include <io.h>
#include <memory.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <sys/timeb.h>

#ifdef _DEBUG
#define new DEBUG_NEW
#undef THIS_FILE
static char THIS_FILE[] = __FILE__;
#endif

#define OP_CONTINUE	1
struct thread_info {
	int	op;		// Opcode defined above
	CPrime95Doc *doc;	// Ptr to main document
} thread_pkt;

/////////////////////////////////////////////////////////////////////////////
// CPrime95Doc

IMPLEMENT_DYNCREATE(CPrime95Doc, CDocument)

BEGIN_MESSAGE_MAP(CPrime95Doc, CDocument)
	//{{AFX_MSG_MAP(CPrime95Doc)
	ON_COMMAND(IDM_CONTINUE, OnContinue)
	ON_UPDATE_COMMAND_UI(IDM_CONTINUE, OnUpdateContinue)
	ON_COMMAND(IDM_STOP, OnStop)
	ON_UPDATE_COMMAND_UI(IDM_STOP, OnUpdateStop)
	ON_COMMAND(IDM_ERRCHK, OnErrchk)
	ON_UPDATE_COMMAND_UI(IDM_ERRCHK, OnUpdateErrchk)
	ON_COMMAND(IDM_CPU, OnCpu)
	ON_COMMAND(IDM_PREFERENCES, OnPreferences)
	ON_COMMAND(IDM_TEST, OnTest)
	ON_COMMAND(IDM_TRAY, OnTray)
	ON_UPDATE_COMMAND_UI(IDM_TRAY, OnUpdateTray)
	ON_COMMAND(IDM_HIDE, OnHide)
	ON_UPDATE_COMMAND_UI(IDM_HIDE, OnUpdateHide)
	ON_COMMAND(IDM_PRIORITY, OnPriority)
	ON_COMMAND(IDM_SERVICE, OnService)
	ON_UPDATE_COMMAND_UI(IDM_SERVICE, OnUpdateService)
	ON_COMMAND(IDM_AFFINITY, OnAffinity)
	ON_UPDATE_COMMAND_UI(IDM_AFFINITY, OnUpdateAffinity)
	ON_UPDATE_COMMAND_UI(IDM_TEST, OnUpdateTest)
	ON_COMMAND(ID_HELP_CONTENTS, OnHelpContents)
	//}}AFX_MSG_MAP
END_MESSAGE_MAP()

/////////////////////////////////////////////////////////////////////////////
// CPrime95Doc construction/destruction

CPrime95Doc *OUTPUT_STR_HACK;

CPrime95Doc::CPrime95Doc()
{
	OUTPUT_STR_HACK = this;
}

CPrime95Doc::~CPrime95Doc()
{
	int	x;
	for (x = 0; x < NumLines; x++)
		if (lines[x] != NULL) free (lines[x]);
}

BOOL CPrime95Doc::OnNewDocument()
{
	if (!CDocument::OnNewDocument())
		return FALSE;

	// TODO: add reinitialization code here
	// (SDI documents will reuse this document)

	return TRUE;
}

/////////////////////////////////////////////////////////////////////////////
// CPrime95Doc serialization

void CPrime95Doc::Serialize(CArchive& ar)
{
	if (ar.IsStoring())
	{
		// TODO: add storing code here
	}
	else
	{
		// TODO: add loading code here
	}
}

/////////////////////////////////////////////////////////////////////////////
// CPrime95Doc diagnostics

#ifdef _DEBUG
void CPrime95Doc::AssertValid() const
{
	CDocument::AssertValid();
}

void CPrime95Doc::Dump(CDumpContext& dc) const
{
	CDocument::Dump(dc);
}
#endif //_DEBUG

/////////////////////////////////////////////////////////////////////////////
// CPrime95Doc commands

void CPrime95Doc::OnCloseDocument() 
{

// Set flag indicating we are exiting.  Needed because setting the icon
// while in the sleep 50 loop causes a hang.

	EXIT_IN_PROGRESS = 1;

// Stop background thread before exiting

	if (THREAD_ACTIVE) {
		OnStop ();
		while (THREAD_STOP) Sleep (50);
	}

// Remember the main window's size and position
// Note: INI_FILE may not be initialized when this instance
// was used to activate another instance

	CWinApp* pApp = AfxGetApp();
	WINDOWPLACEMENT wp;
	if (pApp->m_pMainWnd && INI_FILE[0]) {
		pApp->m_pMainWnd->GetWindowPlacement (&wp);
		IniWriteInt (INI_FILE, "Left", wp.rcNormalPosition.left);
		IniWriteInt (INI_FILE, "Top", wp.rcNormalPosition.top);
		IniWriteInt (INI_FILE, "Right", wp.rcNormalPosition.right);
		IniWriteInt (INI_FILE, "Bottom", wp.rcNormalPosition.bottom);
	}

// Finish closing

	CDocument::OnCloseDocument();
}

void CPrime95Doc::OnUpdateTest(CCmdUI* pCmdUI) 
{
	pCmdUI->Enable (!THREAD_ACTIVE);
}

void CPrime95Doc::OnTest() 
{
	CTestDlg dlg;
	char	buf[2000];

	dlg.m_work = IniGetInt (INI_FILE, "Work", 0);
	IniGetString (INI_FILE, "PgenInputFile", buf, sizeof(buf), NULL);
	dlg.m_pgen_input = buf;
	IniGetString (INI_FILE, "PgenOutputFile", buf, sizeof(buf), NULL);
	dlg.m_pgen_output = buf;
	dlg.m_pgen_line = IniGetInt (INI_FILE, "PgenLine", 1);
	IniGetString (INI_FILE, "Expr", buf, sizeof(buf), NULL);
	dlg.m_expr = buf;
	IniGetString (INI_FILE, "ExprFile", buf, sizeof(buf), NULL);
	dlg.m_exprfile = buf;
	dlg.m_exprfile_line = IniGetInt (INI_FILE, "ExprFileLine", 1);
dlg.m_expr = "Not yet implemented";
dlg.m_exprfile = "Not yet implemented";
	if (dlg.DoModal () == IDOK) {
		IniWriteInt (INI_FILE, "Work", dlg.m_work);
		IniWriteString (INI_FILE, "PgenInputFile", (char *)(LPCTSTR) dlg.m_pgen_input);
		IniWriteString (INI_FILE, "PgenOutputFile", (char *)(LPCTSTR) dlg.m_pgen_output);
		IniWriteInt (INI_FILE, "PgenLine", dlg.m_pgen_line);
		IniWriteString (INI_FILE, "Expr", (char *)(LPCTSTR) dlg.m_expr);
		IniWriteString (INI_FILE, "ExprFile", (char *)(LPCTSTR) dlg.m_exprfile);
		IniWriteInt (INI_FILE, "ExprFileLine", dlg.m_exprfile_line);
		IniWriteInt (INI_FILE, "WorkDone", 0);
		OnContinue ();
	}
}

void CPrime95Doc::OnUpdateContinue(CCmdUI* pCmdUI) 
{
	pCmdUI->Enable (! THREAD_ACTIVE);
}

void CPrime95Doc::OnContinue() 
{
	CWinThread *thread;

	// Set flag to check laptop battery status
	// This ugly 999 code happens when we are auto-restarting
	// from a POWER_BROADCAST message.  By setting STOPPED_ON_BATTERY
	// to one, escapeCheck will not output a message should the
	// POWER_BROADCAST message turn out to be spurious.
	CHECK_BATTERY = 1;
	if (STOPPED_ON_BATTERY == 999)
		STOPPED_ON_BATTERY = 1;
	else
		STOPPED_ON_BATTERY = 0;

	// Start the thread
	thread_pkt.op = OP_CONTINUE;
	thread_pkt.doc = this;
	thread = AfxBeginThread (threadDispatch, NULL);
}

void CPrime95Doc::OnUpdateStop(CCmdUI* pCmdUI) 
{
	pCmdUI->Enable (THREAD_ACTIVE && !THREAD_STOP);
}

void CPrime95Doc::OnStop() 
{
	THREAD_STOP = 1;
	SetPriorityClass (GetCurrentProcess (), NORMAL_PRIORITY_CLASS);
	SetThreadPriority (WORKER_THREAD, THREAD_PRIORITY_ABOVE_NORMAL);
}

// Advanced Menu

void CPrime95Doc::OnUpdateErrchk(CCmdUI* pCmdUI) 
{
	pCmdUI->SetCheck (ERRCHK);
}

void CPrime95Doc::OnErrchk() 
{
	ERRCHK = !ERRCHK;
	IniWriteInt (INI_FILE, "ErrorCheck", ERRCHK);
}

void CPrime95Doc::OnPriority() 
{
	CPriority dlg;

	dlg.m_priority = PRIORITY;
	if (dlg.DoModal () == IDOK) {
		if (PRIORITY != dlg.m_priority) {
			Restart1 ();
			PRIORITY = dlg.m_priority;
			IniWriteInt (INI_FILE, "Priority", PRIORITY);
			Restart2 ();
		}
	}
}

void CPrime95Doc::OnUpdateAffinity(CCmdUI* pCmdUI) 
{
	pCmdUI->Enable (!isWindows95 ());
}

void CPrime95Doc::OnAffinity() 
{
	CAffinityDlg dlg;

	if (CPU_AFFINITY == 99) {
		dlg.m_all_cpus = 1;
		dlg.m_cpu = 0;
	} else {
		dlg.m_all_cpus = 0;
		dlg.m_cpu = CPU_AFFINITY;
	}
	if (dlg.DoModal () == IDOK) {
		if (dlg.m_all_cpus) CPU_AFFINITY = 99;
		else CPU_AFFINITY = dlg.m_cpu;
		IniWriteInt (INI_FILE, "Affinity", CPU_AFFINITY);
	}
}


// Options menu

void CPrime95Doc::OnCpu() 
{
	CCpuDlg dlg;
	char	buf[512];

	getCpuDescription (buf, 0);
	dlg.m_cpu_info = buf;
	if (dlg.DoModal () == IDOK) {
	}
}

void CPrime95Doc::OnPreferences() 
{
	CPreferencesDlg dlg;

	dlg.m_iter = ITER_OUTPUT;
	dlg.m_r_iter = ITER_OUTPUT_RES;
	dlg.m_disk_write_time = DISK_WRITE_TIME;
	dlg.m_cumulative = CUMULATIVE_TIMING;
	dlg.m_backup = TWO_BACKUP_FILES;
	dlg.m_battery = RUN_ON_BATTERY;
	if (dlg.DoModal () == IDOK) {
		ITER_OUTPUT = dlg.m_iter;
		ITER_OUTPUT_RES = dlg.m_r_iter;
		DISK_WRITE_TIME = dlg.m_disk_write_time;
		TWO_BACKUP_FILES = dlg.m_backup;
		RUN_ON_BATTERY = dlg.m_battery;
		CUMULATIVE_TIMING = dlg.m_cumulative;
		IniWriteInt (INI_FILE, "OutputIterations", ITER_OUTPUT);
		IniWriteInt (INI_FILE, "ResultsFileIterations", ITER_OUTPUT_RES);
		IniWriteInt (INI_FILE, "DiskWriteTime", DISK_WRITE_TIME);
		IniWriteInt (INI_FILE, "TwoBackupFiles", TWO_BACKUP_FILES);
		IniWriteInt (INI_FILE, "CumulativeTiming", CUMULATIVE_TIMING);
		IniWriteInt (INI_FILE, "RunOnBattery", RUN_ON_BATTERY);
		CHECK_BATTERY = 1;
	}
}

void CPrime95Doc::OnUpdateTray(CCmdUI* pCmdUI) 
{
	pCmdUI->SetCheck (TRAY_ICON);
}

void CPrime95Doc::OnTray() 
{
	CPrime95App* pApp = (CPrime95App *) AfxGetApp();
	TRAY_ICON = ! TRAY_ICON;
	if (TRAY_ICON) {
		HIDE_ICON = 0;
		pApp->TrayMessage (NIM_ADD, NULL, 0);
		ChangeIcon (-1);
	} else {
		pApp->TrayMessage (NIM_DELETE, NULL, 0);
	}
	IniWriteInt (INI_FILE, "HideIcon", HIDE_ICON);
	IniWriteInt (INI_FILE, "TrayIcon", TRAY_ICON);
}

void CPrime95Doc::OnUpdateHide(CCmdUI* pCmdUI) 
{
	pCmdUI->SetCheck (HIDE_ICON);
}

void CPrime95Doc::OnHide() 
{
	CPrime95App* pApp = (CPrime95App *) AfxGetApp();

	HIDE_ICON = ! HIDE_ICON;
	if (HIDE_ICON) {
		if (TRAY_ICON) pApp->TrayMessage (NIM_DELETE, NULL, 0);
		TRAY_ICON = 0;
	}
	IniWriteInt (INI_FILE, "HideIcon", HIDE_ICON);
	IniWriteInt (INI_FILE, "TrayIcon", TRAY_ICON);
}

// When running as an NT service we can delete the service (it will take
// effect when the service is stopped), but we cannot recreate the service
// until the next time prime95 is run.  Thus, disable this menu choice once
// an NT service has turned this option off.

void CPrime95Doc::OnUpdateService(CCmdUI* pCmdUI) 
{
	pCmdUI->Enable (!NTSERVICENAME[0] || WINDOWS95_SERVICE);
	pCmdUI->SetCheck (WINDOWS95_SERVICE);
}

void CPrime95Doc::OnService() 
{
	WINDOWS95_SERVICE = !WINDOWS95_SERVICE;
	IniWriteInt (INI_FILE, "Windows95Service", WINDOWS95_SERVICE);
	Service95 ();
}

void CPrime95Doc::OnHelpContents() 
{
	ContentsDlg dlg;
	char	buf[100000];
	int	fd, len;

	fd = _open ("readme.txt", _O_BINARY | _O_RDONLY);
	if (fd) {
		len = _read (fd, buf, sizeof (buf));
		_close (fd);
	}
	if (len <= 0) {
		dlg.m_helptext = "No Help Available.";
	} else {
		buf[len] = 0;
		dlg.m_helptext = buf;
	}

	dlg.DoModal ();
}

/////////////////////////////////////////////////////////////////////////////
// CPrime95Doc private routines

// Crude routine to stop and restart the worker thread

int	RESTARTING = FALSE;

void CPrime95Doc::Restart1 ()
{
	if (!THREAD_ACTIVE) return;
	RESTARTING = TRUE;
	OutputStr ("Stopping and restarting using new settings.\n");
	OnStop ();
	while (THREAD_ACTIVE) {
		MSG m;
		while (PeekMessage (&m, NULL, 0, 0, PM_NOREMOVE))
			AfxGetApp()->PumpMessage ();
		Sleep (100);
	}
}

void CPrime95Doc::Restart2 ()
{
	if (RESTARTING) {
		RESTARTING = FALSE;
		OnContinue ();
	}
}

void CPrime95Doc::LineFeed ()
{
	char	*p;

// Scroll the line on the screen

	if (!Replacing) {
		p = lines[NumLines-1];
		memmove (&lines[1], &lines[0], (NumLines-1)*sizeof(char *));
		lines[0] = p;
		*p = 0;
	}
	UpdateAllViews (NULL, Replacing);
	Replacing = FALSE;
}

void OutputStr (char *str)
{
	OUTPUT_STR_HACK->OutputStr (str);
}

void ReplaceableLine (int type)
{
	OUTPUT_STR_HACK->ReplaceableLine (type);
}

void title (char *str)
{
	OUTPUT_STR_HACK->title (str);
}

void flashWindowAndBeep ()
{
	CWinApp* pApp = AfxGetApp();
	pApp->m_pMainWnd->FlashWindow (TRUE);
	MessageBeep (0xFFFFFFFF);
}

void CPrime95Doc::OutputStr (
	char	*str)
{
	if (lines[0] == NULL) {
		for (int i = 0; i < NumLines; i++) {
			lines[i] = (char *) malloc (128);
			*lines[i] = 0;
		}
		Replacing = FALSE;
	}
	char *origp, *p;
	if (Replacing) origp = ReplLine;
	else origp = lines[0];
	p = origp + strlen (origp);
	for ( ; *str; str++) {
		if (*str == '\n') *p = 0, LineFeed (), p = lines[0];
		else if (p - origp < 127) *p++ = *str;
	}
	*p = 0;
}

void CPrime95Doc::ReplaceableLine (
	int	type)
{
	if (type == 1)
		ReplLine = lines[1];
	else {
		Replacing = TRUE;
		ReplLine[0] = 0;
	}
}


/////////////////////////////////////////////////////////////////////////////
// CPrime95Doc public routines

#define PORT	1
#include "prp.c"
#include "prp95.c"

UINT threadDispatch (
	LPVOID stuff)
{

// Set thread-active flags

	THREAD_ACTIVE = TRUE;
	THREAD_STOP = 0;

// Stall if we've just booted (within 5 minutes of Windows starting)

	if (GetTickCount () < 300000) {
		int	delay;
		delay = IniGetInt (INI_FILE, "BootDelay", 90);
		delay -= GetTickCount () / 1000;
		if (delay > 0) {
			char buf[50];
			sprintf (buf, "Waiting %d seconds for boot to complete.\n", delay);
			OutputStr (buf);
			Sleep (delay * 1000);
		}
	}

// Change the icon

	ChangeIcon (WORKING_ICON);

// Dispatch to the correct code

	switch (thread_pkt.op) {
	case OP_CONTINUE:
		primeContinue ();
		break;
	}

// Output informative message

	if (IniGetInt (INI_FILE, "WorkDone", 1)) {
		OutputStr ("Choose Test / Input Data to test more numbers.\n");
	}
	else if (!STOPPED_ON_BATTERY && !RESTARTING) {
		OutputStr ("Execution halted.\n");
		OutputStr ("Choose Test/Continue to restart.\n");
	}

// Clear thread-active flags

	THREAD_ACTIVE = FALSE;
	THREAD_STOP = 0;
	ChangeIcon (IDLE_ICON);
	title ("IDLE");

// Thread complete

	return (0);
}


/* Return TRUE if we should stop calculating */
/* Perform other misceallanous tasks */

int escapeCheck ()
{
	// If we couldn't add the icon to the task bar, then keep
	// trying until we finally succeed!
	if (WINDOWS95_TRAY_ADD) {
		CPrime95App *pApp = (CPrime95App *)AfxGetApp();
		pApp->TrayMessage (NIM_ADD, NULL, 0);
	}
	// If we aren't supposed to run on the battery and we're not
	// connected to the AC power, then stop processing.
	if (CHECK_BATTERY) {
		SYSTEM_POWER_STATUS power;
		if (!RUN_ON_BATTERY &&
		    GetSystemPowerStatus (&power) &&
		    power.ACLineStatus != 1) {
			if (STOPPED_ON_BATTERY == 0) {
				OutputStr ("Processing stopped while on battery power.\n");
				STOPPED_ON_BATTERY = 1;
			}
			return (TRUE);
		}
		CHECK_BATTERY = 0;
		STOPPED_ON_BATTERY = 0;
	}

/* Return TRUE if the thread must stop executing */

	return (THREAD_STOP);
}

/* Put a title on the main window */

void CPrime95Doc::title (
	char	*str)
{
	char buf[80];
	CWinApp* pApp = AfxGetApp();
static	int was_iconic = TRUE;

	sprintf (buf, "PRP - %s", str);
	if (pApp->m_pMainWnd) {
		if (TRAY_ICON)
			((CPrime95App *)pApp)->
				TrayMessage (NIM_MODIFY, buf, 0);
		if (pApp->m_pMainWnd->IsIconic ()) {
			pApp->m_pMainWnd->SetWindowText (buf);
			was_iconic = TRUE;
		} else if (was_iconic) {
			pApp->m_pMainWnd->SetWindowText ("PRP");
			was_iconic = FALSE;
		}
	}
}
