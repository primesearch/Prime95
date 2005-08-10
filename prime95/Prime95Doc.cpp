// Prime95Doc.cpp : implementation of the CPrime95Doc class
//

#include "stdafx.h"
#include "MainFrm.h"
#include "Prime95.h"
#include "Prime95Doc.h"

#include <direct.h>
#include "math.h"

#include "AffinityDlg.h"
#include "CpuDlg.h"
#include "EcmDlg.h"
#include "ManualCommDlg.h"
#include "Password.h"
#include "Pminus1Dlg.h"
#include "PreferencesDlg.h"
#include "PrimeNetDlg.h"
#include "Priority.h"
#include "ServerDlg.h"
#include "TestDlg.h"
#include "TimeDlg.h"
#include "TortureDlg.h"
#include "UnreserveDlg.h"
#include "UserDlg.h"
#include "VacationDlg.h"
#include "WelcomeDlg.h"

#ifdef _DEBUG
#define new DEBUG_NEW
#undef THIS_FILE
static char THIS_FILE[] = __FILE__;
#endif

#define OP_CONTINUE	1
#define OP_TIME		2
#define OP_BENCH	3
#define OP_TORTURE	5
struct thread_info {
	int	op;		// Opcode defined above
	unsigned long p;	// Prime to test
	unsigned long p2;	// End prime
	unsigned short b1;	// Start # bits
	unsigned short b2;	// End # bits
	int	db;		// Search DB
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
	ON_UPDATE_COMMAND_UI(IDM_TEST, OnUpdateTest)
	ON_COMMAND(IDM_TEST, OnTest)
	ON_UPDATE_COMMAND_UI(IDM_TIME, OnUpdateTime)
	ON_COMMAND(IDM_TIME, OnTime)
	ON_COMMAND(ID_RANGE_STATUS, OnRangeStatus)
	ON_COMMAND(IDM_PASSWORD, OnPassword)
	ON_UPDATE_COMMAND_UI(IDM_PASSWORD, OnUpdatePassword)
	ON_UPDATE_COMMAND_UI(ID_HELP_FINDER, OnUpdateHelpFinder)
	ON_COMMAND(IDM_TRAY, OnTray)
	ON_UPDATE_COMMAND_UI(IDM_TRAY, OnUpdateTray)
	ON_COMMAND(IDM_HIDE, OnHide)
	ON_UPDATE_COMMAND_UI(IDM_HIDE, OnUpdateHide)
	ON_COMMAND(IDM_TORTURE, OnTorture)
	ON_UPDATE_COMMAND_UI(IDM_TORTURE, OnUpdateTorture)
	ON_COMMAND(ID_RANGE_USERINFORMATION, OnRangeUserinformation)
	ON_COMMAND(IDM_PRIMENET, OnPrimenet)
	ON_COMMAND(IDM_PRIORITY, OnPriority)
	ON_UPDATE_COMMAND_UI(IDM_PRIORITY, OnUpdatePriority)
	ON_COMMAND(IDM_SERVER, OnServer)
	ON_UPDATE_COMMAND_UI(IDM_SERVER, OnUpdateServer)
	ON_COMMAND(IDM_QUIT, OnQuitGimps)
	ON_COMMAND(IDM_VACATION, OnVacation)
	ON_COMMAND(IDM_SERVICE, OnService)
	ON_UPDATE_COMMAND_UI(IDM_SERVICE, OnUpdateService)
	ON_COMMAND(ID_MANUALCOMM, OnManualcomm)
	ON_UPDATE_COMMAND_UI(ID_MANUALCOMM, OnUpdateManualcomm)
	ON_UPDATE_COMMAND_UI(ID_RANGE_USERINFORMATION, OnUpdateRangeUserinformation)
	ON_COMMAND(IDM_ECM, OnEcm)
	ON_UPDATE_COMMAND_UI(IDM_ECM, OnUpdateEcm)
	ON_COMMAND(IDM_AFFINITY, OnAffinity)
	ON_UPDATE_COMMAND_UI(IDM_AFFINITY, OnUpdateAffinity)
	ON_COMMAND(IDM_PMINUS1, OnPminus1)
	ON_UPDATE_COMMAND_UI(IDM_PMINUS1, OnUpdatePminus1)
	ON_COMMAND(USR_WELCOME, OnWelcome)
	ON_COMMAND(USR_TORTURE, OnUsrTorture)
	ON_COMMAND(USR_BROADCAST, OnBroadcast)
	ON_COMMAND(IDM_UNRESERVE, OnUnreserve)
	ON_UPDATE_COMMAND_UI(IDM_UNRESERVE, OnUpdateUnreserve)
	ON_UPDATE_COMMAND_UI(IDM_VACATION, OnUpdateVacation)
	ON_UPDATE_COMMAND_UI(IDM_QUIT, OnUpdateQuit)
	ON_COMMAND(IDM_BENCHMARK, OnBenchmark)
	ON_UPDATE_COMMAND_UI(IDM_BENCHMARK, OnUpdateBenchmark)
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

// Free the networking library

	UnloadPrimeNet ();

// Finish closing

	CDocument::OnCloseDocument();
}

// Range menu

#define MSG_BIGONES	"A 500 MHz Pentium-III computer will take a full year to test just one 10,000,000 digit number.  Your chance of finding a new prime is roughly 1 in 250,000.  Read the prize rules at http://www.mersenne.org/prize.htm.  If you accept these rules and still want to search for 10,000,000 digit primes, then click Yes."

void CPrime95Doc::OnPrimenet() 
{
	PrimenetDlg dlg;
	short	work_pref;

	dlg.m_primenet = USE_PRIMENET;
	dlg.m_dialup = DIAL_UP;
	dlg.m_work = DAYS_OF_WORK;
	if (WORK_PREFERENCE == 0) {
		dlg.m_work_dflt = 1;
		work_pref = default_work_type ();
	} else {
		dlg.m_work_dflt = 0;
		work_pref = WORK_PREFERENCE;
	}
	dlg.m_bigones = !! (work_pref & PRIMENET_ASSIGN_BIGONES);
	dlg.m_lucas = !! (work_pref & PRIMENET_ASSIGN_TEST);
	dlg.m_pfactor = !! (work_pref & PRIMENET_ASSIGN_PFACTOR);
	dlg.m_factor = !! (work_pref & PRIMENET_ASSIGN_FACTOR);
	dlg.m_dblchk = !! (work_pref & PRIMENET_ASSIGN_DBLCHK);

	if (dlg.DoModal () == IDOK) {
		if (!USE_PRIMENET && dlg.m_primenet) {
			USE_PRIMENET = 1;
			spoolMessage (PRIMENET_MAINTAIN_USER_INFO, NULL);
			spoolMessage (PRIMENET_SET_COMPUTER_INFO, NULL);
			spoolExistingResultsFile ();
		}
		USE_PRIMENET = dlg.m_primenet;
		DIAL_UP = dlg.m_dialup;
		DAYS_OF_WORK = dlg.m_work;
		if (dlg.m_work_dflt)
			WORK_PREFERENCE = 0;
		else {
			work_pref =
				(dlg.m_bigones ? PRIMENET_ASSIGN_BIGONES : 0) +
				(dlg.m_lucas ? PRIMENET_ASSIGN_TEST : 0) +
				(dlg.m_pfactor ? PRIMENET_ASSIGN_PFACTOR : 0) +
				(dlg.m_factor ? PRIMENET_ASSIGN_FACTOR : 0) +
				(dlg.m_dblchk ? PRIMENET_ASSIGN_DBLCHK : 0);
			if (! (work_pref & PRIMENET_ASSIGN_BIGONES) ||
			    WORK_PREFERENCE & PRIMENET_ASSIGN_BIGONES ||
			    AfxMessageBox (MSG_BIGONES, MB_YESNO | MB_ICONQUESTION | MB_DEFBUTTON2) == IDYES)
				WORK_PREFERENCE = work_pref;
		}
		IniWriteInt (INI_FILE, "UsePrimenet", USE_PRIMENET);
		IniWriteInt (INI_FILE, "DialUp", DIAL_UP);
		IniWriteInt (INI_FILE, "DaysOfWork", DAYS_OF_WORK);
		IniWriteInt (INI_FILE, "WorkPreference", WORK_PREFERENCE);
		CHECK_WORK_QUEUE = 1;
		if (!THREAD_ACTIVE &&
		    (STARTUP_IN_PROGRESS || USE_PRIMENET)) {
			STARTUP_IN_PROGRESS = 0;
			OnContinue ();
		}
	} else
		STARTUP_IN_PROGRESS = 0;
}

void CPrime95Doc::OnUpdateQuit(CCmdUI* pCmdUI) 
{
	pCmdUI->Enable (USE_PRIMENET || IniGetNumLines (WORKTODO_FILE));
}

void CPrime95Doc::OnQuitGimps() 
{
	if (!USE_PRIMENET) {
		if (AfxMessageBox (MANUAL_QUIT, MB_YESNO | MB_ICONQUESTION | MB_DEFBUTTON2) == IDYES) {
			OutputBoth ("Quitting GIMPS.\n");
			IniDeleteAllLines (WORKTODO_FILE);
			THREAD_STOP = 1;
			GIMPS_QUIT = 1;
			if (WINDOWS95_SERVICE) OnService ();
		}
	} else {
		int	res;
		res = AfxMessageBox (PRIMENET_QUIT, MB_YESNOCANCEL | MB_ICONQUESTION | MB_DEFBUTTON3);
		if (res == IDYES) {
			OutputBoth ("Quitting GIMPS after current work completes.\n");
			IniWriteInt (INI_FILE, "NoMoreWork", 1);
			if (!THREAD_ACTIVE) OnContinue ();
		}
		if (res == IDNO) {
			OutputBoth ("Quitting GIMPS immediately.\n");
			spoolMessage (999, NULL);
			if (!THREAD_ACTIVE) OnContinue ();
		}
	}
}

void CPrime95Doc::OnUpdateRangeUserinformation(CCmdUI* pCmdUI) 
{
	pCmdUI->Enable (! IniGetInt (INI_FILE, "LockUserInfo", 0));
}

void CPrime95Doc::OnRangeUserinformation() 
{
	CUserDlg dlg;

	dlg.m_name = USER_NAME;
	dlg.m_email = USER_ADDR;
	dlg.m_sendemail = NEWSLETTERS;
	dlg.m_userid = USERID;
	dlg.m_password = USER_PWD;
	dlg.m_compid = COMPID;
	dlg.m_team = 0;
	if (dlg.DoModal () == IDOK) {
		if (strcmp (COMPID, dlg.m_compid) != 0) {
			strcpy (COMPID, dlg.m_compid);
			sanitizeString (COMPID);
			IniWriteString (LOCALINI_FILE, "ComputerID", COMPID);
			spoolMessage (PRIMENET_SET_COMPUTER_INFO, NULL);
		}
		if (OLD_USERID[0] == 0 &&
		    strcmp (USERID, dlg.m_userid) != 0) {
			strcpy (OLD_USERID, USERID);
			strcpy (OLD_USER_PWD, USER_PWD);
			IniWriteString (INI_FILE, "OldUserID", OLD_USERID);
			IniWriteString (INI_FILE, "OldUserPWD", OLD_USER_PWD);
		}
		strcpy (USER_PWD, dlg.m_password);
		sanitizeString (USER_PWD);
		IniWriteString (INI_FILE, "UserPWD", USER_PWD);
		if (strcmp (USER_NAME, dlg.m_name) != 0 ||
		    strcmp (USER_ADDR, dlg.m_email) != 0 ||
		    NEWSLETTERS != dlg.m_sendemail ||
		    dlg.m_team ||
		    strcmp (USERID, dlg.m_userid) != 0) {
			strcpy (USER_NAME, dlg.m_name);
			strcpy (USER_ADDR, dlg.m_email);
			NEWSLETTERS = dlg.m_sendemail;
			strcpy (USERID, dlg.m_userid);
			sanitizeString (USERID);
			IniWriteString (INI_FILE, "UserName", USER_NAME);
			IniWriteString (INI_FILE, "UserEmailAddr", USER_ADDR);
			IniWriteInt (INI_FILE, "Newsletters", NEWSLETTERS);
			IniWriteString (INI_FILE, "UserID", USERID);
			spoolMessage (dlg.m_team ?
				PRIMENET_MAINTAIN_USER_INFO + 0x80 :
				PRIMENET_MAINTAIN_USER_INFO, NULL);
		}
		if (STARTUP_IN_PROGRESS) OnCpu ();
	} else
		STARTUP_IN_PROGRESS = 0;
}


void CPrime95Doc::OnUpdateVacation(CCmdUI* pCmdUI) 
{
	pCmdUI->Enable (USE_PRIMENET);
}

void CPrime95Doc::OnVacation() 
{
	CVacationDlg dlg;

	dlg.m_vacation_days =
		(secondsUntilVacationEnds () + 43200) / 86400;
	dlg.m_computer_on = ON_DURING_VACATION;
	if (dlg.DoModal () == IDOK) {
		if (dlg.m_vacation_days) {
			time (&VACATION_END);
			VACATION_END += dlg.m_vacation_days * 86400;
		} else
			VACATION_END = 0;
		ON_DURING_VACATION = dlg.m_computer_on;
                IniWriteInt (LOCALINI_FILE, "VacationEnd", (long) VACATION_END);
                IniWriteInt (LOCALINI_FILE, "VacationOn", ON_DURING_VACATION);
		if (VACATION_END && !ON_DURING_VACATION)
			IniWriteInt (LOCALINI_FILE, "RollingStartTime", 0);
		next_comm_time = 0;
		UpdateEndDates ();
	}
}

void CPrime95Doc::OnRangeStatus() 
{
	rangeStatus ();
}

void CPrime95Doc::OnUpdateContinue(CCmdUI* pCmdUI) 
{
	pCmdUI->Enable (! THREAD_ACTIVE &&
			(USE_PRIMENET || IniGetNumLines (WORKTODO_FILE)));
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

void CPrime95Doc::OnUpdatePassword(CCmdUI* pCmdUI) 
{
	pCmdUI->Enable (!ADVANCED_ENABLED);
}

void CPrime95Doc::OnPassword() 
{
	Password dlg;

	dlg.m_password = 0;
	if (dlg.DoModal () == IDOK) {
		if (dlg.m_password == 9876) {	// secret code
			ADVANCED_ENABLED = 1;
			IniWriteInt (INI_FILE, "Advanced", 1);
			return;
		}
	}
}

void CPrime95Doc::OnUpdateTest(CCmdUI* pCmdUI) 
{
	pCmdUI->Enable (ADVANCED_ENABLED && ! THREAD_ACTIVE);
}

void CPrime95Doc::OnTest() 
{
	CTestDlg dlg;

	if (dlg.DoModal () == IDOK) {
		IniFileOpen (WORKTODO_FILE, 1);
		IniInsertLineAsInt (WORKTODO_FILE, 1, "AdvancedTest", dlg.m_p);
		OnContinue ();
	}
}

void CPrime95Doc::OnUpdateTime(CCmdUI* pCmdUI) 
{
	pCmdUI->Enable (ADVANCED_ENABLED && ! THREAD_ACTIVE);
}

void CPrime95Doc::OnTime() 
{
	CTimeDlg dlg;

	dlg.m_p = 10000000;
	dlg.m_iter = 10;
	if (dlg.DoModal () == IDOK) {
		thread_pkt.op = OP_TIME;
		thread_pkt.p = dlg.m_p;
		thread_pkt.p2 = dlg.m_iter;
		thread_pkt.doc = this;
		CWinThread *thread = AfxBeginThread (threadDispatch, NULL);
	}
}

void CPrime95Doc::OnUpdatePminus1(CCmdUI* pCmdUI) 
{
	pCmdUI->Enable (ADVANCED_ENABLED && ! THREAD_ACTIVE);
}

void CPrime95Doc::OnPminus1() 
{
	CPminus1Dlg dlg;

	dlg.m_bound1 = 1000000;
	if (dlg.DoModal () == IDOK) {
		work_unit w;
		w.work_type = WORK_PMINUS1;
		w.k = 1.0;
		w.b = 2;
		w.n = dlg.m_p;
		w.c = dlg.m_plus1 ? 1 : -1;
		w.B1 = dlg.m_bound1;
		w.B2_start = 0;
		w.B2_end = dlg.m_bound2;
		addWorkToDoLine (&w);
		OnContinue ();
	}
}

void CPrime95Doc::OnUpdateEcm(CCmdUI* pCmdUI) 
{
	pCmdUI->Enable (ADVANCED_ENABLED && ! THREAD_ACTIVE);
}

void CPrime95Doc::OnEcm() 
{
	CEcmDlg dlg;

	dlg.m_bound1 = 1000000;
	dlg.m_num_curves = 100;
	if (dlg.DoModal () == IDOK) {
		work_unit w;
		w.work_type = WORK_ECM;
		w.k = 1.0;
		w.b = 2;
		w.n = dlg.m_p;
		w.c = dlg.m_plus1 ? 1 : -1;
		w.B1 = dlg.m_bound1;
		w.B2_start = 0;
		w.B2_end = dlg.m_bound2;
		w.curves_to_do = dlg.m_num_curves;
		w.curve = dlg.m_curve;
		addWorkToDoLine (&w);
		OnContinue ();
	}
}

void CPrime95Doc::OnUpdateErrchk(CCmdUI* pCmdUI) 
{
	pCmdUI->SetCheck (ERRCHK);
	pCmdUI->Enable (ADVANCED_ENABLED);
}

void CPrime95Doc::OnErrchk() 
{
	ERRCHK = !ERRCHK;
	IniWriteInt (INI_FILE, "ErrorCheck", ERRCHK);
}

void CPrime95Doc::OnUpdatePriority(CCmdUI* pCmdUI) 
{
	pCmdUI->Enable (ADVANCED_ENABLED);
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
		IniWriteInt (LOCALINI_FILE, "Affinity", CPU_AFFINITY);
	}
}

void CPrime95Doc::OnUpdateManualcomm(CCmdUI* pCmdUI) 
{
	pCmdUI->Enable (ADVANCED_ENABLED && USE_PRIMENET);
}

void CPrime95Doc::OnManualcomm() 
{
	CManualCommDlg dlg;

	dlg.m_manual_comm = MANUAL_COMM & 0x1;
	dlg.m_comm_now = 1;
	dlg.m_new_dates = 0;
	if (dlg.DoModal () == IDOK) {
		MANUAL_COMM = dlg.m_manual_comm;
		IniWriteInt (INI_FILE, "ManualComm", MANUAL_COMM);
		if (dlg.m_new_dates) UpdateEndDates ();
		if (dlg.m_comm_now) {
			MANUAL_COMM |= 0x2;
			CHECK_WORK_QUEUE = 1;
			if (!THREAD_ACTIVE) OnContinue ();
		}
	}
}

void CPrime95Doc::OnUpdateUnreserve(CCmdUI* pCmdUI) 
{
	pCmdUI->Enable (ADVANCED_ENABLED && USE_PRIMENET);
}

void CPrime95Doc::OnUnreserve() 
{
	CUnreserveDlg dlg;

	if (dlg.DoModal () == IDOK) {
		if (THREAD_ACTIVE) {
			Restart1 ();
			unreserve (dlg.m_p);
			Restart2 ();
		} else {
			unreserve (dlg.m_p);
			OnContinue ();
		}
	}
}


// Options menu

void CPrime95Doc::OnCpu() 
{
	CCpuDlg dlg;
	char	buf[512];

	dlg.m_hours = CPU_HOURS;
	dlg.m_day_memory = DAY_MEMORY;
	dlg.m_night_memory = NIGHT_MEMORY;
	minutesToStr (DAY_START_TIME, buf);
	dlg.m_start_time = buf;
	minutesToStr (DAY_END_TIME, buf);
	dlg.m_end_time = buf;
	getCpuDescription (buf, 0);
	dlg.m_cpu_info = buf;
again:	if (dlg.DoModal () == IDOK) {
		unsigned int new_day_start_time, new_day_end_time;

		if (CPU_HOURS != dlg.m_hours) {
			ROLLING_AVERAGE = 1000;
			IniWriteInt (LOCALINI_FILE, "RollingAverage", 1000);
			IniWriteInt (LOCALINI_FILE, "RollingStartTime", 0);
			spoolMessage (PRIMENET_SET_COMPUTER_INFO, NULL);
			next_comm_time = 0;
			UpdateEndDates ();
		}
		new_day_start_time = strToMinutes ((char *)(LPCTSTR) dlg.m_start_time);
		new_day_end_time = strToMinutes ((char *)(LPCTSTR) dlg.m_end_time);
		if (DAY_MEMORY != dlg.m_day_memory ||
		    NIGHT_MEMORY != dlg.m_night_memory ||
		    DAY_START_TIME != new_day_start_time ||
		    DAY_END_TIME != new_day_end_time)
			if (THREAD_ACTIVE) memSettingsChanged ();

/* Save the new information */

		CPU_HOURS = dlg.m_hours;
		DAY_MEMORY = dlg.m_day_memory;
		NIGHT_MEMORY = dlg.m_night_memory;
		DAY_START_TIME = new_day_start_time;
		DAY_END_TIME = new_day_end_time;
		IniWriteInt (LOCALINI_FILE, "CPUHours", CPU_HOURS);
		IniWriteInt (LOCALINI_FILE, "DayMemory", DAY_MEMORY);
		IniWriteInt (LOCALINI_FILE, "NightMemory", NIGHT_MEMORY);
		IniWriteInt (LOCALINI_FILE, "DayStartTime", DAY_START_TIME);
		IniWriteInt (LOCALINI_FILE, "DayEndTime", DAY_END_TIME);

		if (!IniGetInt (INI_FILE, "AskedAboutMemory", 0)) {
			IniWriteInt (INI_FILE, "AskedAboutMemory", 1);
			if (DAY_MEMORY == 8 && NIGHT_MEMORY == 8 &&
			    AfxMessageBox (MSG_MEMORY, MB_YESNO | MB_ICONQUESTION) == IDYES)
				goto again;
		}

		if (STARTUP_IN_PROGRESS) OnPrimenet ();
	} else
		STARTUP_IN_PROGRESS = 0;
}

void CPrime95Doc::OnPreferences() 
{
	CPreferencesDlg dlg;

	dlg.m_iter = ITER_OUTPUT;
	dlg.m_r_iter = ITER_OUTPUT_RES;
	dlg.m_disk_write_time = DISK_WRITE_TIME;
	dlg.m_modem = MODEM_RETRY_TIME;
	dlg.m_retry = NETWORK_RETRY_TIME;
	dlg.m_end_dates = DAYS_BETWEEN_CHECKINS;
	dlg.m_backup = TWO_BACKUP_FILES;
	dlg.m_noise = !SILENT_VICTORY;
	dlg.m_battery = RUN_ON_BATTERY;
	if (dlg.DoModal () == IDOK) {
		ITER_OUTPUT = dlg.m_iter;
		ITER_OUTPUT_RES = dlg.m_r_iter;
		DISK_WRITE_TIME = dlg.m_disk_write_time;
		MODEM_RETRY_TIME = dlg.m_modem;
		NETWORK_RETRY_TIME = dlg.m_retry;
		DAYS_BETWEEN_CHECKINS = dlg.m_end_dates;
		TWO_BACKUP_FILES = dlg.m_backup;
		SILENT_VICTORY = !dlg.m_noise;
		RUN_ON_BATTERY = dlg.m_battery;
		IniWriteInt (INI_FILE, "OutputIterations", ITER_OUTPUT);
		IniWriteInt (INI_FILE, "ResultsFileIterations", ITER_OUTPUT_RES);
		IniWriteInt (INI_FILE, "DiskWriteTime", DISK_WRITE_TIME);
		IniWriteInt (INI_FILE, "NetworkRetryTime", MODEM_RETRY_TIME);
		IniWriteInt (INI_FILE, "NetworkRetryTime2", NETWORK_RETRY_TIME);
		IniWriteInt (INI_FILE, "DaysBetweenCheckins", DAYS_BETWEEN_CHECKINS);
		IniWriteInt (INI_FILE, "TwoBackupFiles", TWO_BACKUP_FILES);
		IniWriteInt (INI_FILE, "SilentVictory", SILENT_VICTORY);
		IniWriteInt (LOCALINI_FILE, "RunOnBattery", RUN_ON_BATTERY);
		CHECK_BATTERY = 1;
	}
}

void CPrime95Doc::OnUpdateBenchmark(CCmdUI* pCmdUI) 
{
	pCmdUI->Enable (! THREAD_ACTIVE);
}

void CPrime95Doc::OnBenchmark() 
{
	CWinThread *thread;

	thread_pkt.op = OP_BENCH;
	thread_pkt.doc = this;
	thread = AfxBeginThread (threadDispatch, NULL);
}

void CPrime95Doc::OnUpdateTorture(CCmdUI* pCmdUI)
{
	pCmdUI->Enable (! THREAD_ACTIVE);
}

void CPrime95Doc::OnTorture() 
{
	CTortureDlg dlg;
	int	mem;

	dlg.m_minfft = 8;
	dlg.m_maxfft = 4096;
	mem = physical_memory ();
	if (mem >= 500) {
		dlg.m_memory = GetSuggestedMemory (mem - 256);
		dlg.m_in_place_fft = FALSE;
	} else if (mem >= 200) {
		dlg.m_memory = GetSuggestedMemory (mem / 2);
		dlg.m_in_place_fft = TRUE;
	} else {
		dlg.m_memory = 8;
		dlg.m_in_place_fft = TRUE;
	}
	dlg.m_timefft = 15;
	if (dlg.DoModal () == IDOK) {
		CWinThread *thread;

		IniWriteInt (INI_FILE, "MinTortureFFT", dlg.m_minfft);
		IniWriteInt (INI_FILE, "MaxTortureFFT", dlg.m_maxfft);
		IniWriteInt (INI_FILE, "TortureMem",
			dlg.m_in_place_fft ? 8 : dlg.m_memory);
		IniWriteInt (INI_FILE, "TortureTime", dlg.m_timefft);

		thread_pkt.op = OP_TORTURE;
		thread_pkt.doc = this;
		thread = AfxBeginThread (threadDispatch, NULL);
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
	pCmdUI->Enable (ADVANCED_ENABLED);
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
// an NT service has turned this option off.  Also, some NT users do not
// have permission to create and delete services.  For those users, change
// the menu text to "Start at logon."

void CPrime95Doc::OnUpdateService(CCmdUI* pCmdUI)
{
	pCmdUI->SetText (canModifyServices () ?
				"Start at Bootup" : "Start at Logon");
	pCmdUI->Enable (!NTSERVICENAME[0] || WINDOWS95_SERVICE);
	pCmdUI->SetCheck (WINDOWS95_SERVICE);
}

void CPrime95Doc::OnService() 
{
	WINDOWS95_SERVICE = !WINDOWS95_SERVICE;
	IniWriteInt (INI_FILE, "Windows95Service", WINDOWS95_SERVICE);
	Service95 ();
}


// Help menu


void CPrime95Doc::OnUpdateHelpFinder(CCmdUI* pCmdUI) 
{
	pCmdUI->Enable (TRUE);
}

void CPrime95Doc::OnUpdateServer(CCmdUI* pCmdUI) 
{
	pCmdUI->Enable (USE_PRIMENET);
}

void CPrime95Doc::OnServer() 
{
	CServerDlg dlg;
	struct primenetPingServerInfo pkt;
	
	memset (&pkt, 0, sizeof (pkt));
	pkt.u.serverInfo.versionNumber = PRIMENET_VERSION;
	if (sendMessage (PRIMENET_PING_SERVER_INFO, &pkt)) {
		AfxMessageBox (PING_ERROR, MB_OK | MB_ICONEXCLAMATION);
	} else {
		dlg.m_line1 = pkt.u.serverInfo.buildID;
		dlg.m_line2 = pkt.u.serverInfo.primenetServerName;
		dlg.m_line3 = pkt.u.serverInfo.adminEmailAddr;
		dlg.DoModal ();
	}
}

// Other internal commands

void CPrime95Doc::OnWelcome() 
{
	CWelcomeDlg dlg;

// After the welcome screen, install prime95 as an auto-start program
// and then go collect the user information.

	if (dlg.DoModal () == IDOK) {
		if (!WINDOWS95_SERVICE) OnService ();
		OnRangeUserinformation();
	} else {
		IniWriteInt (INI_FILE, "StressTester", 1);
		IniWriteInt (INI_FILE, "UsePrimenet", USE_PRIMENET = 0);
		STARTUP_IN_PROGRESS = 0;
	}
}

void CPrime95Doc::OnUsrTorture() 
{
	CWinThread *thread;

	thread_pkt.op = OP_TORTURE;
	thread_pkt.doc = this;
	thread = AfxBeginThread (threadDispatch, NULL);
}

void CPrime95Doc::OnBroadcast() 
{
static	int	msg_box_up = 0;
	char	filename[33];

	if (msg_box_up || BROADCAST_MESSAGE == NULL) return;
	msg_box_up = 1;
	AfxGetApp()->m_pMainWnd->MessageBox (
		BROADCAST_MESSAGE, "Important Message from PrimeNet Server", MB_OK | MB_ICONEXCLAMATION);
	free (BROADCAST_MESSAGE);
	BROADCAST_MESSAGE = NULL;
	strcpy (filename, "bcastmsg");
	strcat (filename, EXTENSION);
	unlink (filename);
	BlinkIcon (-1);
	msg_box_up = 0;
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

	p = lines[NumLines-1];
	memmove (&lines[1], &lines[0], (NumLines-1)*sizeof(char *));
	lines[0] = p;
	*p = 0;
	UpdateAllViews (NULL);
}

void OutputStr (char *str)
{
	OUTPUT_STR_HACK->OutputStr (str);
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
	if (EXIT_IN_PROGRESS) return;
	if (lines[0] == NULL) {
		for (int i = 0; i < NumLines; i++) {
			lines[i] = (char *) malloc (128);
			*lines[i] = 0;
		}
	}
	char *p = lines[0] + strlen (lines[0]);
	for ( ; *str; str++) {
		if (*str == '\r') continue;
		if (*str == '\n') *p = 0, LineFeed (), p = lines[0];
		else if (p - lines[0] < 127) *p++ = *str;
	}
	*p = 0;
}


/////////////////////////////////////////////////////////////////////////////
// CPrime95Doc public routines

#include <ctype.h>
#include <dos.h>
#include <fcntl.h>
#include <io.h>
#include <memory.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <sys/timeb.h>
#include <gwutil.h>

#ifdef X86_64
#define PORT	4
#else
#define PORT	1
#endif
#include "commona.c"
#include "commonb.c"
#include "commonc.c"
#include "ecm.c"
#include "comm95a.c"
#include "comm95b.c"
#include "comm95c.c"
#include "primenet.c"

UINT threadDispatch (
	LPVOID stuff)
{

// Set thread-active flags

	THREAD_ACTIVE = TRUE;
	THREAD_STOP = 0;

// Stall if we've just booted (within 5 minutes of Windows starting)

	if (GetTickCount () < 300000 && thread_pkt.op == OP_CONTINUE) {
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
	case OP_TIME:
		primeTime (thread_pkt.p, thread_pkt.p2);
		break;
	case OP_BENCH:
		primeBench ();
		break;
	case OP_TORTURE:
		selfTest (1);
		break;
	}

// Output informative message

	if (!GIMPS_QUIT && !STOPPED_ON_BATTERY && !RESTARTING) {
		if (thread_pkt.op == OP_CONTINUE || thread_pkt.op == OP_TORTURE)
			OutputStr ("Execution halted.\n");
		if (thread_pkt.op == OP_CONTINUE)
			OutputStr ("Choose Test/Continue to restart.\n");
	}

// Clear thread-active flags

	ChangeIcon (IDLE_ICON);
	THREAD_ACTIVE = FALSE;
	THREAD_STOP = 0;
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
		    (power.ACLineStatus != 1 || 
		     (power.ACLineStatus == 1 &&
		      power.BatteryLifePercent < IniGetInt (INI_FILE, "BatteryPercent", 0)))) {
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

	sprintf (buf, "Prime95 - %s", str);
	if (pApp->m_pMainWnd) {
		if (TRAY_ICON)
			((CPrime95App *)pApp)->
				TrayMessage (NIM_MODIFY, buf, 0);
		if (pApp->m_pMainWnd->IsIconic ()) {
			pApp->m_pMainWnd->SetWindowText (buf);
			was_iconic = TRUE;
		} else if (was_iconic) {
			pApp->m_pMainWnd->SetWindowText ("Prime95");
			was_iconic = FALSE;
		}
	}
}

// Output a status report for the range

void CPrime95Doc::rangeStatus ()
{
	char	buf[2000];

	rangeStatusMessage (buf);

	AfxMessageBox (buf, MB_ICONINFORMATION);
}
