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
#include "Pminus1Dlg.h"
#include "PreferencesDlg.h"
#include "PrimeNetDlg.h"
#include "Priority.h"
#include "ServerDlg.h"
#include "UserDlg.h"
#include "VacationDlg.h"

#ifdef _DEBUG
#define new DEBUG_NEW
#undef THIS_FILE
static char THIS_FILE[] = __FILE__;
#endif

/////////////////////////////////////////////////////////////////////////////
// CPrime95Doc

IMPLEMENT_DYNCREATE(CPrime95Doc, CDocument)

BEGIN_MESSAGE_MAP(CPrime95Doc, CDocument)
	//{{AFX_MSG_MAP(CPrime95Doc)
	ON_COMMAND(IDM_CPU, OnCpu)
	ON_COMMAND(IDM_PREFERENCES, OnPreferences)
	ON_COMMAND(ID_RANGE_STATUS, OnRangeStatus)
	ON_UPDATE_COMMAND_UI(ID_HELP_FINDER, OnUpdateHelpFinder)
	ON_COMMAND(ID_RANGE_USERINFORMATION, OnRangeUserinformation)
	ON_COMMAND(IDM_PRIMENET, OnPrimenet)
	ON_COMMAND(IDM_PRIORITY, OnPriority)
	ON_COMMAND(IDM_SERVER, OnServer)
	ON_UPDATE_COMMAND_UI(IDM_SERVER, OnUpdateServer)
	ON_COMMAND(IDM_QUIT, OnQuitGimps)
	ON_COMMAND(IDM_VACATION, OnVacation)
	ON_COMMAND(IDM_MANUALCOMM, OnManualcomm)
	ON_UPDATE_COMMAND_UI(IDM_MANUALCOMM, OnUpdateManualcomm)
	ON_COMMAND(IDM_AFFINITY, OnAffinity)
	ON_COMMAND(IDM_ECM, OnEcm)
	ON_UPDATE_COMMAND_UI(IDM_ECM, OnUpdateEcm)
	ON_COMMAND(IDM_PMINUS1, OnPminus1)
	ON_UPDATE_COMMAND_UI(IDM_PMINUS1, OnUpdatePminus1)
	ON_COMMAND(IDM_BROADCAST, OnBroadcast)
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

// Remember the main window's size and position

	CWinApp* pApp = AfxGetApp();

// Free the networking library

	if (HLIB) FreeLibrary (HLIB);

// Finish closing

	CDocument::OnCloseDocument();
}

// Test menu

void CPrime95Doc::OnPrimenet() 
{
	PrimenetDlg dlg;
	short	work_pref;

	dlg.m_primenet = USE_PRIMENET;
	dlg.m_dialup = DIAL_UP;
	dlg.m_rpc = USE_HTTP;
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
	dlg.m_factor = !! (work_pref & PRIMENET_ASSIGN_FACTOR);
	dlg.m_dblchk = !! (work_pref & PRIMENET_ASSIGN_DBLCHK);

	if (dlg.DoModal () == IDOK) {
		if (!USE_PRIMENET && dlg.m_primenet) {
			USE_PRIMENET = 1;
			spoolMessage (PRIMENET_MAINTAIN_USER_INFO, NULL);
			spoolMessage (PRIMENET_SET_COMPUTER_INFO, NULL);
			spoolExistingResultsFile ();
		}
		if (dlg.m_rpc != USE_HTTP && HLIB) {
			FreeLibrary (HLIB);
			HLIB = 0;
		}
		USE_PRIMENET = dlg.m_primenet;
		DIAL_UP = dlg.m_dialup;
		USE_HTTP = dlg.m_rpc;
		DAYS_OF_WORK = dlg.m_work;
		if (dlg.m_work_dflt)
			WORK_PREFERENCE = 0;
		else
			WORK_PREFERENCE =
				(dlg.m_bigones ? PRIMENET_ASSIGN_BIGONES : 0) +
				(dlg.m_lucas ? PRIMENET_ASSIGN_TEST : 0) +
				(dlg.m_factor ? PRIMENET_ASSIGN_FACTOR : 0) +
				(dlg.m_dblchk ? PRIMENET_ASSIGN_DBLCHK : 0);
		IniWriteInt (INI_FILE, "UsePrimenet", USE_PRIMENET);
		IniWriteInt (INI_FILE, "DialUp", DIAL_UP);
		IniWriteInt (INI_FILE, "UseHTTP", USE_HTTP);
		IniWriteInt (INI_FILE, "DaysOfWork", DAYS_OF_WORK);
		IniWriteInt (INI_FILE, "WorkPreference", WORK_PREFERENCE);
	}
	STARTUP_IN_PROGRESS = 0;
}

void CPrime95Doc::OnUpdatePminus1(CCmdUI* pCmdUI) 
{
	pCmdUI->Enable (! USE_PRIMENET);
}

void CPrime95Doc::OnPminus1() 
{
	CPminus1Dlg dlg;

	dlg.m_bound1 = 1000000;
	if (dlg.DoModal () == IDOK) {
		work_unit w;
		w.work_type = WORK_PMINUS1;
		w.p = dlg.m_p;
		w.B1 = dlg.m_bound1;
		w.B2_start = 0;
		w.B2_end = dlg.m_bound2;
		w.plus1 = dlg.m_plus1;
		addWorkToDoLine (&w);
	}
}

void CPrime95Doc::OnUpdateEcm(CCmdUI* pCmdUI) 
{
	pCmdUI->Enable (! USE_PRIMENET);
}

void CPrime95Doc::OnEcm() 
{
        CEcmDlg dlg;

        dlg.m_bound1 = 1000000;
        dlg.m_num_curves = 100;
        if (dlg.DoModal () == IDOK) {
                work_unit w;
                w.work_type = WORK_ECM;
                w.p = dlg.m_p;
                w.B1 = dlg.m_bound1;
                w.B2_start = 0;
                w.B2_end = dlg.m_bound2;
                w.curves_to_do = dlg.m_num_curves;
                w.curves_completed = 0;
                w.curve = dlg.m_curve;
                w.plus1 = dlg.m_plus1;
                addWorkToDoLine (&w);
        }
}

void CPrime95Doc::OnQuitGimps() 
{
	if (!USE_PRIMENET) {
		if (AfxMessageBox (MANUAL_QUIT, MB_YESNO | MB_ICONQUESTION | MB_DEFBUTTON2) == IDYES) {
			writeResults ("Quitting GIMPS.\n");
			IniDeleteAllLines (WORKTODO_FILE);
		}
	} else {
		if (AfxMessageBox (PRIMENET_QUIT, MB_YESNO | MB_ICONQUESTION | MB_DEFBUTTON2) == IDYES) {
			writeResults ("Quitting GIMPS.\n");
			spoolMessage (999, NULL);
			MANUAL_COMM |= 0x2;
			CHECK_WORK_QUEUE = 1;
			communicateWithServer ();
		}
	}
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
                IniWriteInt (LOCALINI_FILE, "VacationEnd", VACATION_END);
                IniWriteInt (LOCALINI_FILE, "VacationOn", ON_DURING_VACATION);
		if (VACATION_END && !ON_DURING_VACATION)
			IniWriteInt (LOCALINI_FILE, "RollingStartTime", 0);
		UpdateEndDates ();
		MANUAL_COMM |= 0x2;
		CHECK_WORK_QUEUE = 1;
		communicateWithServer ();
	}
}

void CPrime95Doc::OnRangeStatus() 
{
	rangeStatus ();
}


// Options menu

void CPrime95Doc::OnCpu() 
{
	CCpuDlg dlg;
	char	buf[20];

	dlg.m_speed = CPU_SPEED;
	dlg.m_cpu_type =
		(CPU_TYPE == 10) ? 0 : (CPU_TYPE == 9) ? 1 :
		(CPU_TYPE == 8) ? 2 : (CPU_TYPE == 6) ? 3 :
		(CPU_TYPE == 5) ? 4 : (CPU_TYPE == 4) ? 5 :
		(CPU_TYPE == 11) ? 6 : (CPU_TYPE == 7) ? 7 : 8;
	dlg.m_hours = CPU_HOURS;
	dlg.m_day_memory = DAY_MEMORY;
	dlg.m_night_memory = NIGHT_MEMORY;
	minutesToStr (DAY_START_TIME, buf);
	dlg.m_start_time = buf;
	minutesToStr (DAY_END_TIME, buf);
	dlg.m_end_time = buf;
again:	if (dlg.DoModal () == IDOK) {
		unsigned int new_cpu_type;

		new_cpu_type = (dlg.m_cpu_type == 0) ? 10 :
			       (dlg.m_cpu_type == 1) ? 9 :
			       (dlg.m_cpu_type == 2) ? 8 :
			       (dlg.m_cpu_type == 3) ? 6 :
			       (dlg.m_cpu_type == 4) ? 5 :
			       (dlg.m_cpu_type == 5) ? 4 :
			       (dlg.m_cpu_type == 6) ? 11 :
			       (dlg.m_cpu_type == 7) ? 7 : 3;
		if (CPU_SPEED != dlg.m_speed ||
		    CPU_TYPE != new_cpu_type ||
		    CPU_HOURS != dlg.m_hours) {
			ROLLING_AVERAGE = 1000;
			IniWriteInt (LOCALINI_FILE, "RollingAverage", 1000);
			IniWriteInt (LOCALINI_FILE, "RollingStartTime", 0);
			spoolMessage (PRIMENET_SET_COMPUTER_INFO, NULL);
			UpdateEndDates ();
		}
		CPU_SPEED = dlg.m_speed;
		CPU_TYPE = new_cpu_type;
		CPU_HOURS = dlg.m_hours;
		DAY_MEMORY = dlg.m_day_memory;
		NIGHT_MEMORY = dlg.m_night_memory;
		DAY_START_TIME = strToMinutes ((char *)(LPCTSTR) dlg.m_start_time);
		DAY_END_TIME = strToMinutes ((char *)(LPCTSTR) dlg.m_end_time);
		IniWriteInt (LOCALINI_FILE, "CPUType", CPU_TYPE);
		IniWriteInt (LOCALINI_FILE, "CPUSpeed", CPU_SPEED);
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

	dlg.m_r_iter = ITER_OUTPUT_RES;
	dlg.m_disk_write_time = DISK_WRITE_TIME;
	dlg.m_modem = MODEM_RETRY_TIME;
	dlg.m_retry = NETWORK_RETRY_TIME;
	dlg.m_end_dates = DAYS_BETWEEN_CHECKINS;
	dlg.m_backup = TWO_BACKUP_FILES;
	dlg.m_noise = !SILENT_VICTORY;
	if (dlg.DoModal () == IDOK) {
		ITER_OUTPUT_RES = dlg.m_r_iter;
		DISK_WRITE_TIME = dlg.m_disk_write_time;
		MODEM_RETRY_TIME = dlg.m_modem;
		NETWORK_RETRY_TIME = dlg.m_retry;
		DAYS_BETWEEN_CHECKINS = dlg.m_end_dates;
		TWO_BACKUP_FILES = dlg.m_backup;
		SILENT_VICTORY = !dlg.m_noise;
		IniWriteInt (INI_FILE, "ResultsFileIterations", ITER_OUTPUT_RES);
		IniWriteInt (INI_FILE, "DiskWriteTime", DISK_WRITE_TIME);
		IniWriteInt (INI_FILE, "NetworkRetryTime", MODEM_RETRY_TIME);
		IniWriteInt (INI_FILE, "NetworkRetryTime2", NETWORK_RETRY_TIME);
		IniWriteInt (INI_FILE, "DaysBetweenCheckins", DAYS_BETWEEN_CHECKINS);
		IniWriteInt (INI_FILE, "TwoBackupFiles", TWO_BACKUP_FILES);
		IniWriteInt (INI_FILE, "SilentVictory", SILENT_VICTORY);
	}
}

void CPrime95Doc::OnUpdateManualcomm(CCmdUI* pCmdUI)
{
	pCmdUI->Enable (USE_PRIMENET);
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
			communicateWithServer ();
                }
        }
}

void CPrime95Doc::OnPriority() 
{
	CPriority dlg;

	dlg.m_priority = PRIORITY;
	if (dlg.DoModal () == IDOK) {
		if (PRIORITY != dlg.m_priority) {
			PRIORITY = dlg.m_priority;
			IniWriteInt (INI_FILE, "Priority", PRIORITY);
		}
	}
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


// Help menu


void CPrime95Doc::OnUpdateHelpFinder(CCmdUI* pCmdUI) 
{
	pCmdUI->Enable (FALSE);
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

void CPrime95Doc::OutputStr (
	char	*str)
{
	if (lines[0] == NULL) {
		for (int i = 0; i < NumLines; i++) {
			lines[i] = (char *) malloc (128);
			*lines[i] = 0;
		}
	}
	char *p = lines[0] + strlen (lines[0]);
	for ( ; *str; str++) {
		if (*str == '\n') *p = 0, LineFeed (), p = lines[0];
		else *p++ = *str;
	}
	*p = 0;
}

void CPrime95Doc::OutputNum (
	LONG	num)
{
	char	buf[20];

	sprintf (buf, "%ld", num);
	OutputStr (buf);
}

void CPrime95Doc::OutputClocks (
	double	tickCount)
{
	char	buf[80];

        sprintf (buf, "Clocks: %.0f = %.3f sec.\n",
                 tickCount, tickCount / CPU_SPEED / 1000000);
	OutputStr (buf);
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

#include "cpuid.c"
#include "speed.c"

#include "gwnum.c"
#include "commona.c"
#include "commonc.c"
#include "comm95a.c"
#include "comm95c.c"

/* Return TRUE if we should stop calculating */

int escapeCheck ()
{
	return (FALSE);
}

/* Put a title on the main window */

void CPrime95Doc::title (
	char	*str)
{
	char buf[80];
	CWinApp* pApp = AfxGetApp();

	strcpy (buf, "Setup for NT Service version of Prime95");
	pApp->m_pMainWnd->SetWindowText (buf);
}

// Output a status report for the range

void CPrime95Doc::rangeStatus ()
{
	char	buf[2000];

	rangeStatusMessage (buf);

	AfxMessageBox (buf, MB_ICONINFORMATION);
}
