// PreferencesDlg.cpp : implementation file
//

#include "stdafx.h"
#include "Prime95.h"
#include "PreferencesDlg.h"

#ifdef _DEBUG
#define new DEBUG_NEW
#undef THIS_FILE
static char THIS_FILE[] = __FILE__;
#endif

/////////////////////////////////////////////////////////////////////////////
// CPreferencesDlg dialog


CPreferencesDlg::CPreferencesDlg(CWnd* pParent /*=NULL*/)
	: CDialog(CPreferencesDlg::IDD, pParent)
{
	//{{AFX_DATA_INIT(CPreferencesDlg)
	m_iter = 0;
	m_disk_write_time = 0;
	m_backup = FALSE;
	m_r_iter = 0;
	m_battery = FALSE;
	m_cumulative = FALSE;
	//}}AFX_DATA_INIT
}


void CPreferencesDlg::DoDataExchange(CDataExchange* pDX)
{
	CDialog::DoDataExchange(pDX);
	//{{AFX_DATA_MAP(CPreferencesDlg)
	DDX_Text(pDX, IDC_P, m_iter);
	DDV_MinMaxUInt(pDX, m_iter, 1, 999999999);
	DDX_Text(pDX, IDC_DISK, m_disk_write_time);
	DDV_MinMaxUInt(pDX, m_disk_write_time, 10, 999999);
	DDX_Check(pDX, IDC_BACKUP, m_backup);
	DDX_Text(pDX, IDC_R_ITER, m_r_iter);
	DDV_MinMaxUInt(pDX, m_r_iter, 10000, 999999999);
	DDX_Check(pDX, IDC_BATTERY, m_battery);
	DDX_Check(pDX, IDC_CUMULATIVE, m_cumulative);
	//}}AFX_DATA_MAP
}


BEGIN_MESSAGE_MAP(CPreferencesDlg, CDialog)
	//{{AFX_MSG_MAP(CPreferencesDlg)
		// NOTE: the ClassWizard will add message map macros here
	//}}AFX_MSG_MAP
END_MESSAGE_MAP()

/////////////////////////////////////////////////////////////////////////////
// CPreferencesDlg message handlers
