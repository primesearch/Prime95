// CpuDlg.cpp : implementation file
//

#include "stdafx.h"
#include "Prime95.h"
#include "CpuDlg.h"

#ifdef _DEBUG
#define new DEBUG_NEW
#undef THIS_FILE
static char THIS_FILE[] = __FILE__;
#endif

/////////////////////////////////////////////////////////////////////////////
// CCpuDlg dialog


CCpuDlg::CCpuDlg(CWnd* pParent /*=NULL*/)
	: CDialog(CCpuDlg::IDD, pParent)
{
	//{{AFX_DATA_INIT(CCpuDlg)
	m_speed = 0;
	m_hours = 0;
	m_cpu_type = -1;
	m_day_memory = 0;
	m_night_memory = 0;
	m_end_time = _T("");
	m_start_time = _T("");
	//}}AFX_DATA_INIT
}


void CCpuDlg::DoDataExchange(CDataExchange* pDX)
{
	CDialog::DoDataExchange(pDX);
	//{{AFX_DATA_MAP(CCpuDlg)
	DDX_Control(pDX, IDC_PRO, c_pentium_pro);
	DDX_Control(pDX, IDC_CELERON, c_celeron);
	DDX_Control(pDX, IDC_PII, c_pii);
	DDX_Control(pDX, IDC_PIII, c_piii);
	DDX_Control(pDX, IDC_PENTIUM, c_pentium);
	DDX_Control(pDX, IDC_K6, c_k6);
	DDX_Control(pDX, IDC_K7, c_k7);
	DDX_Text(pDX, IDC_P, m_speed);
	DDV_MinMaxUInt(pDX, m_speed, 25, 10000);
	DDX_Text(pDX, IDC_HOURS, m_hours);
	DDV_MinMaxUInt(pDX, m_hours, 1, 24);
	DDX_Radio(pDX, IDC_PIII, m_cpu_type);
	DDX_Text(pDX, IDC_DAY_MEMORY, m_day_memory);
	DDV_MinMaxUInt(pDX, m_day_memory, 8, (UINT) (0.9 * physical_memory ()));
	DDX_Text(pDX, IDC_NIGHT_MEMORY, m_night_memory);
	DDV_MinMaxUInt(pDX, m_night_memory, 8, (UINT) (0.9 * physical_memory ()));
	DDX_Text(pDX, IDC_END_TIME, m_end_time);
	DDX_Text(pDX, IDC_START_TIME, m_start_time);
	//}}AFX_DATA_MAP
	c_k6.EnableWindow (isPentium ());
	c_k7.EnableWindow (isPentium ());
	c_pentium.EnableWindow (isPentium ());
	c_pentium_pro.EnableWindow (isPentiumPro ());
	c_celeron.EnableWindow (isPentiumPro ());
	c_pii.EnableWindow (isPentiumPro ());
	c_piii.EnableWindow (isPentiumPro ());
}


BEGIN_MESSAGE_MAP(CCpuDlg, CDialog)
	//{{AFX_MSG_MAP(CCpuDlg)
		// NOTE: the ClassWizard will add message map macros here
	//}}AFX_MSG_MAP
END_MESSAGE_MAP()

/////////////////////////////////////////////////////////////////////////////
// CCpuDlg message handlers
