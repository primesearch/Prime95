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
	m_hours = 0;
	m_start_time = _T("");
	m_end_time = _T("");
	m_day_memory = 0;
	m_night_memory = 0;
	m_cpu_info = _T("");
	//}}AFX_DATA_INIT
}


void CCpuDlg::DoDataExchange(CDataExchange* pDX)
{
	CDialog::DoDataExchange(pDX);
	//{{AFX_DATA_MAP(CCpuDlg)
	DDX_Text(pDX, IDC_HOURS, m_hours);
	DDV_MinMaxUInt(pDX, m_hours, 1, 24);
	DDX_Text(pDX, IDC_START_TIME, m_start_time);
	DDX_Text(pDX, IDC_END_TIME, m_end_time);
	DDX_Text(pDX, IDC_DAY_MEMORY, m_day_memory);
	DDV_MinMaxUInt(pDX, m_day_memory, 8, (UINT) (0.9 * physical_memory ()));
	DDX_Text(pDX, IDC_NIGHT_MEMORY, m_night_memory);
	DDV_MinMaxUInt(pDX, m_night_memory, 8, (UINT) (0.9 * physical_memory ()));
	DDX_Text(pDX, IDC_CPU_INFO, m_cpu_info);
	//}}AFX_DATA_MAP
}


BEGIN_MESSAGE_MAP(CCpuDlg, CDialog)
	//{{AFX_MSG_MAP(CCpuDlg)
	//}}AFX_MSG_MAP
END_MESSAGE_MAP()

/////////////////////////////////////////////////////////////////////////////
// CCpuDlg message handlers

