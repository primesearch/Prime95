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
	m_cpu_type = -1;
	//}}AFX_DATA_INIT
}


void CCpuDlg::DoDataExchange(CDataExchange* pDX)
{
	CDialog::DoDataExchange(pDX);
	//{{AFX_DATA_MAP(CCpuDlg)
	DDX_Control(pDX, IDC_P4, c_p4);
	DDX_Control(pDX, IDC_PRO, c_pentium_pro);
	DDX_Control(pDX, IDC_PENTIUM, c_pentium);
	DDX_Control(pDX, IDC_K6, c_k6);
	DDX_Control(pDX, IDC_CELERON, c_celeron);
	DDX_Control(pDX, IDC_PII, c_pii);
	DDX_Control(pDX, IDC_PIII, c_piii);
	DDX_Control(pDX, IDC_K7, c_k7);
	DDX_Text(pDX, IDC_P, m_speed);
	DDV_MinMaxUInt(pDX, m_speed, 25, 10000);
	DDX_Radio(pDX, IDC_P4, m_cpu_type);
	//}}AFX_DATA_MAP
	c_k6.EnableWindow (isPentium ());
	c_k7.EnableWindow (isPentium ());
	c_pentium.EnableWindow (isPentium ());
	c_pentium_pro.EnableWindow (isPentiumPro ());
	c_celeron.EnableWindow (isPentiumPro ());
	c_pii.EnableWindow (isPentiumPro ());
	c_piii.EnableWindow (isPentium3 ());
	c_p4.EnableWindow (isPentium4 ());
}


BEGIN_MESSAGE_MAP(CCpuDlg, CDialog)
	//{{AFX_MSG_MAP(CCpuDlg)
		// NOTE: the ClassWizard will add message map macros here
	//}}AFX_MSG_MAP
END_MESSAGE_MAP()

/////////////////////////////////////////////////////////////////////////////
// CCpuDlg message handlers
