// EcmDlg.cpp : implementation file
//

#include "stdafx.h"
#include "Prime95.h"
#include "EcmDlg.h"

#ifdef _DEBUG
#define new DEBUG_NEW
#undef THIS_FILE
static char THIS_FILE[] = __FILE__;
#endif

/////////////////////////////////////////////////////////////////////////////
// CEcmDlg dialog


CEcmDlg::CEcmDlg(CWnd* pParent /*=NULL*/)
	: CDialog(CEcmDlg::IDD, pParent)
{
	//{{AFX_DATA_INIT(CEcmDlg)
	m_p = 0;
	m_bound1 = 0;
	m_bound2 = 0;
	m_curve = 0.0;
	m_num_curves = 0;
	m_plus1 = FALSE;
	//}}AFX_DATA_INIT
}


void CEcmDlg::DoDataExchange(CDataExchange* pDX)
{
	CDialog::DoDataExchange(pDX);
	//{{AFX_DATA_MAP(CEcmDlg)
	DDX_Text(pDX, IDC_P1, m_p);
	DDV_MinMaxUInt(pDX, m_p, 100, 79300000);
	DDX_Text(pDX, IDC_P2, m_bound1);
	DDV_MinMaxUInt(pDX, m_bound1, 0, 4000000000);
	DDX_Text(pDX, IDC_P3, m_bound2);
	DDV_MinMaxUInt(pDX, m_bound2, 0, 4000000000);
	DDX_Text(pDX, IDC_CURVE, m_curve);
	DDX_Text(pDX, IDC_NUM_CURVES, m_num_curves);
	DDV_MinMaxUInt(pDX, m_num_curves, 1, 999999999);
	DDX_Check(pDX, IDC_PLUS1, m_plus1);
	//}}AFX_DATA_MAP
}


BEGIN_MESSAGE_MAP(CEcmDlg, CDialog)
	//{{AFX_MSG_MAP(CEcmDlg)
		// NOTE: the ClassWizard will add message map macros here
	//}}AFX_MSG_MAP
END_MESSAGE_MAP()

/////////////////////////////////////////////////////////////////////////////
// CEcmDlg message handlers
