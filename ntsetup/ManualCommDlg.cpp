// ManualCommDlg.cpp : implementation file
//

#include "stdafx.h"
#include "Prime95.h"
#include "ManualCommDlg.h"

#ifdef _DEBUG
#define new DEBUG_NEW
#undef THIS_FILE
static char THIS_FILE[] = __FILE__;
#endif

/////////////////////////////////////////////////////////////////////////////
// CManualCommDlg dialog


CManualCommDlg::CManualCommDlg(CWnd* pParent /*=NULL*/)
	: CDialog(CManualCommDlg::IDD, pParent)
{
	//{{AFX_DATA_INIT(CManualCommDlg)
	m_manual_comm = FALSE;
	m_comm_now = FALSE;
	m_new_dates = FALSE;
	//}}AFX_DATA_INIT
}


void CManualCommDlg::DoDataExchange(CDataExchange* pDX)
{
	CDialog::DoDataExchange(pDX);
	//{{AFX_DATA_MAP(CManualCommDlg)
	DDX_Check(pDX, IDC_MANUAL, m_manual_comm);
	DDX_Check(pDX, IDC_NOW, m_comm_now);
	DDX_Check(pDX, IDC_COMPLETION, m_new_dates);
	//}}AFX_DATA_MAP
}


BEGIN_MESSAGE_MAP(CManualCommDlg, CDialog)
	//{{AFX_MSG_MAP(CManualCommDlg)
		// NOTE: the ClassWizard will add message map macros here
	//}}AFX_MSG_MAP
END_MESSAGE_MAP()

/////////////////////////////////////////////////////////////////////////////
// CManualCommDlg message handlers
