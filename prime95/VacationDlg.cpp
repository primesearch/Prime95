// VacationDlg.cpp : implementation file
//

#include "stdafx.h"
#include "Prime95.h"
#include "VacationDlg.h"

#ifdef _DEBUG
#define new DEBUG_NEW
#undef THIS_FILE
static char THIS_FILE[] = __FILE__;
#endif

/////////////////////////////////////////////////////////////////////////////
// CVacationDlg dialog


CVacationDlg::CVacationDlg(CWnd* pParent /*=NULL*/)
	: CDialog(CVacationDlg::IDD, pParent)
{
	//{{AFX_DATA_INIT(CVacationDlg)
	m_computer_on = FALSE;
	m_vacation_days = 0;
	//}}AFX_DATA_INIT
}


void CVacationDlg::DoDataExchange(CDataExchange* pDX)
{
	CDialog::DoDataExchange(pDX);
	//{{AFX_DATA_MAP(CVacationDlg)
	DDX_Check(pDX, IDC_CHECK1, m_computer_on);
	DDX_Text(pDX, IDC_DAYS, m_vacation_days);
	DDV_MinMaxUInt(pDX, m_vacation_days, 0, 120);
	//}}AFX_DATA_MAP
}


BEGIN_MESSAGE_MAP(CVacationDlg, CDialog)
	//{{AFX_MSG_MAP(CVacationDlg)
		// NOTE: the ClassWizard will add message map macros here
	//}}AFX_MSG_MAP
END_MESSAGE_MAP()

/////////////////////////////////////////////////////////////////////////////
// CVacationDlg message handlers
