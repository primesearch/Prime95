// ServerDlg.cpp : implementation file
//

#include "stdafx.h"
#include "Prime95.h"
#include "ServerDlg.h"

#ifdef _DEBUG
#define new DEBUG_NEW
#undef THIS_FILE
static char THIS_FILE[] = __FILE__;
#endif

/////////////////////////////////////////////////////////////////////////////
// CServerDlg dialog


CServerDlg::CServerDlg(CWnd* pParent /*=NULL*/)
	: CDialog(CServerDlg::IDD, pParent)
{
	//{{AFX_DATA_INIT(CServerDlg)
	m_line1 = _T("");
	m_line2 = _T("");
	m_line3 = _T("");
	//}}AFX_DATA_INIT
}


void CServerDlg::DoDataExchange(CDataExchange* pDX)
{
	CDialog::DoDataExchange(pDX);
	//{{AFX_DATA_MAP(CServerDlg)
	DDX_Text(pDX, IDC_LINE1, m_line1);
	DDX_Text(pDX, IDC_LINE2, m_line2);
	DDX_Text(pDX, IDC_LINE3, m_line3);
	//}}AFX_DATA_MAP
}


BEGIN_MESSAGE_MAP(CServerDlg, CDialog)
	//{{AFX_MSG_MAP(CServerDlg)
		// NOTE: the ClassWizard will add message map macros here
	//}}AFX_MSG_MAP
END_MESSAGE_MAP()

/////////////////////////////////////////////////////////////////////////////
// CServerDlg message handlers
