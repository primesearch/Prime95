// ContentsDlg.cpp : implementation file
//

#include "stdafx.h"
#include "Prime95.h"
#include "ContentsDlg.h"

#ifdef _DEBUG
#define new DEBUG_NEW
#undef THIS_FILE
static char THIS_FILE[] = __FILE__;
#endif

/////////////////////////////////////////////////////////////////////////////
// ContentsDlg dialog


ContentsDlg::ContentsDlg(CWnd* pParent /*=NULL*/)
	: CDialog(ContentsDlg::IDD, pParent)
{
	//{{AFX_DATA_INIT(ContentsDlg)
	m_helptext = _T("");
	//}}AFX_DATA_INIT
}


void ContentsDlg::DoDataExchange(CDataExchange* pDX)
{
	CDialog::DoDataExchange(pDX);
	//{{AFX_DATA_MAP(ContentsDlg)
	DDX_Text(pDX, IDC_HELP_TEXT, m_helptext);
	//}}AFX_DATA_MAP
}


BEGIN_MESSAGE_MAP(ContentsDlg, CDialog)
	//{{AFX_MSG_MAP(ContentsDlg)
		// NOTE: the ClassWizard will add message map macros here
	//}}AFX_MSG_MAP
END_MESSAGE_MAP()

/////////////////////////////////////////////////////////////////////////////
// ContentsDlg message handlers
