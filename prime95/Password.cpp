// Password.cpp : implementation file
//

#include "stdafx.h"
#include "Prime95.h"
#include "Password.h"

#ifdef _DEBUG
#define new DEBUG_NEW
#undef THIS_FILE
static char THIS_FILE[] = __FILE__;
#endif

/////////////////////////////////////////////////////////////////////////////
// Password dialog


Password::Password(CWnd* pParent /*=NULL*/)
	: CDialog(Password::IDD, pParent)
{
	//{{AFX_DATA_INIT(Password)
	m_password = 0;
	//}}AFX_DATA_INIT
}


void Password::DoDataExchange(CDataExchange* pDX)
{
	CDialog::DoDataExchange(pDX);
	//{{AFX_DATA_MAP(Password)
	DDX_Text(pDX, IDC_P, m_password);
	//}}AFX_DATA_MAP
}


BEGIN_MESSAGE_MAP(Password, CDialog)
	//{{AFX_MSG_MAP(Password)
		// NOTE: the ClassWizard will add message map macros here
	//}}AFX_MSG_MAP
END_MESSAGE_MAP()

/////////////////////////////////////////////////////////////////////////////
// Password message handlers
