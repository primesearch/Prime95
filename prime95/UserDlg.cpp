// UserDlg.cpp : implementation file
//

#include "stdafx.h"
#include "Prime95.h"
#include "UserDlg.h"

#ifdef _DEBUG
#define new DEBUG_NEW
#undef THIS_FILE
static char THIS_FILE[] = __FILE__;
#endif

/////////////////////////////////////////////////////////////////////////////
// CUserDlg dialog


CUserDlg::CUserDlg(CWnd* pParent /*=NULL*/)
	: CDialog(CUserDlg::IDD, pParent)
{
	//{{AFX_DATA_INIT(CUserDlg)
	m_email = _T("");
	m_name = _T("");
	m_userid = _T("");
	m_compid = _T("");
	m_password = _T("");
	m_sendemail = FALSE;
	m_team = FALSE;
	//}}AFX_DATA_INIT
}

void CUserDlg::DoDataExchange(CDataExchange* pDX)
{
	CDialog::DoDataExchange(pDX);
	//{{AFX_DATA_MAP(CUserDlg)
	DDX_Control(pDX, IDC_PWD_TEXT, c_pwd_text);
	DDX_Control(pDX, IDC_NAME_TEXT, c_name_text);
	DDX_Control(pDX, IDC_ID_TEXT, c_id_text);
	DDX_Control(pDX, IDC_USERID, c_userid);
	DDX_Control(pDX, IDC_PASSWORD, c_password);
	DDX_Text(pDX, IDC_EMAIL, m_email);
	DDV_MaxChars(pDX, m_email, 76);
	DDX_Text(pDX, IDC_NAME, m_name);
	DDV_MaxChars(pDX, m_name, 76);
	DDX_Text(pDX, IDC_USERID, m_userid);
	DDV_MaxChars(pDX, m_userid, 14);
	DDX_Text(pDX, IDC_COMPID, m_compid);
	DDV_MaxChars(pDX, m_compid, 12);
	DDX_Text(pDX, IDC_PASSWORD, m_password);
	DDV_MaxChars(pDX, m_password, 8);
	DDX_Check(pDX, IDC_SENDEMAIL, m_sendemail);
	DDX_Check(pDX, IDC_TEAM, m_team);
	//}}AFX_DATA_MAP
	{
		int	can_change_userid;
		can_change_userid =
			USERID[0] == 0 ||
			IniGetNumLines (WORKTODO_FILE) == 0;
//		c_userid.SetReadOnly (!can_change_userid);
//		c_password.SetReadOnly (!can_change_userid);
	}
	c_name_text.SetWindowText (m_team ? "Team name:" : "Your name:");
	c_id_text.SetWindowText (m_team ? "Team ID:" : "Your user ID:");
	c_pwd_text.SetWindowText (m_team ? "Team password:" : "Your password:");
}


BEGIN_MESSAGE_MAP(CUserDlg, CDialog)
	//{{AFX_MSG_MAP(CUserDlg)
	ON_BN_CLICKED(IDC_TEAM, OnTeam)
	//}}AFX_MSG_MAP
END_MESSAGE_MAP()

/////////////////////////////////////////////////////////////////////////////
// CUserDlg message handlers


void CUserDlg::OnTeam() 
{
	UpdateData ();
}
