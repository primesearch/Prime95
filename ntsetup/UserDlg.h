// UserDlg.h : header file
//

/////////////////////////////////////////////////////////////////////////////
// CUserDlg dialog

class CUserDlg : public CDialog
{
// Construction
public:
	CUserDlg(CWnd* pParent = NULL);   // standard constructor

// Dialog Data
	//{{AFX_DATA(CUserDlg)
	enum { IDD = IDD_USER };
	CStatic	c_pwd_text;
	CStatic	c_name_text;
	CStatic	c_id_text;
	CEdit	c_userid;
	CEdit	c_password;
	CString	m_email;
	CString	m_name;
	CString	m_userid;
	CString	m_compid;
	CString	m_password;
	BOOL	m_sendemail;
	BOOL	m_team;
	//}}AFX_DATA


// Overrides
	// ClassWizard generated virtual function overrides
	//{{AFX_VIRTUAL(CUserDlg)
	protected:
	virtual void DoDataExchange(CDataExchange* pDX);    // DDX/DDV support
	//}}AFX_VIRTUAL

// Implementation
protected:

	// Generated message map functions
	//{{AFX_MSG(CUserDlg)
	afx_msg void OnTeam();
	//}}AFX_MSG
	DECLARE_MESSAGE_MAP()
};
