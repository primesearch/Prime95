// Password.h : header file
//

/////////////////////////////////////////////////////////////////////////////
// Password dialog

class Password : public CDialog
{
// Construction
public:
	Password(CWnd* pParent = NULL);   // standard constructor

// Dialog Data
	//{{AFX_DATA(Password)
	enum { IDD = IDD_PASSWORD };
	UINT	m_password;
	//}}AFX_DATA


// Overrides
	// ClassWizard generated virtual function overrides
	//{{AFX_VIRTUAL(Password)
	protected:
	virtual void DoDataExchange(CDataExchange* pDX);    // DDX/DDV support
	//}}AFX_VIRTUAL

// Implementation
protected:

	// Generated message map functions
	//{{AFX_MSG(Password)
		// NOTE: the ClassWizard will add member functions here
	//}}AFX_MSG
	DECLARE_MESSAGE_MAP()
};
