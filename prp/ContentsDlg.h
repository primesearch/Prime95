// ContentsDlg.h : header file
//

/////////////////////////////////////////////////////////////////////////////
// ContentsDlg dialog

class ContentsDlg : public CDialog
{
// Construction
public:
	ContentsDlg(CWnd* pParent = NULL);   // standard constructor

// Dialog Data
	//{{AFX_DATA(ContentsDlg)
	enum { IDD = IDD_CONTENTS };
	CString	m_helptext;
	//}}AFX_DATA


// Overrides
	// ClassWizard generated virtual function overrides
	//{{AFX_VIRTUAL(ContentsDlg)
	protected:
	virtual void DoDataExchange(CDataExchange* pDX);    // DDX/DDV support
	//}}AFX_VIRTUAL

// Implementation
protected:

	// Generated message map functions
	//{{AFX_MSG(ContentsDlg)
		// NOTE: the ClassWizard will add member functions here
	//}}AFX_MSG
	DECLARE_MESSAGE_MAP()
};
