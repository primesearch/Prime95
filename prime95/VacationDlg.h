// VacationDlg.h : header file
//

/////////////////////////////////////////////////////////////////////////////
// CVacationDlg dialog

class CVacationDlg : public CDialog
{
// Construction
public:
	CVacationDlg(CWnd* pParent = NULL);   // standard constructor

// Dialog Data
	//{{AFX_DATA(CVacationDlg)
	enum { IDD = IDD_VACATION };
	BOOL	m_computer_on;
	UINT	m_vacation_days;
	//}}AFX_DATA


// Overrides
	// ClassWizard generated virtual function overrides
	//{{AFX_VIRTUAL(CVacationDlg)
	protected:
	virtual void DoDataExchange(CDataExchange* pDX);    // DDX/DDV support
	//}}AFX_VIRTUAL

// Implementation
protected:

	// Generated message map functions
	//{{AFX_MSG(CVacationDlg)
		// NOTE: the ClassWizard will add member functions here
	//}}AFX_MSG
	DECLARE_MESSAGE_MAP()
};
