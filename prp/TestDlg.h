// TestDlg.h : header file
//

/////////////////////////////////////////////////////////////////////////////
// CTestDlg dialog

class CTestDlg : public CDialog
{
// Construction
public:
	CTestDlg(CWnd* pParent = NULL);   // standard constructor

// Dialog Data
	//{{AFX_DATA(CTestDlg)
	enum { IDD = IDD_TEST };
	CEditDropFiles	e_pgen_output;
	CEdit	e_pgen_line;
	CEditDropFiles	e_pgen_input;
	CEdit	e_exprfile_line;
	CEditDropFiles	e_exprfile;
	CEdit	e_expr;
	CStatic	c_pgen_output;
	CStatic	c_pgen_line;
	CStatic	c_pgen_input;
	CStatic	c_exprfile_line;
	CStatic	c_exprfile;
	int		m_work;
	CString	m_expr;
	CString	m_exprfile;
	int		m_exprfile_line;
	CString	m_pgen_input;
	int		m_pgen_line;
	CString	m_pgen_output;
	//}}AFX_DATA


// Overrides
	// ClassWizard generated virtual function overrides
	//{{AFX_VIRTUAL(CTestDlg)
	protected:
	virtual void DoDataExchange(CDataExchange* pDX);    // DDX/DDV support
	//}}AFX_VIRTUAL

// Implementation
protected:

	// Generated message map functions
	//{{AFX_MSG(CTestDlg)
	afx_msg void OnWork();
	//}}AFX_MSG
	DECLARE_MESSAGE_MAP()
};
