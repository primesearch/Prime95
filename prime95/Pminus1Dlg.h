// Pminus1Dlg.h : header file
//

/////////////////////////////////////////////////////////////////////////////
// CPminus1Dlg dialog

class CPminus1Dlg : public CDialog
{
// Construction
public:
	CPminus1Dlg(CWnd* pParent = NULL);   // standard constructor

// Dialog Data
	//{{AFX_DATA(CPminus1Dlg)
	enum { IDD = IDD_PMINUS1 };
	UINT	m_p;
	UINT	m_bound1;
	UINT	m_bound2;
	BOOL	m_plus1;
	//}}AFX_DATA


// Overrides
	// ClassWizard generated virtual function overrides
	//{{AFX_VIRTUAL(CPminus1Dlg)
	protected:
	virtual void DoDataExchange(CDataExchange* pDX);    // DDX/DDV support
	//}}AFX_VIRTUAL

// Implementation
protected:

	// Generated message map functions
	//{{AFX_MSG(CPminus1Dlg)
		// NOTE: the ClassWizard will add member functions here
	//}}AFX_MSG
	DECLARE_MESSAGE_MAP()
};
