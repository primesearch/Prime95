// Priority.h : header file
//

/////////////////////////////////////////////////////////////////////////////
// CPriority dialog

class CPriority : public CDialog
{
// Construction
public:
	CPriority(CWnd* pParent = NULL);   // standard constructor

// Dialog Data
	//{{AFX_DATA(CPriority)
	enum { IDD = IDD_PRIORITY };
	UINT	m_priority;
	//}}AFX_DATA


// Overrides
	// ClassWizard generated virtual function overrides
	//{{AFX_VIRTUAL(CPriority)
	protected:
	virtual void DoDataExchange(CDataExchange* pDX);    // DDX/DDV support
	//}}AFX_VIRTUAL

// Implementation
protected:

	// Generated message map functions
	//{{AFX_MSG(CPriority)
		// NOTE: the ClassWizard will add member functions here
	//}}AFX_MSG
	DECLARE_MESSAGE_MAP()
};
