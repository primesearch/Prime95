// AffinityDlg.h : header file
//

/////////////////////////////////////////////////////////////////////////////
// CAffinityDlg dialog

class CAffinityDlg : public CDialog
{
// Construction
public:
	CAffinityDlg(CWnd* pParent = NULL);   // standard constructor

// Dialog Data
	//{{AFX_DATA(CAffinityDlg)
	enum { IDD = IDD_AFFINITY };
	CEdit	c_affinity;
	CStatic	c_affinity_text;
	BOOL	m_all_cpus;
	UINT	m_cpu;
	//}}AFX_DATA


// Overrides
	// ClassWizard generated virtual function overrides
	//{{AFX_VIRTUAL(CAffinityDlg)
	protected:
	virtual void DoDataExchange(CDataExchange* pDX);    // DDX/DDV support
	//}}AFX_VIRTUAL

// Implementation
protected:

	// Generated message map functions
	//{{AFX_MSG(CAffinityDlg)
	afx_msg void OnAllCpus();
	//}}AFX_MSG
	DECLARE_MESSAGE_MAP()
};
