// EcmDlg.h : header file
//

/////////////////////////////////////////////////////////////////////////////
// CEcmDlg dialog

class CEcmDlg : public CDialog
{
// Construction
public:
	CEcmDlg(CWnd* pParent = NULL);   // standard constructor

// Dialog Data
	//{{AFX_DATA(CEcmDlg)
	enum { IDD = IDD_ECM };
	UINT	m_p;
	UINT	m_bound1;
	UINT	m_bound2;
	double	m_curve;
	UINT	m_num_curves;
	BOOL	m_plus1;
	//}}AFX_DATA


// Overrides
	// ClassWizard generated virtual function overrides
	//{{AFX_VIRTUAL(CEcmDlg)
	protected:
	virtual void DoDataExchange(CDataExchange* pDX);    // DDX/DDV support
	//}}AFX_VIRTUAL

// Implementation
protected:

	// Generated message map functions
	//{{AFX_MSG(CEcmDlg)
		// NOTE: the ClassWizard will add member functions here
	//}}AFX_MSG
	DECLARE_MESSAGE_MAP()
};
