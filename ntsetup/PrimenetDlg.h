// PrimenetDlg.h : header file
//

/////////////////////////////////////////////////////////////////////////////
// PrimenetDlg dialog

class PrimenetDlg : public CDialog
{
// Construction
public:
	PrimenetDlg(CWnd* pParent = NULL);   // standard constructor

// Dialog Data
	//{{AFX_DATA(PrimenetDlg)
	enum { IDD = IDD_PRIMENET };
	CButton	c_bigones;
	CStatic	c_protocol_text;
	CButton	c_rpc;
	CButton	c_http;
	CButton	c_dialup;
	CButton	c_work_dflt;
	CStatic	c_work_text;
	CEdit	c_work;
	CButton	c_workgroup;
	CButton	c_ll;
	CButton	c_factor;
	CButton	c_dblchk;
	BOOL	m_dblchk;
	BOOL	m_factor;
	BOOL	m_lucas;
	UINT	m_work;
	BOOL	m_work_dflt;
	BOOL	m_primenet;
	BOOL	m_dialup;
	int	m_rpc;
	BOOL	m_bigones;
	//}}AFX_DATA


// Overrides
	// ClassWizard generated virtual function overrides
	//{{AFX_VIRTUAL(PrimenetDlg)
	protected:
	virtual void DoDataExchange(CDataExchange* pDX);    // DDX/DDV support
	//}}AFX_VIRTUAL

// Implementation
protected:

	// Generated message map functions
	//{{AFX_MSG(PrimenetDlg)
	afx_msg void OnWorkDflt();
	afx_msg void OnPrimenet();
	//}}AFX_MSG
	DECLARE_MESSAGE_MAP()
};
