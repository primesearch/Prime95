// CpuDlg.h : header file
//

/////////////////////////////////////////////////////////////////////////////
// CCpuDlg dialog

class CCpuDlg : public CDialog
{
// Construction
public:
	CCpuDlg(CWnd* pParent = NULL);   // standard constructor

// Dialog Data
	//{{AFX_DATA(CCpuDlg)
	enum { IDD = IDD_CPU };
	CButton	c_pentium_pro;
	CButton	c_pentium;
	CButton	c_k6;
	CButton	c_celeron;
	CButton	c_pii;
	CButton	c_piii;
	CButton	c_k7;
	UINT	m_speed;
	UINT	m_hours;
	int		m_cpu_type;
	CString	m_start_time;
	CString	m_end_time;
	UINT	m_day_memory;
	UINT	m_night_memory;
	//}}AFX_DATA


// Overrides
	// ClassWizard generated virtual function overrides
	//{{AFX_VIRTUAL(CCpuDlg)
	protected:
	virtual void DoDataExchange(CDataExchange* pDX);    // DDX/DDV support
	//}}AFX_VIRTUAL

// Implementation
protected:

	// Generated message map functions
	//{{AFX_MSG(CCpuDlg)
		// NOTE: the ClassWizard will add member functions here
	//}}AFX_MSG
	DECLARE_MESSAGE_MAP()
};
