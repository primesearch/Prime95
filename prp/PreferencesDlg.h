// PreferencesDlg.h : header file
//

/////////////////////////////////////////////////////////////////////////////
// CPreferencesDlg dialog

class CPreferencesDlg : public CDialog
{
// Construction
public:
	CPreferencesDlg(CWnd* pParent = NULL);   // standard constructor

// Dialog Data
	//{{AFX_DATA(CPreferencesDlg)
	enum { IDD = IDD_PREFERENCES };
	UINT	m_iter;
	UINT	m_disk_write_time;
	BOOL	m_backup;
	UINT	m_r_iter;
	BOOL	m_battery;
	BOOL	m_cumulative;
	//}}AFX_DATA


// Overrides
	// ClassWizard generated virtual function overrides
	//{{AFX_VIRTUAL(CPreferencesDlg)
	protected:
	virtual void DoDataExchange(CDataExchange* pDX);    // DDX/DDV support
	//}}AFX_VIRTUAL

// Implementation
protected:

	// Generated message map functions
	//{{AFX_MSG(CPreferencesDlg)
		// NOTE: the ClassWizard will add member functions here
	//}}AFX_MSG
	DECLARE_MESSAGE_MAP()
};
