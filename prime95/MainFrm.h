// MainFrm.h : interface of the CMainFrame class
//
/////////////////////////////////////////////////////////////////////////////
const int MYWM_TRAYMESSAGE = WM_APP + 100;

class CMainFrame : public CFrameWnd
{
protected: // create from serialization only
	CMainFrame();
	DECLARE_DYNCREATE(CMainFrame)

// Attributes
public:

// Operations
public:

// Overrides
	// ClassWizard generated virtual function overrides
	//{{AFX_VIRTUAL(CMainFrame)
	public:
	virtual BOOL PreCreateWindow(CREATESTRUCT& cs);
	virtual BOOL DestroyWindow();
	protected:
	virtual LRESULT WindowProc(UINT message, WPARAM wParam, LPARAM lParam);
	//}}AFX_VIRTUAL
	afx_msg LONG OnPower( UINT, LONG );
	afx_msg LONG OnTrayMessage( UINT, LONG );
	LRESULT OnTaskBarCreated (WPARAM, LPARAM);

// Implementation
public:
	virtual ~CMainFrame();
#ifdef _DEBUG
	virtual void AssertValid() const;
	virtual void Dump(CDumpContext& dc) const;
#endif

// Generated message map functions
protected:
	//{{AFX_MSG(CMainFrame)
	afx_msg void OnSize(UINT nType, int cx, int cy);
	afx_msg void OnEndSession(BOOL bEnding);
	afx_msg void OnActivateApp(BOOL bActive, DWORD hTask);
	afx_msg void OnSysCommand(UINT nID, LPARAM lParam);
	afx_msg void OnTrayOpenWindow();
	afx_msg void OnStopContinue();
	afx_msg LRESULT OnServiceStop(WPARAM wParam, LPARAM lParam);
	//}}AFX_MSG
	DECLARE_MESSAGE_MAP()
public:
	virtual BOOL Create(LPCTSTR lpszClassName, LPCTSTR lpszWindowName, DWORD dwStyle = WS_OVERLAPPEDWINDOW, const RECT& rect = rectDefault, CWnd* pParentWnd = NULL, LPCTSTR lpszMenuName = NULL, DWORD dwExStyle = 0, CCreateContext* pContext = NULL);
};

/////////////////////////////////////////////////////////////////////////////
