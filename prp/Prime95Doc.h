// Prime95Doc.h : interface of the CPrime95Doc class
//
/////////////////////////////////////////////////////////////////////////////

class CPrime95Doc : public CDocument
{
protected: // create from serialization only
	CPrime95Doc();
	DECLARE_DYNCREATE(CPrime95Doc)

// Attributes
public:
	char	*ReplLine;
	int	Replacing;

// Operations
public:
	void title (char *);
	void rangeStatus ();
	void LineFeed ();
	void OutputStr (char *);
	void ReplaceableLine (int);
	void Restart1 ();
	void Restart2 ();

// Overrides
	// ClassWizard generated virtual function overrides
	//{{AFX_VIRTUAL(CPrime95Doc)
	public:
	virtual BOOL OnNewDocument();
	virtual void Serialize(CArchive& ar);
	virtual void OnCloseDocument();
	//}}AFX_VIRTUAL

// Implementation
public:
	virtual ~CPrime95Doc();
#ifdef _DEBUG
	virtual void AssertValid() const;
	virtual void Dump(CDumpContext& dc) const;
#endif

protected:

// Generated message map functions
protected:
	//{{AFX_MSG(CPrime95Doc)
	afx_msg void OnContinue();
	afx_msg void OnUpdateContinue(CCmdUI* pCmdUI);
	afx_msg void OnStop();
	afx_msg void OnUpdateStop(CCmdUI* pCmdUI);
	afx_msg void OnErrchk();
	afx_msg void OnUpdateErrchk(CCmdUI* pCmdUI);
	afx_msg void OnCpu();
	afx_msg void OnPreferences();
	afx_msg void OnTest();
	afx_msg void OnTray();
	afx_msg void OnUpdateTray(CCmdUI* pCmdUI);
	afx_msg void OnHide();
	afx_msg void OnUpdateHide(CCmdUI* pCmdUI);
	afx_msg void OnPriority();
	afx_msg void OnService();
	afx_msg void OnUpdateService(CCmdUI* pCmdUI);
	afx_msg void OnAffinity();
	afx_msg void OnUpdateAffinity(CCmdUI* pCmdUI);
	afx_msg void OnUpdateTest(CCmdUI* pCmdUI);
	afx_msg void OnHelpContents();
	//}}AFX_MSG
	DECLARE_MESSAGE_MAP()
};

/////////////////////////////////////////////////////////////////////////////

extern CPrime95Doc *OUTPUT_STR_HACK;
