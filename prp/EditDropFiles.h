// EditDropFiles.h : header file
//

/////////////////////////////////////////////////////////////////////////////
// CEditDropFiles window

class CEditDropFiles : public CEdit
{
// Construction
public:
	CEditDropFiles();

// Attributes
public:

// Operations
public:

// Overrides
	// ClassWizard generated virtual function overrides
	//{{AFX_VIRTUAL(CEditDropFiles)
	//}}AFX_VIRTUAL

// Implementation
public:
	virtual ~CEditDropFiles();

	// Generated message map functions
protected:
	//{{AFX_MSG(CEditDropFiles)
	afx_msg void OnDropFiles(HDROP hDropInfo);
	//}}AFX_MSG

	DECLARE_MESSAGE_MAP()
};

/////////////////////////////////////////////////////////////////////////////
