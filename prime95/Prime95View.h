// Prime95View.h : interface of the CPrime95View class
//
/////////////////////////////////////////////////////////////////////////////

class CPrime95View : public CView
{
protected: // create from serialization only
	CPrime95View();
	DECLARE_DYNCREATE(CPrime95View)

// Attributes
public:
	CPrime95Doc* GetDocument();

// Operations
public:

// Overrides
	// ClassWizard generated virtual function overrides
	//{{AFX_VIRTUAL(CPrime95View)
	public:
	virtual void OnDraw(CDC* pDC);  // overridden to draw this view
	virtual BOOL PreCreateWindow(CREATESTRUCT& cs);
	protected:
	virtual void OnUpdate(CView* pSender, LPARAM lHint, CObject* pHint);
	//}}AFX_VIRTUAL

// Implementation
public:
	virtual ~CPrime95View();
#ifdef _DEBUG
	virtual void AssertValid() const;
	virtual void Dump(CDumpContext& dc) const;
#endif

protected:

// Generated message map functions
protected:
	//{{AFX_MSG(CPrime95View)
		// NOTE - the ClassWizard will add and remove member functions here.
		//    DO NOT EDIT what you see in these blocks of generated code !
	//}}AFX_MSG
	DECLARE_MESSAGE_MAP()
};

#ifndef _DEBUG  // debug version in Prime95View.cpp
inline CPrime95Doc* CPrime95View::GetDocument()
   { return (CPrime95Doc*)m_pDocument; }
#endif

/////////////////////////////////////////////////////////////////////////////
