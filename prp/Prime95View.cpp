// Prime95View.cpp : implementation of the CPrime95View class
//

#include "stdafx.h"
#include "Prime95.h"

#include "Prime95Doc.h"
#include "Prime95View.h"

#ifdef _DEBUG
#define new DEBUG_NEW
#undef THIS_FILE
static char THIS_FILE[] = __FILE__;
#endif

/////////////////////////////////////////////////////////////////////////////
// CPrime95View

IMPLEMENT_DYNCREATE(CPrime95View, CView)

BEGIN_MESSAGE_MAP(CPrime95View, CView)
	//{{AFX_MSG_MAP(CPrime95View)
		// NOTE - the ClassWizard will add and remove mapping macros here.
		//    DO NOT EDIT what you see in these blocks of generated code!
	//}}AFX_MSG_MAP
END_MESSAGE_MAP()

/////////////////////////////////////////////////////////////////////////////
// CPrime95View construction/destruction

CPrime95View::CPrime95View()
{
	// TODO: add construction code here

}

CPrime95View::~CPrime95View()
{
}

BOOL CPrime95View::PreCreateWindow(CREATESTRUCT& cs)
{
	// TODO: Modify the Window class or styles here by modifying
	//  the CREATESTRUCT cs

	return CView::PreCreateWindow(cs);
}

/////////////////////////////////////////////////////////////////////////////
// CPrime95View drawing

void CPrime95View::OnDraw(CDC* pDC)
{
	if (lines[0] == NULL) {
		CPrime95Doc* pDoc = GetDocument();
		pDoc->title ("IDLE");
		return;
	}

	RECT	r;
	int	ypos;
	int	i;

	pDC->SetBkMode (TRANSPARENT);
	pDC->SetTextColor (GetSysColor (COLOR_WINDOWTEXT));

	GetClientRect (&r);
	for (ypos = r.bottom, i = 1;
	     ypos > 0 && i < NumLines;
	     i++) {
		ypos -= charHeight;
		pDC->TextOut (0, ypos, lines[i], strlen (lines[i]));
	}
}

/////////////////////////////////////////////////////////////////////////////
// CPrime95View diagnostics

#ifdef _DEBUG
void CPrime95View::AssertValid() const
{
//	CView::AssertValid();
}

void CPrime95View::Dump(CDumpContext& dc) const
{
	CView::Dump(dc);
}

CPrime95Doc* CPrime95View::GetDocument() // non-debug version is inline
{
	ASSERT(m_pDocument->IsKindOf(RUNTIME_CLASS(CPrime95Doc)));
	return (CPrime95Doc*)m_pDocument;
}
#endif //_DEBUG

/////////////////////////////////////////////////////////////////////////////
// CPrime95View message handlers

void CPrime95View::OnUpdate(CView* pSender, LPARAM lHint, CObject* pHint) 
{
	if (charHeight == 0) {
		CDC *dc = GetDC ();
		CSize size;
		size = dc->GetTextExtent ("A", 1);
		charHeight = size.cy;
	}
	if (lHint == 0) ScrollWindow (0, -charHeight, NULL, NULL);
	else Invalidate ();
	UpdateWindow ();
}
