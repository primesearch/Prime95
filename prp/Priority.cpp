// Priority.cpp : implementation file
//

#include "stdafx.h"
#include "Prime95.h"
#include "Priority.h"

#ifdef _DEBUG
#define new DEBUG_NEW
#undef THIS_FILE
static char THIS_FILE[] = __FILE__;
#endif

/////////////////////////////////////////////////////////////////////////////
// CPriority dialog


CPriority::CPriority(CWnd* pParent /*=NULL*/)
	: CDialog(CPriority::IDD, pParent)
{
	//{{AFX_DATA_INIT(CPriority)
	m_priority = 0;
	//}}AFX_DATA_INIT
}


void CPriority::DoDataExchange(CDataExchange* pDX)
{
	CDialog::DoDataExchange(pDX);
	//{{AFX_DATA_MAP(CPriority)
	DDX_Text(pDX, IDC_PRIORITY, m_priority);
	DDV_MinMaxUInt(pDX, m_priority, 1, 10);
	//}}AFX_DATA_MAP
}


BEGIN_MESSAGE_MAP(CPriority, CDialog)
	//{{AFX_MSG_MAP(CPriority)
		// NOTE: the ClassWizard will add message map macros here
	//}}AFX_MSG_MAP
END_MESSAGE_MAP()

/////////////////////////////////////////////////////////////////////////////
// CPriority message handlers
