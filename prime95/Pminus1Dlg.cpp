// Pminus1Dlg.cpp : implementation file
//

#include "stdafx.h"
#include "Prime95.h"
#include "Pminus1Dlg.h"

#ifdef _DEBUG
#define new DEBUG_NEW
#undef THIS_FILE
static char THIS_FILE[] = __FILE__;
#endif

/////////////////////////////////////////////////////////////////////////////
// CPminus1Dlg dialog


CPminus1Dlg::CPminus1Dlg(CWnd* pParent /*=NULL*/)
	: CDialog(CPminus1Dlg::IDD, pParent)
{
	//{{AFX_DATA_INIT(CPminus1Dlg)
	m_p = 0;
	m_bound1 = 0;
	m_bound2 = 0;
	m_plus1 = FALSE;
	//}}AFX_DATA_INIT
}


void CPminus1Dlg::DoDataExchange(CDataExchange* pDX)
{
	CDialog::DoDataExchange(pDX);
	//{{AFX_DATA_MAP(CPminus1Dlg)
	DDX_Text(pDX, IDC_P1, m_p);
	DDV_MinMaxUInt(pDX, m_p, 100, 20500000);
	DDX_Text(pDX, IDC_P2, m_bound1);
	DDV_MinMaxUInt(pDX, m_bound1, 0, 4000000000);
	DDX_Text(pDX, IDC_P3, m_bound2);
	DDV_MinMaxUInt(pDX, m_bound2, 0, 4000000000);
	DDX_Check(pDX, IDC_PLUS1, m_plus1);
	//}}AFX_DATA_MAP
}


BEGIN_MESSAGE_MAP(CPminus1Dlg, CDialog)
	//{{AFX_MSG_MAP(CPminus1Dlg)
		// NOTE: the ClassWizard will add message map macros here
	//}}AFX_MSG_MAP
END_MESSAGE_MAP()

/////////////////////////////////////////////////////////////////////////////
// CPminus1Dlg message handlers
