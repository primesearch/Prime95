// TestDlg.cpp : implementation file
//

#include "stdafx.h"
#include "Prime95.h"
#include "TestDlg.h"

#ifdef _DEBUG
#define new DEBUG_NEW
#undef THIS_FILE
static char THIS_FILE[] = __FILE__;
#endif

/////////////////////////////////////////////////////////////////////////////
// CTestDlg dialog


CTestDlg::CTestDlg(CWnd* pParent /*=NULL*/)
	: CDialog(CTestDlg::IDD, pParent)
{
	//{{AFX_DATA_INIT(CTestDlg)
	m_work = -1;
	m_expr = _T("");
	m_exprfile = _T("");
	m_exprfile_line = 0;
	m_pgen_input = _T("");
	m_pgen_line = 0;
	m_pgen_output = _T("");
	//}}AFX_DATA_INIT
}

void CTestDlg::DoDataExchange(CDataExchange* pDX)
{
	CDialog::DoDataExchange(pDX);
	//{{AFX_DATA_MAP(CTestDlg)
	DDX_Control(pDX, IDC_PGEN_OUTPUT, e_pgen_output);
	DDX_Control(pDX, IDC_PGEN_LINE, e_pgen_line);
	DDX_Control(pDX, IDC_PGEN_INPUT, e_pgen_input);
	DDX_Control(pDX, IDC_EXPRFILE_LINE, e_exprfile_line);
	DDX_Control(pDX, IDC_EXPRFILE, e_exprfile);
	DDX_Control(pDX, IDC_EXPR, e_expr);
	DDX_Control(pDX, IDS_PGEN_OUTPUT, c_pgen_output);
	DDX_Control(pDX, IDS_PGEN_LINE, c_pgen_line);
	DDX_Control(pDX, IDS_PGEN_INPUT, c_pgen_input);
	DDX_Control(pDX, IDS_EXPRFILE_LINE, c_exprfile_line);
	DDX_Control(pDX, IDS_EXPRFILE, c_exprfile);
	DDX_Radio(pDX, IDC_WORK_PGEN, m_work);
	DDX_Text(pDX, IDC_EXPR, m_expr);
	DDX_Text(pDX, IDC_EXPRFILE, m_exprfile);
	DDX_Text(pDX, IDC_EXPRFILE_LINE, m_exprfile_line);
	DDX_Text(pDX, IDC_PGEN_INPUT, m_pgen_input);
	DDX_Text(pDX, IDC_PGEN_LINE, m_pgen_line);
	DDX_Text(pDX, IDC_PGEN_OUTPUT, m_pgen_output);
	//}}AFX_DATA_MAP
	c_pgen_input.EnableWindow (m_work == 0);
	c_pgen_output.EnableWindow (m_work == 0);
	c_pgen_line.EnableWindow (m_work == 0);
	c_exprfile.EnableWindow (m_work == 2);
	c_exprfile_line.EnableWindow (m_work == 2);
	e_pgen_input.EnableWindow (m_work == 0);
	e_pgen_output.EnableWindow (m_work == 0);
	e_pgen_line.EnableWindow (m_work == 0);
	e_expr.EnableWindow (m_work == 1);
	e_exprfile.EnableWindow (m_work == 2);
	e_exprfile_line.EnableWindow (m_work == 2);

c_exprfile.EnableWindow (0);
c_exprfile_line.EnableWindow (0);
e_expr.EnableWindow (0);
e_exprfile.EnableWindow (0);
e_exprfile_line.EnableWindow (0);
}


BEGIN_MESSAGE_MAP(CTestDlg, CDialog)
	//{{AFX_MSG_MAP(CTestDlg)
	ON_BN_CLICKED(IDC_WORK_EXPR, OnWork)
	ON_BN_CLICKED(IDC_WORK_EXPRFILE, OnWork)
	ON_BN_CLICKED(IDC_WORK_PGEN, OnWork)
	//}}AFX_MSG_MAP
END_MESSAGE_MAP()

/////////////////////////////////////////////////////////////////////////////
// CTestDlg message handlers

void CTestDlg::OnWork() 
{
	UpdateData ();
}
