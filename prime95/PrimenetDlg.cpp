// PrimenetDlg.cpp : implementation file
//

#include "stdafx.h"
#include "Prime95.h"
#include "PrimenetDlg.h"

#ifdef _DEBUG
#define new DEBUG_NEW
#undef THIS_FILE
static char THIS_FILE[] = __FILE__;
#endif

/////////////////////////////////////////////////////////////////////////////
// PrimenetDlg dialog


PrimenetDlg::PrimenetDlg(CWnd* pParent /*=NULL*/)
	: CDialog(PrimenetDlg::IDD, pParent)
{
	//{{AFX_DATA_INIT(PrimenetDlg)
	m_dblchk = FALSE;
	m_factor = FALSE;
	m_lucas = FALSE;
	m_work = 0;
	m_work_dflt = FALSE;
	m_primenet = FALSE;
	m_dialup = FALSE;
	m_bigones = FALSE;
	m_pfactor = FALSE;
	//}}AFX_DATA_INIT
}


void PrimenetDlg::DoDataExchange(CDataExchange* pDX)
{
	CDialog::DoDataExchange(pDX);
	//{{AFX_DATA_MAP(PrimenetDlg)
	DDX_Control(pDX, IDC_BIGONES, c_bigones);
	DDX_Control(pDX, IDC_DIALUP, c_dialup);
	DDX_Control(pDX, IDC_WORK_DFLT, c_work_dflt);
	DDX_Control(pDX, IDC_WORK_TEXT, c_work_text);
	DDX_Control(pDX, IDC_WORK, c_work);
	DDX_Control(pDX, IDC_WORKGROUP, c_workgroup);
	DDX_Control(pDX, IDC_LL, c_ll);
	DDX_Control(pDX, IDC_FACTOR, c_factor);
	DDX_Control(pDX, IDC_DBLCHK, c_dblchk);
	DDX_Check(pDX, IDC_DBLCHK, m_dblchk);
	DDX_Check(pDX, IDC_FACTOR, m_factor);
	DDX_Check(pDX, IDC_LL, m_lucas);
	DDX_Text(pDX, IDC_WORK, m_work);
	DDV_MinMaxUInt(pDX, m_work, 1, 90);
	DDX_Check(pDX, IDC_WORK_DFLT, m_work_dflt);
	DDX_Check(pDX, IDC_PRIMENET, m_primenet);
	DDX_Check(pDX, IDC_DIALUP, m_dialup);
	DDX_Check(pDX, IDC_BIGONES, m_bigones);
	//}}AFX_DATA_MAP
	c_dialup.EnableWindow (m_primenet);
	c_work_text.EnableWindow (m_primenet);
	c_work.EnableWindow (m_primenet);
	c_work_dflt.EnableWindow (m_primenet);
	c_workgroup.EnableWindow (m_primenet && !m_work_dflt);
	c_bigones.EnableWindow (m_primenet && !m_work_dflt);
	c_dblchk.EnableWindow (m_primenet && !m_work_dflt);
	c_factor.EnableWindow (m_primenet && !m_work_dflt);
	c_ll.EnableWindow (m_primenet && !m_work_dflt);
}


BEGIN_MESSAGE_MAP(PrimenetDlg, CDialog)
	//{{AFX_MSG_MAP(PrimenetDlg)
	ON_BN_CLICKED(IDC_WORK_DFLT, OnWorkDflt)
	ON_BN_CLICKED(IDC_PRIMENET, OnPrimenet)
	//}}AFX_MSG_MAP
END_MESSAGE_MAP()

/////////////////////////////////////////////////////////////////////////////
// PrimenetDlg message handlers

void PrimenetDlg::OnWorkDflt() 
{
	UpdateData ();		// Get the values from the dialog box
}

void PrimenetDlg::OnPrimenet() 
{
	UpdateData ();		// Get the values from the dialog box
}
