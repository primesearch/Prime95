// AffinityDlg.cpp : implementation file
//

#include "stdafx.h"
#include "Prime95.h"
#include "AffinityDlg.h"

#ifdef _DEBUG
#define new DEBUG_NEW
#undef THIS_FILE
static char THIS_FILE[] = __FILE__;
#endif

/////////////////////////////////////////////////////////////////////////////
// CAffinityDlg dialog


CAffinityDlg::CAffinityDlg(CWnd* pParent /*=NULL*/)
	: CDialog(CAffinityDlg::IDD, pParent)
{
	//{{AFX_DATA_INIT(CAffinityDlg)
	m_all_cpus = FALSE;
	m_cpu = 0;
	//}}AFX_DATA_INIT
}


void CAffinityDlg::DoDataExchange(CDataExchange* pDX)
{
	CDialog::DoDataExchange(pDX);
	//{{AFX_DATA_MAP(CAffinityDlg)
	DDX_Control(pDX, IDC_AFFINITY, c_affinity);
	DDX_Control(pDX, IDC_AFFINITY_TEXT, c_affinity_text);
	DDX_Check(pDX, IDC_ALL_CPUS, m_all_cpus);
	DDX_Text(pDX, IDC_AFFINITY, m_cpu);
	DDV_MinMaxUInt(pDX, m_cpu, 0, 31);
	//}}AFX_DATA_MAP
	c_affinity_text.EnableWindow (!m_all_cpus);
	c_affinity.EnableWindow (!m_all_cpus);
}


BEGIN_MESSAGE_MAP(CAffinityDlg, CDialog)
	//{{AFX_MSG_MAP(CAffinityDlg)
	ON_BN_CLICKED(IDC_ALL_CPUS, OnAllCpus)
	//}}AFX_MSG_MAP
END_MESSAGE_MAP()

/////////////////////////////////////////////////////////////////////////////
// CAffinityDlg message handlers

void CAffinityDlg::OnAllCpus() 
{
	UpdateData ();
}
