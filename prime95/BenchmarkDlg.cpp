// BenchmarkDlg.cpp : implementation file
//
// Copyright 2017 Mersenne Research, Inc.  All rights reserved
//

#include "stdafx.h"
#include "Prime95.h"
#include "BenchmarkDlg.h"

#ifdef _DEBUG
#define new DEBUG_NEW
#undef THIS_FILE
static char THIS_FILE[] = __FILE__;
#endif

/////////////////////////////////////////////////////////////////////////////
// CBenchmarkDlg dialog


CBenchmarkDlg::CBenchmarkDlg(CWnd* pParent /*=NULL*/)
	: CDialog(CBenchmarkDlg::IDD, pParent)
{
	//{{AFX_DATA_INIT(CBenchmarkDlg)
	m_bench_type = 0;
	m_minFFT = 0;
	m_maxFFT = 0;
	m_all_FFT_sizes = 0;
	m_bench_time = 0;
	m_bench_one_core = 0;
	m_bench_all_cores = 0;
	m_bench_in_between_cores = 0;
	m_hyperthreading = 0;
	//}}AFX_DATA_INIT
}


void CBenchmarkDlg::DoDataExchange(CDataExchange* pDX)
{
	CDialog::DoDataExchange(pDX);
	//{{AFX_DATA_MAP(CBenchmarkDlg)
	DDX_Control(pDX, IDC_COMBO1, c_bench_type);

	DDX_Control(pDX, IDC_MINFFT_TEXT, c_minFFT_text);
	DDX_Control(pDX, IDC_MINFFT, c_minFFT);
	DDX_Text(pDX, IDC_MINFFT, m_minFFT);
	DDV_MinMaxUInt(pDX, m_minFFT, 1, 32768);
	DDX_Control(pDX, IDC_MAXFFT_TEXT, c_maxFFT_text);
	DDX_Control(pDX, IDC_MAXFFT, c_maxFFT);
	DDX_Text(pDX, IDC_MAXFFT, m_maxFFT);
	DDV_MinMaxUInt(pDX, m_maxFFT, 1, 32768);
	DDX_Control(pDX, IDC_ALL_FFT_SIZES, c_all_FFT_sizes);
	DDX_Check(pDX, IDC_ALL_FFT_SIZES, m_all_FFT_sizes);

	DDX_Control(pDX, IDC_ONE_CORE, c_bench_one_core);
	DDX_Check(pDX, IDC_ONE_CORE, m_bench_one_core);
	DDX_Control(pDX, IDC_ALL_CORES, c_bench_all_cores);
	DDX_Check(pDX, IDC_ALL_CORES, m_bench_all_cores);
	DDX_Control(pDX, IDC_OTHER_CORES, c_bench_in_between_cores);
	DDX_Check(pDX, IDC_OTHER_CORES, m_bench_in_between_cores);
	DDX_Control(pDX, IDC_HYPERTHREADING, c_hyperthreading);
	DDX_Check(pDX, IDC_HYPERTHREADING, m_hyperthreading);

	DDX_Control(pDX, IDC_ONE_WORKER, c_bench_one_worker);
	DDX_Check(pDX, IDC_ONE_WORKER, m_bench_one_worker);
	DDX_Control(pDX, IDC_MAX_WORKERS, c_bench_max_workers);
	DDX_Check(pDX, IDC_MAX_WORKERS, m_bench_max_workers);
	DDX_Control(pDX, IDC_OTHER_WORKERS, c_bench_in_between_workers);
	DDX_Check(pDX, IDC_OTHER_WORKERS, m_bench_in_between_workers);
	DDX_Control(pDX, IDC_TIMEFFT_TEXT, c_bench_time_text);
	DDX_Control(pDX, IDC_TIMEFFT, c_bench_time);
	DDX_Text(pDX, IDC_TIMEFFT, m_bench_time);
	DDV_MinMaxUInt(pDX, m_bench_time, 5, 60);
	//}}AFX_DATA_MAP
	c_minFFT_text.EnableWindow (m_bench_type != 2);
	c_minFFT.EnableWindow (m_bench_type != 2);
	c_maxFFT_text.EnableWindow (m_bench_type != 2);
	c_maxFFT.EnableWindow (m_bench_type != 2);
	c_all_FFT_sizes.EnableWindow (m_bench_type != 2 && m_minFFT != m_maxFFT);
	c_bench_one_core.EnableWindow (NUM_CPUS > 1);
	c_bench_all_cores.EnableWindow (NUM_CPUS > 1);
	c_bench_in_between_cores.EnableWindow (NUM_CPUS > 2);
	c_hyperthreading.EnableWindow (CPU_HYPERTHREADS > 1);
	c_bench_one_worker.EnableWindow (m_bench_type == 0 && NUM_CPUS > 1);
	c_bench_max_workers.EnableWindow (m_bench_type == 0 && NUM_CPUS > 1);
	c_bench_in_between_workers.EnableWindow (m_bench_type == 0 && NUM_CPUS > 2);
	c_bench_time_text.EnableWindow (m_bench_type == 0);
	c_bench_time.EnableWindow (m_bench_type == 0);
}


BEGIN_MESSAGE_MAP(CBenchmarkDlg, CDialog)
	//{{AFX_MSG_MAP(CBenchmarkDlg)
	ON_EN_CHANGE(IDC_COMBO1, &CBenchmarkDlg::OnEnChangeBenchType)
	ON_CBN_KILLFOCUS(IDC_COMBO1, &CBenchmarkDlg::OnCbnKillfocusBenchType)
	ON_EN_KILLFOCUS(IDC_MINFFT, &CBenchmarkDlg::OnEnKillfocusMinFFT)
	ON_EN_KILLFOCUS(IDC_MAXFFT, &CBenchmarkDlg::OnEnKillfocusMaxFFT)
	//}}AFX_MSG_MAP
END_MESSAGE_MAP()

/////////////////////////////////////////////////////////////////////////////
// CBenchmarkDlg message handlers

BOOL CBenchmarkDlg::OnInitDialog()
{
	CDialog::OnInitDialog();

// Populate the benchmark type combo box

	c_bench_type.ResetContent ();
	c_bench_type.AddString ("Throughput benchmark");
	c_bench_type.AddString ("FFT timings benchmark");
	c_bench_type.AddString ("Trial factoring benchmark");
	c_bench_type.SetCurSel (0);

	return TRUE;  // return TRUE  unless you set the focus to a control
}

void CBenchmarkDlg::OnEnChangeBenchType() 
{
	m_bench_type = c_bench_type.GetCurSel ();
	UpdateData ();		// Get the values from the dialog box
}

void CBenchmarkDlg::OnCbnKillfocusBenchType()
{
	m_bench_type = c_bench_type.GetCurSel ();
	UpdateData ();		// Get the values from the dialog box
}

void CBenchmarkDlg::OnEnKillfocusMinFFT()
{
	UpdateData ();		// Get the values from the dialog box
}

void CBenchmarkDlg::OnEnKillfocusMaxFFT()
{
	UpdateData ();		// Get the values from the dialog box
}

