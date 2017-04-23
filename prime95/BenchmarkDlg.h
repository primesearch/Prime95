// BenchmarkDlg.h : header file
//
// Copyright 2017 Mersenne Research, Inc.  All rights reserved
//

/////////////////////////////////////////////////////////////////////////////
// CBenchmarkDlg dialog

class CBenchmarkDlg : public CDialog
{
// Construction
public:
	CBenchmarkDlg(CWnd* pParent = NULL);   // standard constructor

// Dialog Data
	//{{AFX_DATA(CBenchmarkDlg)
	enum { IDD = IDD_BENCHMARK };
	CComboBox c_bench_type;
	int	m_bench_type;

	CStatic	c_minFFT_text;
	CEdit	c_minFFT;
	UINT	m_minFFT;
	CStatic	c_maxFFT_text;
	CEdit	c_maxFFT;
	UINT	m_maxFFT;
	CButton	c_all_FFT_sizes;
	BOOL	m_all_FFT_sizes;

	CButton	c_bench_one_core;
	BOOL	m_bench_one_core;
	CButton	c_bench_all_cores;
	BOOL	m_bench_all_cores;
	CButton	c_bench_in_between_cores;
	BOOL	m_bench_in_between_cores;
	CButton	c_hyperthreading;
	BOOL	m_hyperthreading;

	CButton	c_bench_one_worker;
	BOOL	m_bench_one_worker;
	CButton	c_bench_max_workers;
	BOOL	m_bench_max_workers;
	CButton	c_bench_in_between_workers;
	BOOL	m_bench_in_between_workers;
	CStatic	c_bench_time_text;
	CEdit	c_bench_time;
	UINT	m_bench_time;
	//}}AFX_DATA

// Overrides
	// ClassWizard generated virtual function overrides
	//{{AFX_VIRTUAL(CBenchmarkDlg)
	protected:
	virtual void DoDataExchange(CDataExchange* pDX);    // DDX/DDV support
	//}}AFX_VIRTUAL

// Implementation
protected:

	// Generated message map functions
	//{{AFX_MSG(CBenchmarkDlg)
	virtual BOOL OnInitDialog();
	afx_msg void OnEnChangeBenchType();
	afx_msg void OnCbnKillfocusBenchType();
	afx_msg void OnEnKillfocusMinFFT();
	afx_msg void OnEnKillfocusMaxFFT();
	//}}AFX_MSG
	DECLARE_MESSAGE_MAP()
};
