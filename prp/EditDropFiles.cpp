// EditDropFiles.cpp : implementation file
//

#include "stdafx.h"
#include "Prime95.h"

#ifdef _DEBUG
#define new DEBUG_NEW
#undef THIS_FILE
static char THIS_FILE[] = __FILE__;
#endif

/////////////////////////////////////////////////////////////////////////////
// CEditDropFiles

CEditDropFiles::CEditDropFiles()
{
}

CEditDropFiles::~CEditDropFiles()
{
}


BEGIN_MESSAGE_MAP(CEditDropFiles, CEdit)
	//{{AFX_MSG_MAP(CEditDropFiles)
	ON_WM_DROPFILES()
	//}}AFX_MSG_MAP
END_MESSAGE_MAP()

/////////////////////////////////////////////////////////////////////////////
// CEditDropFiles message handlers

void CEditDropFiles::OnDropFiles(HDROP hDropInfo) 
{
	char buf[200];
	
	DragQueryFile (hDropInfo, 0, buf, sizeof (buf));
	SetWindowText (buf);
}
