# Microsoft Developer Studio Generated NMAKE File, Format Version 4.00
# ** DO NOT EDIT **

# TARGTYPE "Win32 (x86) Application" 0x0101

!IF "$(CFG)" == ""
CFG=Prime95 - Win32 Debug
!MESSAGE No configuration specified.  Defaulting to Prime95 - Win32 Debug.
!ENDIF 

!IF "$(CFG)" != "Prime95 - Win32 Release" && "$(CFG)" !=\
 "Prime95 - Win32 Debug"
!MESSAGE Invalid configuration "$(CFG)" specified.
!MESSAGE You can specify a configuration when running NMAKE on this makefile
!MESSAGE by defining the macro CFG on the command line.  For example:
!MESSAGE 
!MESSAGE NMAKE /f "Prime95.mak" CFG="Prime95 - Win32 Debug"
!MESSAGE 
!MESSAGE Possible choices for configuration are:
!MESSAGE 
!MESSAGE "Prime95 - Win32 Release" (based on "Win32 (x86) Application")
!MESSAGE "Prime95 - Win32 Debug" (based on "Win32 (x86) Application")
!MESSAGE 
!ERROR An invalid configuration is specified.
!ENDIF 

!IF "$(OS)" == "Windows_NT"
NULL=
!ELSE 
NULL=nul
!ENDIF 
################################################################################
# Begin Project
# PROP Target_Last_Scanned "Prime95 - Win32 Debug"
CPP=cl.exe
RSC=rc.exe
MTL=mktyplib.exe

!IF  "$(CFG)" == "Prime95 - Win32 Release"

# PROP BASE Use_MFC 6
# PROP BASE Use_Debug_Libraries 0
# PROP BASE Output_Dir "Release"
# PROP BASE Intermediate_Dir "Release"
# PROP BASE Target_Dir ""
# PROP Use_MFC 5
# PROP Use_Debug_Libraries 0
# PROP Output_Dir "Release"
# PROP Intermediate_Dir "Release"
# PROP Target_Dir ""
OUTDIR=.\Release
INTDIR=.\Release

ALL : "$(OUTDIR)\Prime95.exe" "$(OUTDIR)\Prime95.hlp"

CLEAN : 
	-@erase ".\Release\Prime95.exe"
	-@erase ".\Release\AffinityDlg.obj"
	-@erase ".\Release\Prime95.pch"
	-@erase ".\Release\PreferencesDlg.obj"
	-@erase ".\Release\Prime95.obj"
	-@erase ".\Release\Prime95Doc.obj"
	-@erase ".\Release\StdAfx.obj"
	-@erase ".\Release\Priority.obj"
	-@erase ".\Release\EditDropFiles.obj"
	-@erase ".\Release\CpuDlg.obj"
	-@erase ".\Release\Prime95View.obj"
	-@erase ".\Release\TestDlg.obj"
	-@erase ".\Release\MainFrm.obj"
	-@erase ".\Release\ContentsDlg.obj"
	-@erase ".\Release\Prime95.res"
	-@erase ".\Release\Prime95.map"
	-@erase ".\Release\Prime95.hlp"

"$(OUTDIR)" :
    if not exist "$(OUTDIR)/$(NULL)" mkdir "$(OUTDIR)"

# ADD BASE CPP /nologo /MD /W3 /GX /O2 /D "WIN32" /D "NDEBUG" /D "_WINDOWS" /D "_AFXDLL" /D "_MBCS" /Yu"stdafx.h" /c
# ADD CPP /nologo /MT /W3 /GX /Ot /Oi /Ob1 /Gf /Gy /I ".." /D "WIN32" /D "NDEBUG" /D "_WINDOWS" /D "_MBCS" /Yu"stdafx.h" /Gs /Gs /c
CPP_PROJ=/nologo /MT /W3 /GX /Ot /Oi /Ob1 /Gf /Gy /I ".." /D "WIN32" /D\
 "NDEBUG" /D "_WINDOWS" /D "_MBCS" /Fp"$(INTDIR)/Prime95.pch" /Yu"stdafx.h"\
 /Fo"$(INTDIR)/" /Gs /Gs /c 
CPP_OBJS=.\Release/
CPP_SBRS=
# ADD BASE MTL /nologo /D "NDEBUG" /win32
# ADD MTL /nologo /D "NDEBUG" /win32
MTL_PROJ=/nologo /D "NDEBUG" /win32 
# ADD BASE RSC /l 0x409 /d "NDEBUG" /d "_AFXDLL"
# ADD RSC /l 0x409 /d "NDEBUG"
RSC_PROJ=/l 0x409 /fo"$(INTDIR)/Prime95.res" /d "NDEBUG" 
BSC32=bscmake.exe
# ADD BASE BSC32 /nologo
# ADD BSC32 /nologo
BSC32_FLAGS=/nologo /o"$(OUTDIR)/Prime95.bsc" 
BSC32_SBRS=
LINK32=link.exe
# ADD BASE LINK32 /nologo /subsystem:windows /machine:I386
# ADD LINK32 ../prime95/mult.obj ../prime95/mult1.obj ../prime95/mult1aux.obj ../prime95/mult2.obj ../prime95/mult2a.obj ../prime95/mult2aux.obj ../prime95/mult3.obj ../prime95/mult3a.obj ../prime95/mult3aux.obj ../prime95/mult3auq.obj ../prime95/mult4.obj ../prime95/mult4a.obj ../prime95/mult4b.obj ../prime95/mult1p.obj ../prime95/mult2p.obj ../prime95/mult2ap.obj ../prime95/mult3p.obj ../prime95/mult3ap.obj ../prime95/mult4p.obj ../prime95/mult4ap.obj ../prime95/mult4bp.obj ../prime95/mult4aux.obj ../prime95/mult4auq.obj ../prime95/mult3q.obj ../prime95/mult3aq.obj ../prime95/mult4q.obj ../prime95/mult4aq.obj ../prime95/mult4bq.obj ../prime95/ecmhelp.obj ../prime95/xmult1.obj ../prime95/xmult1ax.obj ../prime95/xmult2.obj ../prime95/xmult2a.obj ../prime95/xmult2ax.obj ../prime95/xmult3.obj ../prime95/xmult3a.obj ../prime95/xmult3ax.obj /nologo /subsystem:windows /map /machine:I386
# SUBTRACT LINK32 /debug
LINK32_FLAGS=../prime95/mult.obj ../prime95/mult1.obj ../prime95/mult1aux.obj\
 ../prime95/mult2.obj ../prime95/mult2a.obj ../prime95/mult2aux.obj\
 ../prime95/mult3.obj ../prime95/mult3a.obj ../prime95/mult3aux.obj\
 ../prime95/mult3auq.obj ../prime95/mult4.obj ../prime95/mult4a.obj\
 ../prime95/mult4b.obj ../prime95/mult1p.obj ../prime95/mult2p.obj\
 ../prime95/mult2ap.obj ../prime95/mult3p.obj ../prime95/mult3ap.obj\
 ../prime95/mult4p.obj ../prime95/mult4ap.obj ../prime95/mult4bp.obj\
 ../prime95/mult4aux.obj ../prime95/mult4auq.obj ../prime95/mult3q.obj\
 ../prime95/mult3aq.obj ../prime95/mult4q.obj ../prime95/mult4aq.obj\
 ../prime95/mult4bq.obj ../prime95/ecmhelp.obj ../prime95/xmult1.obj\
 ../prime95/xmult1ax.obj ../prime95/xmult2.obj ../prime95/xmult2a.obj\
 ../prime95/xmult2ax.obj ../prime95/xmult3.obj ../prime95/xmult3a.obj\
 ../prime95/xmult3ax.obj /nologo /subsystem:windows /incremental:no\
 /pdb:"$(OUTDIR)/Prime95.pdb" /map:"$(INTDIR)/Prime95.map" /machine:I386\
 /out:"$(OUTDIR)/Prime95.exe" 
LINK32_OBJS= \
	"$(INTDIR)/AffinityDlg.obj" \
	"$(INTDIR)/PreferencesDlg.obj" \
	"$(INTDIR)/Prime95.obj" \
	"$(INTDIR)/Prime95Doc.obj" \
	"$(INTDIR)/StdAfx.obj" \
	"$(INTDIR)/Priority.obj" \
	"$(INTDIR)/EditDropFiles.obj" \
	"$(INTDIR)/CpuDlg.obj" \
	"$(INTDIR)/Prime95View.obj" \
	"$(INTDIR)/TestDlg.obj" \
	"$(INTDIR)/MainFrm.obj" \
	"$(INTDIR)/ContentsDlg.obj" \
	"$(INTDIR)/Prime95.res"

"$(OUTDIR)\Prime95.exe" : "$(OUTDIR)" $(DEF_FILE) $(LINK32_OBJS)
    $(LINK32) @<<
  $(LINK32_FLAGS) $(LINK32_OBJS)
<<

!ELSEIF  "$(CFG)" == "Prime95 - Win32 Debug"

# PROP BASE Use_MFC 6
# PROP BASE Use_Debug_Libraries 1
# PROP BASE Output_Dir "Debug"
# PROP BASE Intermediate_Dir "Debug"
# PROP BASE Target_Dir ""
# PROP Use_MFC 6
# PROP Use_Debug_Libraries 1
# PROP Output_Dir "Debug"
# PROP Intermediate_Dir "Debug"
# PROP Target_Dir ""
OUTDIR=.\Debug
INTDIR=.\Debug

ALL : "$(OUTDIR)\Prime95.exe" "$(OUTDIR)\Prime95.hlp"

CLEAN : 
	-@erase ".\Debug\vc40.pdb"
	-@erase ".\Debug\Prime95.pch"
	-@erase ".\Debug\vc40.idb"
	-@erase ".\Debug\Prime95.exe"
	-@erase ".\Debug\Priority.obj"
	-@erase ".\Debug\CpuDlg.obj"
	-@erase ".\Debug\PreferencesDlg.obj"
	-@erase ".\Debug\TestDlg.obj"
	-@erase ".\Debug\Prime95.obj"
	-@erase ".\Debug\ContentsDlg.obj"
	-@erase ".\Debug\MainFrm.obj"
	-@erase ".\Debug\EditDropFiles.obj"
	-@erase ".\Debug\Prime95Doc.obj"
	-@erase ".\Debug\Prime95View.obj"
	-@erase ".\Debug\StdAfx.obj"
	-@erase ".\Debug\AffinityDlg.obj"
	-@erase ".\Debug\Prime95.res"
	-@erase ".\Debug\Prime95.ilk"
	-@erase ".\Debug\Prime95.pdb"
	-@erase ".\Debug\Prime95.hlp"

"$(OUTDIR)" :
    if not exist "$(OUTDIR)/$(NULL)" mkdir "$(OUTDIR)"

# ADD BASE CPP /nologo /MDd /W3 /Gm /GX /Zi /Od /D "WIN32" /D "_DEBUG" /D "_WINDOWS" /D "_AFXDLL" /D "_MBCS" /Yu"stdafx.h" /c
# ADD CPP /nologo /MDd /W3 /Gm /GX /Zi /Ot /Oi /I ".." /D "WIN32" /D "_DEBUG" /D "_WINDOWS" /D "_AFXDLL" /D "_MBCS" /D "GDEBUG" /Yu"stdafx.h" /c
CPP_PROJ=/nologo /MDd /W3 /Gm /GX /Zi /Ot /Oi /I ".." /D "WIN32" /D "_DEBUG" /D\
 "_WINDOWS" /D "_AFXDLL" /D "_MBCS" /D "GDEBUG" /Fp"$(INTDIR)/Prime95.pch"\
 /Yu"stdafx.h" /Fo"$(INTDIR)/" /Fd"$(INTDIR)/" /c 
CPP_OBJS=.\Debug/
CPP_SBRS=
# ADD BASE MTL /nologo /D "_DEBUG" /win32
# ADD MTL /nologo /D "_DEBUG" /win32
MTL_PROJ=/nologo /D "_DEBUG" /win32 
# ADD BASE RSC /l 0x409 /d "_DEBUG" /d "_AFXDLL"
# ADD RSC /l 0x409 /d "_DEBUG" /d "_AFXDLL"
RSC_PROJ=/l 0x409 /fo"$(INTDIR)/Prime95.res" /d "_DEBUG" /d "_AFXDLL" 
BSC32=bscmake.exe
# ADD BASE BSC32 /nologo
# ADD BSC32 /nologo
BSC32_FLAGS=/nologo /o"$(OUTDIR)/Prime95.bsc" 
BSC32_SBRS=
LINK32=link.exe
# ADD BASE LINK32 /nologo /subsystem:windows /debug /machine:I386
# ADD LINK32 ../prime95/mult.obj ../prime95/mult1.obj ../prime95/mult1aux.obj ../prime95/mult2.obj ../prime95/mult2a.obj ../prime95/mult2aux.obj ../prime95/mult3.obj ../prime95/mult3a.obj ../prime95/mult3aux.obj ../prime95/mult3auq.obj ../prime95/mult4.obj ../prime95/mult4a.obj ../prime95/mult4b.obj ../prime95/mult1p.obj ../prime95/mult2p.obj ../prime95/mult2ap.obj ../prime95/mult3p.obj ../prime95/mult3ap.obj ../prime95/mult4p.obj ../prime95/mult4ap.obj ../prime95/mult4bp.obj ../prime95/mult4aux.obj ../prime95/mult4auq.obj ../prime95/mult3q.obj ../prime95/mult3aq.obj ../prime95/mult4q.obj ../prime95/mult4aq.obj ../prime95/mult4bq.obj ../prime95/ecmhelp.obj ../prime95/xmult1.obj ../prime95/xmult1ax.obj ../prime95/xmult2.obj ../prime95/xmult2a.obj ../prime95/xmult2ax.obj ../prime95/xmult3.obj ../prime95/xmult3a.obj ../prime95/xmult3ax.obj /nologo /subsystem:windows /debug /machine:I386
LINK32_FLAGS=../prime95/mult.obj ../prime95/mult1.obj ../prime95/mult1aux.obj\
 ../prime95/mult2.obj ../prime95/mult2a.obj ../prime95/mult2aux.obj\
 ../prime95/mult3.obj ../prime95/mult3a.obj ../prime95/mult3aux.obj\
 ../prime95/mult3auq.obj ../prime95/mult4.obj ../prime95/mult4a.obj\
 ../prime95/mult4b.obj ../prime95/mult1p.obj ../prime95/mult2p.obj\
 ../prime95/mult2ap.obj ../prime95/mult3p.obj ../prime95/mult3ap.obj\
 ../prime95/mult4p.obj ../prime95/mult4ap.obj ../prime95/mult4bp.obj\
 ../prime95/mult4aux.obj ../prime95/mult4auq.obj ../prime95/mult3q.obj\
 ../prime95/mult3aq.obj ../prime95/mult4q.obj ../prime95/mult4aq.obj\
 ../prime95/mult4bq.obj ../prime95/ecmhelp.obj ../prime95/xmult1.obj\
 ../prime95/xmult1ax.obj ../prime95/xmult2.obj ../prime95/xmult2a.obj\
 ../prime95/xmult2ax.obj ../prime95/xmult3.obj ../prime95/xmult3a.obj\
 ../prime95/xmult3ax.obj /nologo /subsystem:windows /incremental:yes\
 /pdb:"$(OUTDIR)/Prime95.pdb" /debug /machine:I386 /out:"$(OUTDIR)/Prime95.exe" 
LINK32_OBJS= \
	"$(INTDIR)/Priority.obj" \
	"$(INTDIR)/CpuDlg.obj" \
	"$(INTDIR)/PreferencesDlg.obj" \
	"$(INTDIR)/TestDlg.obj" \
	"$(INTDIR)/Prime95.obj" \
	"$(INTDIR)/ContentsDlg.obj" \
	"$(INTDIR)/MainFrm.obj" \
	"$(INTDIR)/EditDropFiles.obj" \
	"$(INTDIR)/Prime95Doc.obj" \
	"$(INTDIR)/Prime95View.obj" \
	"$(INTDIR)/StdAfx.obj" \
	"$(INTDIR)/AffinityDlg.obj" \
	"$(INTDIR)/Prime95.res"

"$(OUTDIR)\Prime95.exe" : "$(OUTDIR)" $(DEF_FILE) $(LINK32_OBJS)
    $(LINK32) @<<
  $(LINK32_FLAGS) $(LINK32_OBJS)
<<

!ENDIF 

.c{$(CPP_OBJS)}.obj:
   $(CPP) $(CPP_PROJ) $<  

.cpp{$(CPP_OBJS)}.obj:
   $(CPP) $(CPP_PROJ) $<  

.cxx{$(CPP_OBJS)}.obj:
   $(CPP) $(CPP_PROJ) $<  

.c{$(CPP_SBRS)}.sbr:
   $(CPP) $(CPP_PROJ) $<  

.cpp{$(CPP_SBRS)}.sbr:
   $(CPP) $(CPP_PROJ) $<  

.cxx{$(CPP_SBRS)}.sbr:
   $(CPP) $(CPP_PROJ) $<  

################################################################################
# Begin Target

# Name "Prime95 - Win32 Release"
# Name "Prime95 - Win32 Debug"

!IF  "$(CFG)" == "Prime95 - Win32 Release"

!ELSEIF  "$(CFG)" == "Prime95 - Win32 Debug"

!ENDIF 

################################################################################
# Begin Source File

SOURCE=.\ReadMe.txt

!IF  "$(CFG)" == "Prime95 - Win32 Release"

!ELSEIF  "$(CFG)" == "Prime95 - Win32 Debug"

!ENDIF 

# End Source File
################################################################################
# Begin Source File

SOURCE=.\Prime95.cpp

!IF  "$(CFG)" == "Prime95 - Win32 Release"

DEP_CPP_PRIME=\
	".\StdAfx.h"\
	".\Prime95.h"\
	".\MainFrm.h"\
	".\Prime95Doc.h"\
	".\Prime95View.h"\
	

"$(INTDIR)\Prime95.obj" : $(SOURCE) $(DEP_CPP_PRIME) "$(INTDIR)"\
 "$(INTDIR)\Prime95.pch"


!ELSEIF  "$(CFG)" == "Prime95 - Win32 Debug"

DEP_CPP_PRIME=\
	".\StdAfx.h"\
	".\Prime95.h"\
	".\MainFrm.h"\
	".\Prime95Doc.h"\
	".\Prime95View.h"\
	".\EditDropFiles.h"\
	".\..\cpuid.h"\
	".\..\speed.h"\
	".\..\gwnum.h"\
	".\..\prp.h"\
	".\..\prp95.h"\
	

"$(INTDIR)\Prime95.obj" : $(SOURCE) $(DEP_CPP_PRIME) "$(INTDIR)"\
 "$(INTDIR)\Prime95.pch"


!ENDIF 

# End Source File
################################################################################
# Begin Source File

SOURCE=.\StdAfx.cpp
DEP_CPP_STDAF=\
	".\StdAfx.h"\
	

!IF  "$(CFG)" == "Prime95 - Win32 Release"

# ADD CPP /Yc"stdafx.h"

BuildCmds= \
	$(CPP) /nologo /MT /W3 /GX /Ot /Oi /Ob1 /Gf /Gy /I ".." /D "WIN32" /D "NDEBUG"\
 /D "_WINDOWS" /D "_MBCS" /Fp"$(INTDIR)/Prime95.pch" /Yc"stdafx.h"\
 /Fo"$(INTDIR)/" /Gs /Gs /c $(SOURCE) \
	

"$(INTDIR)\StdAfx.obj" : $(SOURCE) $(DEP_CPP_STDAF) "$(INTDIR)"
   $(BuildCmds)

"$(INTDIR)\Prime95.pch" : $(SOURCE) $(DEP_CPP_STDAF) "$(INTDIR)"
   $(BuildCmds)

!ELSEIF  "$(CFG)" == "Prime95 - Win32 Debug"

# ADD CPP /Yc"stdafx.h"

BuildCmds= \
	$(CPP) /nologo /MDd /W3 /Gm /GX /Zi /Ot /Oi /I ".." /D "WIN32" /D "_DEBUG" /D\
 "_WINDOWS" /D "_AFXDLL" /D "_MBCS" /D "GDEBUG" /Fp"$(INTDIR)/Prime95.pch"\
 /Yc"stdafx.h" /Fo"$(INTDIR)/" /Fd"$(INTDIR)/" /c $(SOURCE) \
	

"$(INTDIR)\StdAfx.obj" : $(SOURCE) $(DEP_CPP_STDAF) "$(INTDIR)"
   $(BuildCmds)

"$(INTDIR)\Prime95.pch" : $(SOURCE) $(DEP_CPP_STDAF) "$(INTDIR)"
   $(BuildCmds)

!ENDIF 

# End Source File
################################################################################
# Begin Source File

SOURCE=.\MainFrm.cpp

!IF  "$(CFG)" == "Prime95 - Win32 Release"

DEP_CPP_MAINF=\
	".\StdAfx.h"\
	".\Prime95.h"\
	".\MainFrm.h"\
	

"$(INTDIR)\MainFrm.obj" : $(SOURCE) $(DEP_CPP_MAINF) "$(INTDIR)"\
 "$(INTDIR)\Prime95.pch"


!ELSEIF  "$(CFG)" == "Prime95 - Win32 Debug"

DEP_CPP_MAINF=\
	".\StdAfx.h"\
	".\Prime95.h"\
	".\MainFrm.h"\
	".\EditDropFiles.h"\
	".\..\cpuid.h"\
	".\..\speed.h"\
	".\..\gwnum.h"\
	".\..\prp.h"\
	".\..\prp95.h"\
	

"$(INTDIR)\MainFrm.obj" : $(SOURCE) $(DEP_CPP_MAINF) "$(INTDIR)"\
 "$(INTDIR)\Prime95.pch"


!ENDIF 

# End Source File
################################################################################
# Begin Source File

SOURCE=.\Prime95Doc.cpp

!IF  "$(CFG)" == "Prime95 - Win32 Release"

DEP_CPP_PRIME9=\
	".\StdAfx.h"\
	".\MainFrm.h"\
	".\Prime95.h"\
	".\Prime95Doc.h"\
	".\AffinityDlg.h"\
	".\ContentsDlg.h"\
	".\CpuDlg.h"\
	".\PreferencesDlg.h"\
	".\Priority.h"\
	".\TestDlg.h"\
	{$(INCLUDE)}"\sys\Timeb.h"\
	".\..\cpuid.c"\
	".\..\speed.c"\
	".\..\giants.h"\
	".\..\giants.c"\
	".\..\gwnum.c"\
	".\..\prp.c"\
	".\..\prp95.c"\
	

"$(INTDIR)\Prime95Doc.obj" : $(SOURCE) $(DEP_CPP_PRIME9) "$(INTDIR)"\
 "$(INTDIR)\Prime95.pch"


!ELSEIF  "$(CFG)" == "Prime95 - Win32 Debug"

DEP_CPP_PRIME9=\
	".\StdAfx.h"\
	".\MainFrm.h"\
	".\Prime95.h"\
	".\Prime95Doc.h"\
	".\AffinityDlg.h"\
	".\ContentsDlg.h"\
	".\CpuDlg.h"\
	".\PreferencesDlg.h"\
	".\Priority.h"\
	".\TestDlg.h"\
	{$(INCLUDE)}"\sys\Timeb.h"\
	".\..\cpuid.c"\
	".\..\speed.c"\
	".\..\giants.h"\
	".\..\giants.c"\
	".\..\gwnum.c"\
	".\..\prp.c"\
	".\..\prp95.c"\
	".\EditDropFiles.h"\
	".\..\cpuid.h"\
	".\..\speed.h"\
	".\..\gwnum.h"\
	".\..\prp.h"\
	".\..\prp95.h"\
	

"$(INTDIR)\Prime95Doc.obj" : $(SOURCE) $(DEP_CPP_PRIME9) "$(INTDIR)"\
 "$(INTDIR)\Prime95.pch"


!ENDIF 

# End Source File
################################################################################
# Begin Source File

SOURCE=.\Prime95View.cpp

!IF  "$(CFG)" == "Prime95 - Win32 Release"

DEP_CPP_PRIME95=\
	".\StdAfx.h"\
	".\Prime95.h"\
	".\Prime95Doc.h"\
	".\Prime95View.h"\
	

"$(INTDIR)\Prime95View.obj" : $(SOURCE) $(DEP_CPP_PRIME95) "$(INTDIR)"\
 "$(INTDIR)\Prime95.pch"


!ELSEIF  "$(CFG)" == "Prime95 - Win32 Debug"

DEP_CPP_PRIME95=\
	".\StdAfx.h"\
	".\Prime95.h"\
	".\Prime95Doc.h"\
	".\Prime95View.h"\
	".\EditDropFiles.h"\
	".\..\cpuid.h"\
	".\..\speed.h"\
	".\..\gwnum.h"\
	".\..\prp.h"\
	".\..\prp95.h"\
	

"$(INTDIR)\Prime95View.obj" : $(SOURCE) $(DEP_CPP_PRIME95) "$(INTDIR)"\
 "$(INTDIR)\Prime95.pch"


!ENDIF 

# End Source File
################################################################################
# Begin Source File

SOURCE=.\Prime95.rc
DEP_RSC_PRIME95_=\
	".\res\Prime95.ico"\
	".\res\yellow_i.ico"\
	".\res\cursor1.cur"\
	".\res\cursor2.cur"\
	".\res\Prime95.rc2"\
	

"$(INTDIR)\Prime95.res" : $(SOURCE) $(DEP_RSC_PRIME95_) "$(INTDIR)"
   $(RSC) $(RSC_PROJ) $(SOURCE)


# End Source File
################################################################################
# Begin Source File

SOURCE=.\hlp\Prime95.hpj

!IF  "$(CFG)" == "Prime95 - Win32 Release"

# Begin Custom Build - Making help file...
OutDir=.\Release
ProjDir=.
TargetName=Prime95
InputPath=.\hlp\Prime95.hpj

"$(OutDir)\$(TargetName).hlp" : $(SOURCE) "$(INTDIR)" "$(OUTDIR)"
   "$(ProjDir)\makehelp.bat"

# End Custom Build

!ELSEIF  "$(CFG)" == "Prime95 - Win32 Debug"

# Begin Custom Build - Making help file...
OutDir=.\Debug
ProjDir=.
TargetName=Prime95
InputPath=.\hlp\Prime95.hpj

"$(OutDir)\$(TargetName).hlp" : $(SOURCE) "$(INTDIR)" "$(OUTDIR)"
   "$(ProjDir)\makehelp.bat"

# End Custom Build

!ENDIF 

# End Source File
################################################################################
# Begin Source File

SOURCE=.\CpuDlg.cpp

!IF  "$(CFG)" == "Prime95 - Win32 Release"

DEP_CPP_CPUDL=\
	".\StdAfx.h"\
	".\Prime95.h"\
	".\CpuDlg.h"\
	

"$(INTDIR)\CpuDlg.obj" : $(SOURCE) $(DEP_CPP_CPUDL) "$(INTDIR)"\
 "$(INTDIR)\Prime95.pch"


!ELSEIF  "$(CFG)" == "Prime95 - Win32 Debug"

DEP_CPP_CPUDL=\
	".\StdAfx.h"\
	".\Prime95.h"\
	".\CpuDlg.h"\
	".\EditDropFiles.h"\
	".\..\cpuid.h"\
	".\..\speed.h"\
	".\..\gwnum.h"\
	".\..\prp.h"\
	".\..\prp95.h"\
	

"$(INTDIR)\CpuDlg.obj" : $(SOURCE) $(DEP_CPP_CPUDL) "$(INTDIR)"\
 "$(INTDIR)\Prime95.pch"


!ENDIF 

# End Source File
################################################################################
# Begin Source File

SOURCE=.\PreferencesDlg.cpp

!IF  "$(CFG)" == "Prime95 - Win32 Release"

DEP_CPP_PREFE=\
	".\StdAfx.h"\
	".\Prime95.h"\
	".\PreferencesDlg.h"\
	

"$(INTDIR)\PreferencesDlg.obj" : $(SOURCE) $(DEP_CPP_PREFE) "$(INTDIR)"\
 "$(INTDIR)\Prime95.pch"


!ELSEIF  "$(CFG)" == "Prime95 - Win32 Debug"

DEP_CPP_PREFE=\
	".\StdAfx.h"\
	".\Prime95.h"\
	".\PreferencesDlg.h"\
	".\EditDropFiles.h"\
	".\..\cpuid.h"\
	".\..\speed.h"\
	".\..\gwnum.h"\
	".\..\prp.h"\
	".\..\prp95.h"\
	

"$(INTDIR)\PreferencesDlg.obj" : $(SOURCE) $(DEP_CPP_PREFE) "$(INTDIR)"\
 "$(INTDIR)\Prime95.pch"


!ENDIF 

# End Source File
################################################################################
# Begin Source File

SOURCE=.\TestDlg.cpp

!IF  "$(CFG)" == "Prime95 - Win32 Release"

DEP_CPP_TESTD=\
	".\StdAfx.h"\
	".\Prime95.h"\
	".\TestDlg.h"\
	

"$(INTDIR)\TestDlg.obj" : $(SOURCE) $(DEP_CPP_TESTD) "$(INTDIR)"\
 "$(INTDIR)\Prime95.pch"


!ELSEIF  "$(CFG)" == "Prime95 - Win32 Debug"

DEP_CPP_TESTD=\
	".\StdAfx.h"\
	".\Prime95.h"\
	".\TestDlg.h"\
	".\EditDropFiles.h"\
	".\..\cpuid.h"\
	".\..\speed.h"\
	".\..\gwnum.h"\
	".\..\prp.h"\
	".\..\prp95.h"\
	

"$(INTDIR)\TestDlg.obj" : $(SOURCE) $(DEP_CPP_TESTD) "$(INTDIR)"\
 "$(INTDIR)\Prime95.pch"


!ENDIF 

# End Source File
################################################################################
# Begin Source File

SOURCE=.\Priority.cpp

!IF  "$(CFG)" == "Prime95 - Win32 Release"

DEP_CPP_PRIOR=\
	".\StdAfx.h"\
	".\Prime95.h"\
	".\Priority.h"\
	

"$(INTDIR)\Priority.obj" : $(SOURCE) $(DEP_CPP_PRIOR) "$(INTDIR)"\
 "$(INTDIR)\Prime95.pch"


!ELSEIF  "$(CFG)" == "Prime95 - Win32 Debug"

DEP_CPP_PRIOR=\
	".\StdAfx.h"\
	".\Prime95.h"\
	".\Priority.h"\
	".\EditDropFiles.h"\
	".\..\cpuid.h"\
	".\..\speed.h"\
	".\..\gwnum.h"\
	".\..\prp.h"\
	".\..\prp95.h"\
	

"$(INTDIR)\Priority.obj" : $(SOURCE) $(DEP_CPP_PRIOR) "$(INTDIR)"\
 "$(INTDIR)\Prime95.pch"


!ENDIF 

# End Source File
################################################################################
# Begin Source File

SOURCE=.\AffinityDlg.cpp

!IF  "$(CFG)" == "Prime95 - Win32 Release"

DEP_CPP_AFFIN=\
	".\StdAfx.h"\
	".\Prime95.h"\
	".\AffinityDlg.h"\
	

"$(INTDIR)\AffinityDlg.obj" : $(SOURCE) $(DEP_CPP_AFFIN) "$(INTDIR)"\
 "$(INTDIR)\Prime95.pch"


!ELSEIF  "$(CFG)" == "Prime95 - Win32 Debug"

DEP_CPP_AFFIN=\
	".\StdAfx.h"\
	".\Prime95.h"\
	".\AffinityDlg.h"\
	".\EditDropFiles.h"\
	".\..\cpuid.h"\
	".\..\speed.h"\
	".\..\gwnum.h"\
	".\..\prp.h"\
	".\..\prp95.h"\
	

"$(INTDIR)\AffinityDlg.obj" : $(SOURCE) $(DEP_CPP_AFFIN) "$(INTDIR)"\
 "$(INTDIR)\Prime95.pch"


!ENDIF 

# End Source File
################################################################################
# Begin Source File

SOURCE=.\EditDropFiles.cpp

!IF  "$(CFG)" == "Prime95 - Win32 Release"

DEP_CPP_EDITD=\
	".\StdAfx.h"\
	".\Prime95.h"\
	

"$(INTDIR)\EditDropFiles.obj" : $(SOURCE) $(DEP_CPP_EDITD) "$(INTDIR)"\
 "$(INTDIR)\Prime95.pch"


!ELSEIF  "$(CFG)" == "Prime95 - Win32 Debug"

DEP_CPP_EDITD=\
	".\StdAfx.h"\
	".\Prime95.h"\
	".\EditDropFiles.h"\
	".\..\cpuid.h"\
	".\..\speed.h"\
	".\..\gwnum.h"\
	".\..\prp.h"\
	".\..\prp95.h"\
	

"$(INTDIR)\EditDropFiles.obj" : $(SOURCE) $(DEP_CPP_EDITD) "$(INTDIR)"\
 "$(INTDIR)\Prime95.pch"


!ENDIF 

# End Source File
################################################################################
# Begin Source File

SOURCE=.\ContentsDlg.cpp

!IF  "$(CFG)" == "Prime95 - Win32 Release"

DEP_CPP_CONTE=\
	".\StdAfx.h"\
	".\Prime95.h"\
	".\ContentsDlg.h"\
	

"$(INTDIR)\ContentsDlg.obj" : $(SOURCE) $(DEP_CPP_CONTE) "$(INTDIR)"\
 "$(INTDIR)\Prime95.pch"


!ELSEIF  "$(CFG)" == "Prime95 - Win32 Debug"

DEP_CPP_CONTE=\
	".\StdAfx.h"\
	".\Prime95.h"\
	".\ContentsDlg.h"\
	".\EditDropFiles.h"\
	".\..\cpuid.h"\
	".\..\speed.h"\
	".\..\gwnum.h"\
	".\..\prp.h"\
	".\..\prp95.h"\
	

"$(INTDIR)\ContentsDlg.obj" : $(SOURCE) $(DEP_CPP_CONTE) "$(INTDIR)"\
 "$(INTDIR)\Prime95.pch"


!ENDIF 

# End Source File
# End Target
# End Project
################################################################################
