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
	-@erase ".\Release\Prime95View.obj"
	-@erase ".\Release\Prime95.pch"
	-@erase ".\Release\Prime95.obj"
	-@erase ".\Release\AffinityDlg.obj"
	-@erase ".\Release\Priority.obj"
	-@erase ".\Release\EcmDlg.obj"
	-@erase ".\Release\Prime95Doc.obj"
	-@erase ".\Release\PreferencesDlg.obj"
	-@erase ".\Release\Pminus1Dlg.obj"
	-@erase ".\Release\PrimenetDlg.obj"
	-@erase ".\Release\StdAfx.obj"
	-@erase ".\Release\ServerDlg.obj"
	-@erase ".\Release\MainFrm.obj"
	-@erase ".\Release\VacationDlg.obj"
	-@erase ".\Release\UserDlg.obj"
	-@erase ".\Release\ManualCommDlg.obj"
	-@erase ".\Release\CpuDlg.obj"
	-@erase ".\Release\Prime95.res"
	-@erase ".\Release\Prime95.hlp"

"$(OUTDIR)" :
    if not exist "$(OUTDIR)/$(NULL)" mkdir "$(OUTDIR)"

# ADD BASE CPP /nologo /MD /W3 /GX /O2 /D "WIN32" /D "NDEBUG" /D "_WINDOWS" /D "_AFXDLL" /D "_MBCS" /Yu"stdafx.h" /c
# ADD CPP /nologo /MT /W3 /GX /Oi /Os /Ob1 /Gf /Gy /I ".." /D "WIN32" /D "NDEBUG" /D "_WINDOWS" /D "_MBCS" /Yu"stdafx.h" /Gs /Gs /c
# SUBTRACT CPP /Og /Fr
CPP_PROJ=/nologo /MT /W3 /GX /Oi /Os /Ob1 /Gf /Gy /I ".." /D "WIN32" /D\
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
# ADD LINK32 ..\prime95\mult.obj mult.obj /nologo /subsystem:windows /machine:I386
# SUBTRACT LINK32 /map /debug
LINK32_FLAGS=..\prime95\mult.obj mult.obj /nologo /subsystem:windows\
 /incremental:no /pdb:"$(OUTDIR)/Prime95.pdb" /machine:I386\
 /out:"$(OUTDIR)/Prime95.exe" 
LINK32_OBJS= \
	"$(INTDIR)/Prime95View.obj" \
	"$(INTDIR)/Prime95.obj" \
	"$(INTDIR)/AffinityDlg.obj" \
	"$(INTDIR)/Priority.obj" \
	"$(INTDIR)/EcmDlg.obj" \
	"$(INTDIR)/Prime95Doc.obj" \
	"$(INTDIR)/PreferencesDlg.obj" \
	"$(INTDIR)/Pminus1Dlg.obj" \
	"$(INTDIR)/PrimenetDlg.obj" \
	"$(INTDIR)/StdAfx.obj" \
	"$(INTDIR)/ServerDlg.obj" \
	"$(INTDIR)/MainFrm.obj" \
	"$(INTDIR)/VacationDlg.obj" \
	"$(INTDIR)/UserDlg.obj" \
	"$(INTDIR)/ManualCommDlg.obj" \
	"$(INTDIR)/CpuDlg.obj" \
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
	-@erase ".\Debug\vc40.idb"
	-@erase ".\Debug\Prime95.exe"
	-@erase ".\Debug\Pminus1Dlg.obj"
	-@erase ".\Debug\Prime95.pch"
	-@erase ".\Debug\Prime95View.obj"
	-@erase ".\Debug\AffinityDlg.obj"
	-@erase ".\Debug\CpuDlg.obj"
	-@erase ".\Debug\MainFrm.obj"
	-@erase ".\Debug\ManualCommDlg.obj"
	-@erase ".\Debug\EcmDlg.obj"
	-@erase ".\Debug\ServerDlg.obj"
	-@erase ".\Debug\Prime95Doc.obj"
	-@erase ".\Debug\StdAfx.obj"
	-@erase ".\Debug\Prime95.obj"
	-@erase ".\Debug\PrimenetDlg.obj"
	-@erase ".\Debug\Priority.obj"
	-@erase ".\Debug\PreferencesDlg.obj"
	-@erase ".\Debug\UserDlg.obj"
	-@erase ".\Debug\VacationDlg.obj"
	-@erase ".\Debug\Prime95.res"
	-@erase ".\Debug\Prime95.ilk"
	-@erase ".\Debug\Prime95.pdb"
	-@erase ".\Debug\Prime95.hlp"

"$(OUTDIR)" :
    if not exist "$(OUTDIR)/$(NULL)" mkdir "$(OUTDIR)"

# ADD BASE CPP /nologo /MDd /W3 /Gm /GX /Zi /Od /D "WIN32" /D "_DEBUG" /D "_WINDOWS" /D "_AFXDLL" /D "_MBCS" /Yu"stdafx.h" /c
# ADD CPP /nologo /MDd /W3 /Gm /GX /Zi /Od /I ".." /D "WIN32" /D "_DEBUG" /D "_WINDOWS" /D "_AFXDLL" /D "_MBCS" /Yu"stdafx.h" /c
CPP_PROJ=/nologo /MDd /W3 /Gm /GX /Zi /Od /I ".." /D "WIN32" /D "_DEBUG" /D\
 "_WINDOWS" /D "_AFXDLL" /D "_MBCS" /Fp"$(INTDIR)/Prime95.pch" /Yu"stdafx.h"\
 /Fo"$(INTDIR)/" /Fd"$(INTDIR)/" /c 
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
# ADD LINK32 ..\prime95\mult.obj mult.obj /nologo /subsystem:windows /debug /machine:I386
LINK32_FLAGS=..\prime95\mult.obj mult.obj /nologo /subsystem:windows\
 /incremental:yes /pdb:"$(OUTDIR)/Prime95.pdb" /debug /machine:I386\
 /out:"$(OUTDIR)/Prime95.exe" 
LINK32_OBJS= \
	"$(INTDIR)/Pminus1Dlg.obj" \
	"$(INTDIR)/Prime95View.obj" \
	"$(INTDIR)/AffinityDlg.obj" \
	"$(INTDIR)/CpuDlg.obj" \
	"$(INTDIR)/MainFrm.obj" \
	"$(INTDIR)/ManualCommDlg.obj" \
	"$(INTDIR)/EcmDlg.obj" \
	"$(INTDIR)/ServerDlg.obj" \
	"$(INTDIR)/Prime95Doc.obj" \
	"$(INTDIR)/StdAfx.obj" \
	"$(INTDIR)/Prime95.obj" \
	"$(INTDIR)/PrimenetDlg.obj" \
	"$(INTDIR)/Priority.obj" \
	"$(INTDIR)/PreferencesDlg.obj" \
	"$(INTDIR)/UserDlg.obj" \
	"$(INTDIR)/VacationDlg.obj" \
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

SOURCE=.\Prime95.cpp
DEP_CPP_PRIME=\
	".\StdAfx.h"\
	".\Prime95.h"\
	".\MainFrm.h"\
	".\Prime95Doc.h"\
	".\Prime95View.h"\
	".\..\cpuid.h"\
	".\..\speed.h"\
	".\..\gwnum.h"\
	".\..\commona.h"\
	".\..\commonc.h"\
	".\..\comm95a.h"\
	".\..\comm95c.h"\
	".\..\primenet.h"\
	

"$(INTDIR)\Prime95.obj" : $(SOURCE) $(DEP_CPP_PRIME) "$(INTDIR)"\
 "$(INTDIR)\Prime95.pch"


# End Source File
################################################################################
# Begin Source File

SOURCE=.\StdAfx.cpp
DEP_CPP_STDAF=\
	".\StdAfx.h"\
	

!IF  "$(CFG)" == "Prime95 - Win32 Release"

# ADD CPP /Yc"stdafx.h"

BuildCmds= \
	$(CPP) /nologo /MT /W3 /GX /Oi /Os /Ob1 /Gf /Gy /I ".." /D "WIN32" /D "NDEBUG"\
 /D "_WINDOWS" /D "_MBCS" /Fp"$(INTDIR)/Prime95.pch" /Yc"stdafx.h"\
 /Fo"$(INTDIR)/" /Gs /Gs /c $(SOURCE) \
	

"$(INTDIR)\StdAfx.obj" : $(SOURCE) $(DEP_CPP_STDAF) "$(INTDIR)"
   $(BuildCmds)

"$(INTDIR)\Prime95.pch" : $(SOURCE) $(DEP_CPP_STDAF) "$(INTDIR)"
   $(BuildCmds)

!ELSEIF  "$(CFG)" == "Prime95 - Win32 Debug"

# ADD CPP /Yc"stdafx.h"

BuildCmds= \
	$(CPP) /nologo /MDd /W3 /Gm /GX /Zi /Od /I ".." /D "WIN32" /D "_DEBUG" /D\
 "_WINDOWS" /D "_AFXDLL" /D "_MBCS" /Fp"$(INTDIR)/Prime95.pch" /Yc"stdafx.h"\
 /Fo"$(INTDIR)/" /Fd"$(INTDIR)/" /c $(SOURCE) \
	

"$(INTDIR)\StdAfx.obj" : $(SOURCE) $(DEP_CPP_STDAF) "$(INTDIR)"
   $(BuildCmds)

"$(INTDIR)\Prime95.pch" : $(SOURCE) $(DEP_CPP_STDAF) "$(INTDIR)"
   $(BuildCmds)

!ENDIF 

# End Source File
################################################################################
# Begin Source File

SOURCE=.\MainFrm.cpp
DEP_CPP_MAINF=\
	".\StdAfx.h"\
	".\Prime95.h"\
	".\MainFrm.h"\
	".\..\cpuid.h"\
	".\..\speed.h"\
	".\..\gwnum.h"\
	".\..\commona.h"\
	".\..\commonc.h"\
	".\..\comm95a.h"\
	".\..\comm95c.h"\
	".\..\primenet.h"\
	

"$(INTDIR)\MainFrm.obj" : $(SOURCE) $(DEP_CPP_MAINF) "$(INTDIR)"\
 "$(INTDIR)\Prime95.pch"


# End Source File
################################################################################
# Begin Source File

SOURCE=.\Prime95Doc.cpp
DEP_CPP_PRIME9=\
	".\StdAfx.h"\
	".\MainFrm.h"\
	".\Prime95.h"\
	".\Prime95Doc.h"\
	".\AffinityDlg.h"\
	".\CpuDlg.h"\
	".\EcmDlg.h"\
	".\ManualCommDlg.h"\
	".\Pminus1Dlg.h"\
	".\PreferencesDlg.h"\
	".\PrimenetDlg.h"\
	".\Priority.h"\
	".\ServerDlg.h"\
	".\UserDlg.h"\
	".\VacationDlg.h"\
	".\..\cpuid.c"\
	".\..\speed.c"\
	".\..\gwnum.c"\
	".\..\commona.c"\
	".\..\commonc.c"\
	".\..\comm95a.c"\
	".\..\comm95c.c"\
	".\..\cpuid.h"\
	".\..\speed.h"\
	".\..\gwnum.h"\
	".\..\commona.h"\
	".\..\commonc.h"\
	".\..\comm95a.h"\
	".\..\comm95c.h"\
	".\..\primenet.h"\
	".\..\security.c"\
	

"$(INTDIR)\Prime95Doc.obj" : $(SOURCE) $(DEP_CPP_PRIME9) "$(INTDIR)"\
 "$(INTDIR)\Prime95.pch"


# End Source File
################################################################################
# Begin Source File

SOURCE=.\Prime95View.cpp
DEP_CPP_PRIME95=\
	".\StdAfx.h"\
	".\Prime95.h"\
	".\Prime95Doc.h"\
	".\Prime95View.h"\
	".\..\cpuid.h"\
	".\..\speed.h"\
	".\..\gwnum.h"\
	".\..\commona.h"\
	".\..\commonc.h"\
	".\..\comm95a.h"\
	".\..\comm95c.h"\
	".\..\primenet.h"\
	

"$(INTDIR)\Prime95View.obj" : $(SOURCE) $(DEP_CPP_PRIME95) "$(INTDIR)"\
 "$(INTDIR)\Prime95.pch"


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
DEP_CPP_CPUDL=\
	".\StdAfx.h"\
	".\Prime95.h"\
	".\CpuDlg.h"\
	".\..\cpuid.h"\
	".\..\speed.h"\
	".\..\gwnum.h"\
	".\..\commona.h"\
	".\..\commonc.h"\
	".\..\comm95a.h"\
	".\..\comm95c.h"\
	".\..\primenet.h"\
	

"$(INTDIR)\CpuDlg.obj" : $(SOURCE) $(DEP_CPP_CPUDL) "$(INTDIR)"\
 "$(INTDIR)\Prime95.pch"


# End Source File
################################################################################
# Begin Source File

SOURCE=.\PreferencesDlg.cpp
DEP_CPP_PREFE=\
	".\StdAfx.h"\
	".\Prime95.h"\
	".\PreferencesDlg.h"\
	".\..\cpuid.h"\
	".\..\speed.h"\
	".\..\gwnum.h"\
	".\..\commona.h"\
	".\..\commonc.h"\
	".\..\comm95a.h"\
	".\..\comm95c.h"\
	".\..\primenet.h"\
	

"$(INTDIR)\PreferencesDlg.obj" : $(SOURCE) $(DEP_CPP_PREFE) "$(INTDIR)"\
 "$(INTDIR)\Prime95.pch"


# End Source File
################################################################################
# Begin Source File

SOURCE=.\UserDlg.cpp
DEP_CPP_USERD=\
	".\StdAfx.h"\
	".\Prime95.h"\
	".\UserDlg.h"\
	".\..\cpuid.h"\
	".\..\speed.h"\
	".\..\gwnum.h"\
	".\..\commona.h"\
	".\..\commonc.h"\
	".\..\comm95a.h"\
	".\..\comm95c.h"\
	".\..\primenet.h"\
	

"$(INTDIR)\UserDlg.obj" : $(SOURCE) $(DEP_CPP_USERD) "$(INTDIR)"\
 "$(INTDIR)\Prime95.pch"


# End Source File
################################################################################
# Begin Source File

SOURCE=.\PrimenetDlg.cpp
DEP_CPP_PRIMEN=\
	".\StdAfx.h"\
	".\Prime95.h"\
	".\PrimenetDlg.h"\
	".\..\cpuid.h"\
	".\..\speed.h"\
	".\..\gwnum.h"\
	".\..\commona.h"\
	".\..\commonc.h"\
	".\..\comm95a.h"\
	".\..\comm95c.h"\
	".\..\primenet.h"\
	

"$(INTDIR)\PrimenetDlg.obj" : $(SOURCE) $(DEP_CPP_PRIMEN) "$(INTDIR)"\
 "$(INTDIR)\Prime95.pch"


# End Source File
################################################################################
# Begin Source File

SOURCE=.\Priority.cpp
DEP_CPP_PRIOR=\
	".\StdAfx.h"\
	".\Prime95.h"\
	".\Priority.h"\
	".\..\cpuid.h"\
	".\..\speed.h"\
	".\..\gwnum.h"\
	".\..\commona.h"\
	".\..\commonc.h"\
	".\..\comm95a.h"\
	".\..\comm95c.h"\
	".\..\primenet.h"\
	

"$(INTDIR)\Priority.obj" : $(SOURCE) $(DEP_CPP_PRIOR) "$(INTDIR)"\
 "$(INTDIR)\Prime95.pch"


# End Source File
################################################################################
# Begin Source File

SOURCE=.\ServerDlg.cpp
DEP_CPP_SERVE=\
	".\StdAfx.h"\
	".\Prime95.h"\
	".\ServerDlg.h"\
	".\..\cpuid.h"\
	".\..\speed.h"\
	".\..\gwnum.h"\
	".\..\commona.h"\
	".\..\commonc.h"\
	".\..\comm95a.h"\
	".\..\comm95c.h"\
	".\..\primenet.h"\
	

"$(INTDIR)\ServerDlg.obj" : $(SOURCE) $(DEP_CPP_SERVE) "$(INTDIR)"\
 "$(INTDIR)\Prime95.pch"


# End Source File
################################################################################
# Begin Source File

SOURCE=.\VacationDlg.cpp
DEP_CPP_VACAT=\
	".\StdAfx.h"\
	".\Prime95.h"\
	".\VacationDlg.h"\
	".\..\cpuid.h"\
	".\..\speed.h"\
	".\..\gwnum.h"\
	".\..\commona.h"\
	".\..\commonc.h"\
	".\..\comm95a.h"\
	".\..\comm95c.h"\
	".\..\primenet.h"\
	

"$(INTDIR)\VacationDlg.obj" : $(SOURCE) $(DEP_CPP_VACAT) "$(INTDIR)"\
 "$(INTDIR)\Prime95.pch"


# End Source File
################################################################################
# Begin Source File

SOURCE=.\ManualCommDlg.cpp
DEP_CPP_MANUA=\
	".\StdAfx.h"\
	".\Prime95.h"\
	".\ManualCommDlg.h"\
	".\..\cpuid.h"\
	".\..\speed.h"\
	".\..\gwnum.h"\
	".\..\commona.h"\
	".\..\commonc.h"\
	".\..\comm95a.h"\
	".\..\comm95c.h"\
	".\..\primenet.h"\
	

"$(INTDIR)\ManualCommDlg.obj" : $(SOURCE) $(DEP_CPP_MANUA) "$(INTDIR)"\
 "$(INTDIR)\Prime95.pch"


# End Source File
################################################################################
# Begin Source File

SOURCE=.\AffinityDlg.cpp
DEP_CPP_AFFIN=\
	".\StdAfx.h"\
	".\Prime95.h"\
	".\AffinityDlg.h"\
	".\..\cpuid.h"\
	".\..\speed.h"\
	".\..\gwnum.h"\
	".\..\commona.h"\
	".\..\commonc.h"\
	".\..\comm95a.h"\
	".\..\comm95c.h"\
	".\..\primenet.h"\
	

"$(INTDIR)\AffinityDlg.obj" : $(SOURCE) $(DEP_CPP_AFFIN) "$(INTDIR)"\
 "$(INTDIR)\Prime95.pch"


# End Source File
################################################################################
# Begin Source File

SOURCE=.\EcmDlg.cpp
DEP_CPP_ECMDL=\
	".\StdAfx.h"\
	".\Prime95.h"\
	".\EcmDlg.h"\
	".\..\cpuid.h"\
	".\..\speed.h"\
	".\..\gwnum.h"\
	".\..\commona.h"\
	".\..\commonc.h"\
	".\..\comm95a.h"\
	".\..\comm95c.h"\
	".\..\primenet.h"\
	

"$(INTDIR)\EcmDlg.obj" : $(SOURCE) $(DEP_CPP_ECMDL) "$(INTDIR)"\
 "$(INTDIR)\Prime95.pch"


# End Source File
################################################################################
# Begin Source File

SOURCE=.\Pminus1Dlg.cpp
DEP_CPP_PMINU=\
	".\StdAfx.h"\
	".\Prime95.h"\
	".\Pminus1Dlg.h"\
	".\..\cpuid.h"\
	".\..\speed.h"\
	".\..\gwnum.h"\
	".\..\commona.h"\
	".\..\commonc.h"\
	".\..\comm95a.h"\
	".\..\comm95c.h"\
	".\..\primenet.h"\
	

"$(INTDIR)\Pminus1Dlg.obj" : $(SOURCE) $(DEP_CPP_PMINU) "$(INTDIR)"\
 "$(INTDIR)\Prime95.pch"


# End Source File
# End Target
# End Project
################################################################################
