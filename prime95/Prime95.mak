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
RSC=rc.exe
CPP=cl.exe
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
	-@erase ".\Release\ServerDlg.obj"
	-@erase ".\Release\Prime95.pch"
	-@erase ".\Release\AffinityDlg.obj"
	-@erase ".\Release\TestDlg.obj"
	-@erase ".\Release\Prime95.obj"
	-@erase ".\Release\Priority.obj"
	-@erase ".\Release\MainFrm.obj"
	-@erase ".\Release\VacationDlg.obj"
	-@erase ".\Release\Prime95Doc.obj"
	-@erase ".\Release\PrimenetDlg.obj"
	-@erase ".\Release\StdAfx.obj"
	-@erase ".\Release\PreferencesDlg.obj"
	-@erase ".\Release\ManualCommDlg.obj"
	-@erase ".\Release\Prime95View.obj"
	-@erase ".\Release\CpuDlg.obj"
	-@erase ".\Release\Password.obj"
	-@erase ".\Release\UserDlg.obj"
	-@erase ".\Release\TimeDlg.obj"
	-@erase ".\Release\Pminus1Dlg.obj"
	-@erase ".\Release\EcmDlg.obj"
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
# ADD LINK32 factor64.obj mult.obj mult1.obj mult2.obj mult2a.obj mult3.obj mult3a.obj mult4.obj mult4a.obj mult4b.obj mult1p.obj mult2p.obj mult2ap.obj mult3p.obj mult3ap.obj mult4p.obj mult4ap.obj mult4bp.obj mult4aux.obj ecmhelp.obj /nologo /subsystem:windows /map /machine:I386
# SUBTRACT LINK32 /debug
LINK32_FLAGS=factor64.obj mult.obj mult1.obj mult2.obj mult2a.obj mult3.obj\
 mult3a.obj mult4.obj mult4a.obj mult4b.obj mult1p.obj mult2p.obj mult2ap.obj\
 mult3p.obj mult3ap.obj mult4p.obj mult4ap.obj mult4bp.obj mult4aux.obj\
 ecmhelp.obj /nologo /subsystem:windows /incremental:no\
 /pdb:"$(OUTDIR)/Prime95.pdb" /map:"$(INTDIR)/Prime95.map" /machine:I386\
 /out:"$(OUTDIR)/Prime95.exe" 
LINK32_OBJS= \
	"$(INTDIR)/ServerDlg.obj" \
	"$(INTDIR)/AffinityDlg.obj" \
	"$(INTDIR)/TestDlg.obj" \
	"$(INTDIR)/Prime95.obj" \
	"$(INTDIR)/Priority.obj" \
	"$(INTDIR)/MainFrm.obj" \
	"$(INTDIR)/VacationDlg.obj" \
	"$(INTDIR)/Prime95Doc.obj" \
	"$(INTDIR)/PrimenetDlg.obj" \
	"$(INTDIR)/StdAfx.obj" \
	"$(INTDIR)/PreferencesDlg.obj" \
	"$(INTDIR)/ManualCommDlg.obj" \
	"$(INTDIR)/Prime95View.obj" \
	"$(INTDIR)/CpuDlg.obj" \
	"$(INTDIR)/Password.obj" \
	"$(INTDIR)/UserDlg.obj" \
	"$(INTDIR)/TimeDlg.obj" \
	"$(INTDIR)/Pminus1Dlg.obj" \
	"$(INTDIR)/EcmDlg.obj" \
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
	-@erase ".\Debug\ManualCommDlg.obj"
	-@erase ".\Debug\AffinityDlg.obj"
	-@erase ".\Debug\Pminus1Dlg.obj"
	-@erase ".\Debug\StdAfx.obj"
	-@erase ".\Debug\VacationDlg.obj"
	-@erase ".\Debug\TestDlg.obj"
	-@erase ".\Debug\PrimenetDlg.obj"
	-@erase ".\Debug\Prime95View.obj"
	-@erase ".\Debug\CpuDlg.obj"
	-@erase ".\Debug\Priority.obj"
	-@erase ".\Debug\MainFrm.obj"
	-@erase ".\Debug\ServerDlg.obj"
	-@erase ".\Debug\UserDlg.obj"
	-@erase ".\Debug\TimeDlg.obj"
	-@erase ".\Debug\Password.obj"
	-@erase ".\Debug\Prime95Doc.obj"
	-@erase ".\Debug\PreferencesDlg.obj"
	-@erase ".\Debug\EcmDlg.obj"
	-@erase ".\Debug\Prime95.obj"
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
# ADD LINK32 factor64.obj mult.obj mult1.obj mult2.obj mult2a.obj mult3.obj mult3a.obj mult4.obj mult4a.obj mult4b.obj mult1p.obj mult2p.obj mult2ap.obj mult3p.obj mult3ap.obj mult4p.obj mult4ap.obj mult4bp.obj mult4aux.obj ecmhelp.obj /nologo /subsystem:windows /debug /machine:I386
LINK32_FLAGS=factor64.obj mult.obj mult1.obj mult2.obj mult2a.obj mult3.obj\
 mult3a.obj mult4.obj mult4a.obj mult4b.obj mult1p.obj mult2p.obj mult2ap.obj\
 mult3p.obj mult3ap.obj mult4p.obj mult4ap.obj mult4bp.obj mult4aux.obj\
 ecmhelp.obj /nologo /subsystem:windows /incremental:yes\
 /pdb:"$(OUTDIR)/Prime95.pdb" /debug /machine:I386 /out:"$(OUTDIR)/Prime95.exe" 
LINK32_OBJS= \
	"$(INTDIR)/ManualCommDlg.obj" \
	"$(INTDIR)/AffinityDlg.obj" \
	"$(INTDIR)/Pminus1Dlg.obj" \
	"$(INTDIR)/StdAfx.obj" \
	"$(INTDIR)/VacationDlg.obj" \
	"$(INTDIR)/TestDlg.obj" \
	"$(INTDIR)/PrimenetDlg.obj" \
	"$(INTDIR)/Prime95View.obj" \
	"$(INTDIR)/CpuDlg.obj" \
	"$(INTDIR)/Priority.obj" \
	"$(INTDIR)/MainFrm.obj" \
	"$(INTDIR)/ServerDlg.obj" \
	"$(INTDIR)/UserDlg.obj" \
	"$(INTDIR)/TimeDlg.obj" \
	"$(INTDIR)/Password.obj" \
	"$(INTDIR)/Prime95Doc.obj" \
	"$(INTDIR)/PreferencesDlg.obj" \
	"$(INTDIR)/EcmDlg.obj" \
	"$(INTDIR)/Prime95.obj" \
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
	".\..\cpuid.h"\
	

"$(INTDIR)\Prime95.obj" : $(SOURCE) $(DEP_CPP_PRIME) "$(INTDIR)"\
 "$(INTDIR)\Prime95.pch"


!ELSEIF  "$(CFG)" == "Prime95 - Win32 Debug"

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
	".\..\commonb.h"\
	".\..\commonc.h"\
	".\..\comm95a.h"\
	".\..\comm95b.h"\
	".\..\comm95c.h"\
	".\..\primenet.h"\
	".\..\security.h"\
	

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
	".\..\cpuid.h"\
	

"$(INTDIR)\MainFrm.obj" : $(SOURCE) $(DEP_CPP_MAINF) "$(INTDIR)"\
 "$(INTDIR)\Prime95.pch"


!ELSEIF  "$(CFG)" == "Prime95 - Win32 Debug"

DEP_CPP_MAINF=\
	".\StdAfx.h"\
	".\Prime95.h"\
	".\MainFrm.h"\
	".\..\cpuid.h"\
	".\..\speed.h"\
	".\..\gwnum.h"\
	".\..\commona.h"\
	".\..\commonb.h"\
	".\..\commonc.h"\
	".\..\comm95a.h"\
	".\..\comm95b.h"\
	".\..\comm95c.h"\
	".\..\primenet.h"\
	".\..\security.h"\
	

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
	".\CpuDlg.h"\
	".\EcmDlg.h"\
	".\ManualCommDlg.h"\
	".\Password.h"\
	".\Pminus1Dlg.h"\
	".\PreferencesDlg.h"\
	".\PrimenetDlg.h"\
	".\Priority.h"\
	".\ServerDlg.h"\
	".\TestDlg.h"\
	".\TimeDlg.h"\
	".\UserDlg.h"\
	".\VacationDlg.h"\
	{$(INCLUDE)}"\sys\TIMEB.H"\
	".\..\cpuid.c"\
	".\..\speed.c"\
	".\..\giants.h"\
	".\..\giants.c"\
	".\..\gwnum.c"\
	".\..\commona.c"\
	".\..\commonb.c"\
	".\..\commonc.c"\
	".\..\ecm.c"\
	".\..\comm95a.c"\
	".\..\comm95b.c"\
	".\..\comm95c.c"\
	".\..\cpuid.h"\
	".\..\security.c"\
	

"$(INTDIR)\Prime95Doc.obj" : $(SOURCE) $(DEP_CPP_PRIME9) "$(INTDIR)"\
 "$(INTDIR)\Prime95.pch"


!ELSEIF  "$(CFG)" == "Prime95 - Win32 Debug"

DEP_CPP_PRIME9=\
	".\StdAfx.h"\
	".\MainFrm.h"\
	".\Prime95.h"\
	".\Prime95Doc.h"\
	".\AffinityDlg.h"\
	".\CpuDlg.h"\
	".\EcmDlg.h"\
	".\ManualCommDlg.h"\
	".\Password.h"\
	".\Pminus1Dlg.h"\
	".\PreferencesDlg.h"\
	".\PrimenetDlg.h"\
	".\Priority.h"\
	".\ServerDlg.h"\
	".\TestDlg.h"\
	".\TimeDlg.h"\
	".\UserDlg.h"\
	".\VacationDlg.h"\
	{$(INCLUDE)}"\sys\TIMEB.H"\
	".\..\cpuid.c"\
	".\..\speed.c"\
	".\..\giants.h"\
	".\..\giants.c"\
	".\..\gwnum.c"\
	".\..\commona.c"\
	".\..\commonb.c"\
	".\..\commonc.c"\
	".\..\ecm.c"\
	".\..\comm95a.c"\
	".\..\comm95b.c"\
	".\..\comm95c.c"\
	".\..\cpuid.h"\
	".\..\speed.h"\
	".\..\gwnum.h"\
	".\..\commona.h"\
	".\..\commonb.h"\
	".\..\commonc.h"\
	".\..\comm95a.h"\
	".\..\comm95b.h"\
	".\..\comm95c.h"\
	".\..\primenet.h"\
	".\..\security.h"\
	".\..\security.c"\
	

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
	".\..\cpuid.h"\
	

"$(INTDIR)\Prime95View.obj" : $(SOURCE) $(DEP_CPP_PRIME95) "$(INTDIR)"\
 "$(INTDIR)\Prime95.pch"


!ELSEIF  "$(CFG)" == "Prime95 - Win32 Debug"

DEP_CPP_PRIME95=\
	".\StdAfx.h"\
	".\Prime95.h"\
	".\Prime95Doc.h"\
	".\Prime95View.h"\
	".\..\cpuid.h"\
	".\..\speed.h"\
	".\..\gwnum.h"\
	".\..\commona.h"\
	".\..\commonb.h"\
	".\..\commonc.h"\
	".\..\comm95a.h"\
	".\..\comm95b.h"\
	".\..\comm95c.h"\
	".\..\primenet.h"\
	".\..\security.h"\
	

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
	".\..\cpuid.h"\
	

"$(INTDIR)\CpuDlg.obj" : $(SOURCE) $(DEP_CPP_CPUDL) "$(INTDIR)"\
 "$(INTDIR)\Prime95.pch"


!ELSEIF  "$(CFG)" == "Prime95 - Win32 Debug"

DEP_CPP_CPUDL=\
	".\StdAfx.h"\
	".\Prime95.h"\
	".\CpuDlg.h"\
	".\..\cpuid.h"\
	".\..\speed.h"\
	".\..\gwnum.h"\
	".\..\commona.h"\
	".\..\commonb.h"\
	".\..\commonc.h"\
	".\..\comm95a.h"\
	".\..\comm95b.h"\
	".\..\comm95c.h"\
	".\..\primenet.h"\
	".\..\security.h"\
	

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
	".\..\cpuid.h"\
	

"$(INTDIR)\PreferencesDlg.obj" : $(SOURCE) $(DEP_CPP_PREFE) "$(INTDIR)"\
 "$(INTDIR)\Prime95.pch"


!ELSEIF  "$(CFG)" == "Prime95 - Win32 Debug"

DEP_CPP_PREFE=\
	".\StdAfx.h"\
	".\Prime95.h"\
	".\PreferencesDlg.h"\
	".\..\cpuid.h"\
	".\..\speed.h"\
	".\..\gwnum.h"\
	".\..\commona.h"\
	".\..\commonb.h"\
	".\..\commonc.h"\
	".\..\comm95a.h"\
	".\..\comm95b.h"\
	".\..\comm95c.h"\
	".\..\primenet.h"\
	".\..\security.h"\
	

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
	".\..\cpuid.h"\
	

"$(INTDIR)\TestDlg.obj" : $(SOURCE) $(DEP_CPP_TESTD) "$(INTDIR)"\
 "$(INTDIR)\Prime95.pch"


!ELSEIF  "$(CFG)" == "Prime95 - Win32 Debug"

DEP_CPP_TESTD=\
	".\StdAfx.h"\
	".\Prime95.h"\
	".\TestDlg.h"\
	".\..\cpuid.h"\
	".\..\speed.h"\
	".\..\gwnum.h"\
	".\..\commona.h"\
	".\..\commonb.h"\
	".\..\commonc.h"\
	".\..\comm95a.h"\
	".\..\comm95b.h"\
	".\..\comm95c.h"\
	".\..\primenet.h"\
	".\..\security.h"\
	

"$(INTDIR)\TestDlg.obj" : $(SOURCE) $(DEP_CPP_TESTD) "$(INTDIR)"\
 "$(INTDIR)\Prime95.pch"


!ENDIF 

# End Source File
################################################################################
# Begin Source File

SOURCE=.\TimeDlg.cpp

!IF  "$(CFG)" == "Prime95 - Win32 Release"

DEP_CPP_TIMED=\
	".\StdAfx.h"\
	".\Prime95.h"\
	".\TimeDlg.h"\
	".\..\cpuid.h"\
	

"$(INTDIR)\TimeDlg.obj" : $(SOURCE) $(DEP_CPP_TIMED) "$(INTDIR)"\
 "$(INTDIR)\Prime95.pch"


!ELSEIF  "$(CFG)" == "Prime95 - Win32 Debug"

DEP_CPP_TIMED=\
	".\StdAfx.h"\
	".\Prime95.h"\
	".\TimeDlg.h"\
	".\..\cpuid.h"\
	".\..\speed.h"\
	".\..\gwnum.h"\
	".\..\commona.h"\
	".\..\commonb.h"\
	".\..\commonc.h"\
	".\..\comm95a.h"\
	".\..\comm95b.h"\
	".\..\comm95c.h"\
	".\..\primenet.h"\
	".\..\security.h"\
	

"$(INTDIR)\TimeDlg.obj" : $(SOURCE) $(DEP_CPP_TIMED) "$(INTDIR)"\
 "$(INTDIR)\Prime95.pch"


!ENDIF 

# End Source File
################################################################################
# Begin Source File

SOURCE=.\Password.cpp

!IF  "$(CFG)" == "Prime95 - Win32 Release"

DEP_CPP_PASSW=\
	".\StdAfx.h"\
	".\Prime95.h"\
	".\Password.h"\
	".\..\cpuid.h"\
	

"$(INTDIR)\Password.obj" : $(SOURCE) $(DEP_CPP_PASSW) "$(INTDIR)"\
 "$(INTDIR)\Prime95.pch"


!ELSEIF  "$(CFG)" == "Prime95 - Win32 Debug"

DEP_CPP_PASSW=\
	".\StdAfx.h"\
	".\Prime95.h"\
	".\Password.h"\
	".\..\cpuid.h"\
	".\..\speed.h"\
	".\..\gwnum.h"\
	".\..\commona.h"\
	".\..\commonb.h"\
	".\..\commonc.h"\
	".\..\comm95a.h"\
	".\..\comm95b.h"\
	".\..\comm95c.h"\
	".\..\primenet.h"\
	".\..\security.h"\
	

"$(INTDIR)\Password.obj" : $(SOURCE) $(DEP_CPP_PASSW) "$(INTDIR)"\
 "$(INTDIR)\Prime95.pch"


!ENDIF 

# End Source File
################################################################################
# Begin Source File

SOURCE=.\UserDlg.cpp

!IF  "$(CFG)" == "Prime95 - Win32 Release"

DEP_CPP_USERD=\
	".\StdAfx.h"\
	".\Prime95.h"\
	".\UserDlg.h"\
	".\..\cpuid.h"\
	

"$(INTDIR)\UserDlg.obj" : $(SOURCE) $(DEP_CPP_USERD) "$(INTDIR)"\
 "$(INTDIR)\Prime95.pch"


!ELSEIF  "$(CFG)" == "Prime95 - Win32 Debug"

DEP_CPP_USERD=\
	".\StdAfx.h"\
	".\Prime95.h"\
	".\UserDlg.h"\
	".\..\cpuid.h"\
	".\..\speed.h"\
	".\..\gwnum.h"\
	".\..\commona.h"\
	".\..\commonb.h"\
	".\..\commonc.h"\
	".\..\comm95a.h"\
	".\..\comm95b.h"\
	".\..\comm95c.h"\
	".\..\primenet.h"\
	".\..\security.h"\
	

"$(INTDIR)\UserDlg.obj" : $(SOURCE) $(DEP_CPP_USERD) "$(INTDIR)"\
 "$(INTDIR)\Prime95.pch"


!ENDIF 

# End Source File
################################################################################
# Begin Source File

SOURCE=.\PrimenetDlg.cpp

!IF  "$(CFG)" == "Prime95 - Win32 Release"

DEP_CPP_PRIMEN=\
	".\StdAfx.h"\
	".\Prime95.h"\
	".\PrimenetDlg.h"\
	".\..\cpuid.h"\
	

"$(INTDIR)\PrimenetDlg.obj" : $(SOURCE) $(DEP_CPP_PRIMEN) "$(INTDIR)"\
 "$(INTDIR)\Prime95.pch"


!ELSEIF  "$(CFG)" == "Prime95 - Win32 Debug"

DEP_CPP_PRIMEN=\
	".\StdAfx.h"\
	".\Prime95.h"\
	".\PrimenetDlg.h"\
	".\..\cpuid.h"\
	".\..\speed.h"\
	".\..\gwnum.h"\
	".\..\commona.h"\
	".\..\commonb.h"\
	".\..\commonc.h"\
	".\..\comm95a.h"\
	".\..\comm95b.h"\
	".\..\comm95c.h"\
	".\..\primenet.h"\
	".\..\security.h"\
	

"$(INTDIR)\PrimenetDlg.obj" : $(SOURCE) $(DEP_CPP_PRIMEN) "$(INTDIR)"\
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
	".\..\cpuid.h"\
	

"$(INTDIR)\Priority.obj" : $(SOURCE) $(DEP_CPP_PRIOR) "$(INTDIR)"\
 "$(INTDIR)\Prime95.pch"


!ELSEIF  "$(CFG)" == "Prime95 - Win32 Debug"

DEP_CPP_PRIOR=\
	".\StdAfx.h"\
	".\Prime95.h"\
	".\Priority.h"\
	".\..\cpuid.h"\
	".\..\speed.h"\
	".\..\gwnum.h"\
	".\..\commona.h"\
	".\..\commonb.h"\
	".\..\commonc.h"\
	".\..\comm95a.h"\
	".\..\comm95b.h"\
	".\..\comm95c.h"\
	".\..\primenet.h"\
	".\..\security.h"\
	

"$(INTDIR)\Priority.obj" : $(SOURCE) $(DEP_CPP_PRIOR) "$(INTDIR)"\
 "$(INTDIR)\Prime95.pch"


!ENDIF 

# End Source File
################################################################################
# Begin Source File

SOURCE=.\ServerDlg.cpp

!IF  "$(CFG)" == "Prime95 - Win32 Release"

DEP_CPP_SERVE=\
	".\StdAfx.h"\
	".\Prime95.h"\
	".\ServerDlg.h"\
	".\..\cpuid.h"\
	

"$(INTDIR)\ServerDlg.obj" : $(SOURCE) $(DEP_CPP_SERVE) "$(INTDIR)"\
 "$(INTDIR)\Prime95.pch"


!ELSEIF  "$(CFG)" == "Prime95 - Win32 Debug"

DEP_CPP_SERVE=\
	".\StdAfx.h"\
	".\Prime95.h"\
	".\ServerDlg.h"\
	".\..\cpuid.h"\
	".\..\speed.h"\
	".\..\gwnum.h"\
	".\..\commona.h"\
	".\..\commonb.h"\
	".\..\commonc.h"\
	".\..\comm95a.h"\
	".\..\comm95b.h"\
	".\..\comm95c.h"\
	".\..\primenet.h"\
	".\..\security.h"\
	

"$(INTDIR)\ServerDlg.obj" : $(SOURCE) $(DEP_CPP_SERVE) "$(INTDIR)"\
 "$(INTDIR)\Prime95.pch"


!ENDIF 

# End Source File
################################################################################
# Begin Source File

SOURCE=.\VacationDlg.cpp

!IF  "$(CFG)" == "Prime95 - Win32 Release"

DEP_CPP_VACAT=\
	".\StdAfx.h"\
	".\Prime95.h"\
	".\VacationDlg.h"\
	".\..\cpuid.h"\
	

"$(INTDIR)\VacationDlg.obj" : $(SOURCE) $(DEP_CPP_VACAT) "$(INTDIR)"\
 "$(INTDIR)\Prime95.pch"


!ELSEIF  "$(CFG)" == "Prime95 - Win32 Debug"

DEP_CPP_VACAT=\
	".\StdAfx.h"\
	".\Prime95.h"\
	".\VacationDlg.h"\
	".\..\cpuid.h"\
	".\..\speed.h"\
	".\..\gwnum.h"\
	".\..\commona.h"\
	".\..\commonb.h"\
	".\..\commonc.h"\
	".\..\comm95a.h"\
	".\..\comm95b.h"\
	".\..\comm95c.h"\
	".\..\primenet.h"\
	".\..\security.h"\
	

"$(INTDIR)\VacationDlg.obj" : $(SOURCE) $(DEP_CPP_VACAT) "$(INTDIR)"\
 "$(INTDIR)\Prime95.pch"


!ENDIF 

# End Source File
################################################################################
# Begin Source File

SOURCE=.\ManualCommDlg.cpp

!IF  "$(CFG)" == "Prime95 - Win32 Release"

DEP_CPP_MANUA=\
	".\StdAfx.h"\
	".\Prime95.h"\
	".\ManualCommDlg.h"\
	".\..\cpuid.h"\
	

"$(INTDIR)\ManualCommDlg.obj" : $(SOURCE) $(DEP_CPP_MANUA) "$(INTDIR)"\
 "$(INTDIR)\Prime95.pch"


!ELSEIF  "$(CFG)" == "Prime95 - Win32 Debug"

DEP_CPP_MANUA=\
	".\StdAfx.h"\
	".\Prime95.h"\
	".\ManualCommDlg.h"\
	".\..\cpuid.h"\
	".\..\speed.h"\
	".\..\gwnum.h"\
	".\..\commona.h"\
	".\..\commonb.h"\
	".\..\commonc.h"\
	".\..\comm95a.h"\
	".\..\comm95b.h"\
	".\..\comm95c.h"\
	".\..\primenet.h"\
	".\..\security.h"\
	

"$(INTDIR)\ManualCommDlg.obj" : $(SOURCE) $(DEP_CPP_MANUA) "$(INTDIR)"\
 "$(INTDIR)\Prime95.pch"


!ENDIF 

# End Source File
################################################################################
# Begin Source File

SOURCE=.\EcmDlg.cpp

!IF  "$(CFG)" == "Prime95 - Win32 Release"

DEP_CPP_ECMDL=\
	".\StdAfx.h"\
	".\Prime95.h"\
	".\EcmDlg.h"\
	".\..\cpuid.h"\
	

"$(INTDIR)\EcmDlg.obj" : $(SOURCE) $(DEP_CPP_ECMDL) "$(INTDIR)"\
 "$(INTDIR)\Prime95.pch"


!ELSEIF  "$(CFG)" == "Prime95 - Win32 Debug"

DEP_CPP_ECMDL=\
	".\StdAfx.h"\
	".\Prime95.h"\
	".\EcmDlg.h"\
	".\..\cpuid.h"\
	".\..\speed.h"\
	".\..\gwnum.h"\
	".\..\commona.h"\
	".\..\commonb.h"\
	".\..\commonc.h"\
	".\..\comm95a.h"\
	".\..\comm95b.h"\
	".\..\comm95c.h"\
	".\..\primenet.h"\
	".\..\security.h"\
	

"$(INTDIR)\EcmDlg.obj" : $(SOURCE) $(DEP_CPP_ECMDL) "$(INTDIR)"\
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
	".\..\cpuid.h"\
	

"$(INTDIR)\AffinityDlg.obj" : $(SOURCE) $(DEP_CPP_AFFIN) "$(INTDIR)"\
 "$(INTDIR)\Prime95.pch"


!ELSEIF  "$(CFG)" == "Prime95 - Win32 Debug"

DEP_CPP_AFFIN=\
	".\StdAfx.h"\
	".\Prime95.h"\
	".\AffinityDlg.h"\
	".\..\cpuid.h"\
	".\..\speed.h"\
	".\..\gwnum.h"\
	".\..\commona.h"\
	".\..\commonb.h"\
	".\..\commonc.h"\
	".\..\comm95a.h"\
	".\..\comm95b.h"\
	".\..\comm95c.h"\
	".\..\primenet.h"\
	".\..\security.h"\
	

"$(INTDIR)\AffinityDlg.obj" : $(SOURCE) $(DEP_CPP_AFFIN) "$(INTDIR)"\
 "$(INTDIR)\Prime95.pch"


!ENDIF 

# End Source File
################################################################################
# Begin Source File

SOURCE=.\Pminus1Dlg.cpp

!IF  "$(CFG)" == "Prime95 - Win32 Release"

DEP_CPP_PMINU=\
	".\StdAfx.h"\
	".\Prime95.h"\
	".\Pminus1Dlg.h"\
	".\..\cpuid.h"\
	

"$(INTDIR)\Pminus1Dlg.obj" : $(SOURCE) $(DEP_CPP_PMINU) "$(INTDIR)"\
 "$(INTDIR)\Prime95.pch"


!ELSEIF  "$(CFG)" == "Prime95 - Win32 Debug"

DEP_CPP_PMINU=\
	".\StdAfx.h"\
	".\Prime95.h"\
	".\Pminus1Dlg.h"\
	".\..\cpuid.h"\
	".\..\speed.h"\
	".\..\gwnum.h"\
	".\..\commona.h"\
	".\..\commonb.h"\
	".\..\commonc.h"\
	".\..\comm95a.h"\
	".\..\comm95b.h"\
	".\..\comm95c.h"\
	".\..\primenet.h"\
	".\..\security.h"\
	

"$(INTDIR)\Pminus1Dlg.obj" : $(SOURCE) $(DEP_CPP_PMINU) "$(INTDIR)"\
 "$(INTDIR)\Prime95.pch"


!ENDIF 

# End Source File
# End Target
# End Project
################################################################################
