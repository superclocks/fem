# Microsoft Developer Studio Project File - Name="check one element" - Package Owner=<4>
# Microsoft Developer Studio Generated Build File, Format Version 6.00
# ** DO NOT EDIT **

# TARGTYPE "Win32 (x86) Console Application" 0x0103

CFG=check one element - Win32 Debug
!MESSAGE This is not a valid makefile. To build this project using NMAKE,
!MESSAGE use the Export Makefile command and run
!MESSAGE 
!MESSAGE NMAKE /f "check one element.mak".
!MESSAGE 
!MESSAGE You can specify a configuration when running NMAKE
!MESSAGE by defining the macro CFG on the command line. For example:
!MESSAGE 
!MESSAGE NMAKE /f "check one element.mak" CFG="check one element - Win32 Debug"
!MESSAGE 
!MESSAGE Possible choices for configuration are:
!MESSAGE 
!MESSAGE "check one element - Win32 Release" (based on "Win32 (x86) Console Application")
!MESSAGE "check one element - Win32 Debug" (based on "Win32 (x86) Console Application")
!MESSAGE 

# Begin Project
# PROP AllowPerConfigDependencies 0
# PROP Scc_ProjName ""
# PROP Scc_LocalPath ""
CPP=cl.exe
F90=df.exe
RSC=rc.exe

!IF  "$(CFG)" == "check one element - Win32 Release"

# PROP BASE Use_MFC 0
# PROP BASE Use_Debug_Libraries 0
# PROP BASE Output_Dir "Release"
# PROP BASE Intermediate_Dir "Release"
# PROP BASE Target_Dir ""
# PROP Use_MFC 0
# PROP Use_Debug_Libraries 0
# PROP Output_Dir "Release"
# PROP Intermediate_Dir "Release"
# PROP Target_Dir ""
# ADD BASE F90 /compile_only /nologo /warn:nofileopt
# ADD F90 /compile_only /nologo /warn:nofileopt
# ADD BASE CPP /nologo /W3 /GX /O2 /D "WIN32" /D "NDEBUG" /D "_CONSOLE" /D "_MBCS" /YX /FD /c
# ADD CPP /nologo /W3 /GX /O2 /D "WIN32" /D "NDEBUG" /D "_CONSOLE" /D "_MBCS" /YX /FD /c
# ADD BASE RSC /l 0x804 /d "NDEBUG"
# ADD RSC /l 0x804 /d "NDEBUG"
BSC32=bscmake.exe
# ADD BASE BSC32 /nologo
# ADD BSC32 /nologo
LINK32=link.exe
# ADD BASE LINK32 kernel32.lib user32.lib gdi32.lib winspool.lib comdlg32.lib advapi32.lib shell32.lib ole32.lib oleaut32.lib uuid.lib odbc32.lib odbccp32.lib /nologo /subsystem:console /machine:I386
# ADD LINK32 kernel32.lib user32.lib gdi32.lib winspool.lib comdlg32.lib advapi32.lib shell32.lib ole32.lib oleaut32.lib uuid.lib odbc32.lib odbccp32.lib /nologo /subsystem:console /machine:I386

!ELSEIF  "$(CFG)" == "check one element - Win32 Debug"

# PROP BASE Use_MFC 0
# PROP BASE Use_Debug_Libraries 1
# PROP BASE Output_Dir "Debug"
# PROP BASE Intermediate_Dir "Debug"
# PROP BASE Target_Dir ""
# PROP Use_MFC 0
# PROP Use_Debug_Libraries 1
# PROP Output_Dir "Debug"
# PROP Intermediate_Dir "Debug"
# PROP Ignore_Export_Lib 0
# PROP Target_Dir ""
# ADD BASE F90 /check:bounds /compile_only /debug:full /nologo /traceback /warn:argument_checking /warn:nofileopt
# ADD F90 /check:bounds /compile_only /debug:full /nologo /traceback /warn:argument_checking /warn:nofileopt
# ADD BASE CPP /nologo /W3 /Gm /GX /ZI /Od /D "WIN32" /D "_DEBUG" /D "_CONSOLE" /D "_MBCS" /YX /FD /GZ /c
# ADD CPP /nologo /W3 /Gm /GX /ZI /Od /D "WIN32" /D "_DEBUG" /D "_CONSOLE" /D "_MBCS" /FR /YX /FD /GZ /c
# ADD BASE RSC /l 0x804 /d "_DEBUG"
# ADD RSC /l 0x804 /d "_DEBUG"
BSC32=bscmake.exe
# ADD BASE BSC32 /nologo
# ADD BSC32 /nologo
LINK32=link.exe
# ADD BASE LINK32 kernel32.lib user32.lib gdi32.lib winspool.lib comdlg32.lib advapi32.lib shell32.lib ole32.lib oleaut32.lib uuid.lib odbc32.lib odbccp32.lib /nologo /subsystem:console /debug /machine:I386 /pdbtype:sept
# ADD LINK32 kernel32.lib user32.lib gdi32.lib winspool.lib comdlg32.lib advapi32.lib shell32.lib ole32.lib oleaut32.lib uuid.lib odbc32.lib odbccp32.lib GSS_DLL_6.lib /nologo /stack:0xffffffff,0xffffffff /subsystem:console /incremental:no /debug /machine:I386 /pdbtype:sept

!ENDIF 

# Begin Target

# Name "check one element - Win32 Release"
# Name "check one element - Win32 Debug"
# Begin Group "Source Files"

# PROP Default_Filter "cpp;c;cxx;rc;def;r;odl;idl;hpj;bat;f90;for;f;fpp"
# Begin Source File

SOURCE=".\featp\ASSEM(3).FOR"
# End Source File
# Begin Source File

SOURCE=.\featp\ASSEM_FLUID30.FOR
# End Source File
# Begin Source File

SOURCE=.\featp\ASSEM_SOLID45.FOR
# End Source File
# Begin Source File

SOURCE=.\featp\BEND_MASS_TRI_I.FOR
# End Source File
# Begin Source File

SOURCE=.\featp\BEND_MASS_TRI_II.FOR
# End Source File
# Begin Source File

SOURCE=.\featp\BEND_VD.FOR
# End Source File
# Begin Source File

SOURCE=.\featp\CACUL_OUT_NORMAL_QUAD.FOR
# End Source File
# Begin Source File

SOURCE=.\featp\CACUL_OUT_NORMAL_TRI.FOR
# End Source File
# Begin Source File

SOURCE=.\featp\CHECK_SHAPE.FOR
# End Source File
# Begin Source File

SOURCE=.\featp\DKT_I.FOR
# End Source File
# Begin Source File

SOURCE=.\featp\DKT_II.FOR
# End Source File
# Begin Source File

SOURCE=".\featp\ELEMENT_JACOBI(3_1_7).FOR"
# End Source File
# Begin Source File

SOURCE=.\featp\ELEMENT_MATRIX_QUAD.FOR
# End Source File
# Begin Source File

SOURCE=.\featp\ELEMENT_MATRIX_SOLID.FOR
# End Source File
# Begin Source File

SOURCE=.\featp\ELEMENT_MATRIX_TRI.FOR
# End Source File
# Begin Source File

SOURCE=".\featp\ELEMENT_VB(3_1_8).FOR"
# End Source File
# Begin Source File

SOURCE=".\featp\ELEMENT_VD(3_1_6).FOR"
# End Source File
# Begin Source File

SOURCE=.\featp\FEATP.FOR
# End Source File
# Begin Source File

SOURCE=.\featp\FLUID30_6_I.FOR
DEP_F90_FLUID=\
	{$(INCLUDE)}"IMSL.mod"\
	
# End Source File
# Begin Source File

SOURCE=.\featp\FLUID30_6_II.FOR
# End Source File
# Begin Source File

SOURCE=.\featp\FLUID30_8_I.FOR
DEP_F90_FLUID3=\
	{$(INCLUDE)}"IMSL.mod"\
	
# End Source File
# Begin Source File

SOURCE=.\featp\FLUID30_8_II.FOR
# End Source File
# Begin Source File

SOURCE=".\featp\GAUSS_INTEGRATION(3_1_1).FOR"
# End Source File
# Begin Source File

SOURCE=.\solver\gss\gss_6_demo.for
# End Source File
# Begin Source File

SOURCE=.\featp\MEMBRANE_I.FOR
DEP_F90_MEMBR=\
	{$(INCLUDE)}"IMSL.mod"\
	
# End Source File
# Begin Source File

SOURCE=.\featp\MEMBRANE_II.FOR
# End Source File
# Begin Source File

SOURCE=.\featp\MEMBRANE_MASS.FOR
# End Source File
# Begin Source File

SOURCE=".\featp\MEMBRANE_STIFF_QUAD(0)_I.FOR"
# End Source File
# Begin Source File

SOURCE=".\featp\MEMBRANE_STIFF_QUAD(0)_II.FOR"
DEP_F90_MEMBRA=\
	{$(INCLUDE)}"IMSL.mod"\
	
# End Source File
# Begin Source File

SOURCE=".\featp\MEMBRANE_STIFF_QUAD(1)_I.FOR"
# End Source File
# Begin Source File

SOURCE=".\featp\MEMBRANE_STIFF_QUAD(1)_II.FOR"
# End Source File
# Begin Source File

SOURCE=".\featp\Allman\MEMBRANE_STIFF_QUAD(2)_I.FOR"
# End Source File
# Begin Source File

SOURCE=".\featp\Allman\MEMBRANE_STIFF_QUAD(2)_II.FOR"
DEP_F90_MEMBRAN=\
	{$(INCLUDE)}"IMSL.mod"\
	
# End Source File
# Begin Source File

SOURCE=".\featp\MEMBRANE_STIFF_TRI(0&1)_I.FOR"
# End Source File
# Begin Source File

SOURCE=".\featp\MEMBRANE_STIFF_TRI(0&1)_II.FOR"
# End Source File
# Begin Source File

SOURCE=".\featp\MEMBRANE_STIFF_TRI(2)_I.FOR"
# End Source File
# Begin Source File

SOURCE=".\featp\MEMBRANE_STIFF_TRI(2)_II.FOR"
# End Source File
# Begin Source File

SOURCE=.\featp\MEMBRANE_VD.FOR
# End Source File
# Begin Source File

SOURCE=.\featp\QUAD_COUP_I.FOR
# End Source File
# Begin Source File

SOURCE=.\featp\QUAD_COUP_II.FOR
# End Source File
# Begin Source File

SOURCE=.\featp\QUAD_TRANS.FOR
# End Source File
# Begin Source File

SOURCE=.\featp\REORDER_SOLID_NODE_NUM.FOR
# End Source File
# Begin Source File

SOURCE=.\featp\SHAPE_QUADRANGLE_4.FOR
# End Source File
# Begin Source File

SOURCE=".\featp\SOLID(0)_I.FOR"
DEP_F90_SOLID=\
	{$(INCLUDE)}"IMSL.mod"\
	
# End Source File
# Begin Source File

SOURCE=".\featp\SOLID(0)_II.FOR"
DEP_F90_SOLID_=\
	{$(INCLUDE)}"IMSL.mod"\
	
# End Source File
# Begin Source File

SOURCE=".\featp\SOLID(1)_I.FOR"
DEP_F90_SOLID_1=\
	{$(INCLUDE)}"IMSL.mod"\
	
# End Source File
# Begin Source File

SOURCE=".\featp\SOLID(1)_II.FOR"
# End Source File
# Begin Source File

SOURCE=.\featp\SOLID_TRANS.FOR
# End Source File
# Begin Source File

SOURCE=.\solver\gss\solver.for
# End Source File
# Begin Source File

SOURCE=.\featp\TRI_COUP_I.FOR
# End Source File
# Begin Source File

SOURCE=.\featp\TRI_COUP_II.FOR
# End Source File
# Begin Source File

SOURCE=.\featp\TRI_TRANS.FOR
# End Source File
# End Group
# Begin Group "Header Files"

# PROP Default_Filter "h;hpp;hxx;hm;inl;fi;fd"
# End Group
# Begin Group "Resource Files"

# PROP Default_Filter "ico;cur;bmp;dlg;rc2;rct;bin;rgs;gif;jpg;jpeg;jpe"
# End Group
# End Target
# End Project
