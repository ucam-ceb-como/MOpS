# Microsoft Developer Studio Project File - Name="limex" - Package Owner=<4>
# Microsoft Developer Studio Generated Build File, Format Version 6.00
# ** DO NOT EDIT **

# TARGTYPE "Win32 (x86) Console Application" 0x0103

CFG=limex - Win32 Debug
!MESSAGE This is not a valid makefile. To build this project using NMAKE,
!MESSAGE use the Export Makefile command and run
!MESSAGE 
!MESSAGE NMAKE /f "limex.mak".
!MESSAGE 
!MESSAGE You can specify a configuration when running NMAKE
!MESSAGE by defining the macro CFG on the command line. For example:
!MESSAGE 
!MESSAGE NMAKE /f "limex.mak" CFG="limex - Win32 Debug"
!MESSAGE 
!MESSAGE Possible choices for configuration are:
!MESSAGE 
!MESSAGE "limex - Win32 Release" (based on "Win32 (x86) Console Application")
!MESSAGE "limex - Win32 Debug" (based on "Win32 (x86) Console Application")
!MESSAGE 

# Begin Project
# PROP AllowPerConfigDependencies 0
# PROP Scc_ProjName ""
# PROP Scc_LocalPath ""
CPP=cl.exe
F90=df.exe
RSC=rc.exe

!IF  "$(CFG)" == "limex - Win32 Release"

# PROP BASE Use_MFC 0
# PROP BASE Use_Debug_Libraries 0
# PROP BASE Output_Dir "Release"
# PROP BASE Intermediate_Dir "Release"
# PROP BASE Target_Dir ""
# PROP Use_MFC 0
# PROP Use_Debug_Libraries 0
# PROP Output_Dir "Release"
# PROP Intermediate_Dir "Release"
# PROP Ignore_Export_Lib 0
# PROP Target_Dir ""
# ADD BASE F90 /compile_only /nologo /warn:nofileopt
# ADD F90 /compile_only /nologo /warn:nofileopt
# ADD BASE CPP /nologo /W3 /GX /O2 /D "WIN32" /D "NDEBUG" /D "_CONSOLE" /D "_MBCS" /YX /FD /c
# ADD CPP /nologo /MT /W3 /GX /O2 /D "WIN32" /D "NDEBUG" /D "_CONSOLE" /D "_MBCS" /YX /FD /c
# ADD BASE RSC /l 0x809 /d "NDEBUG"
# ADD RSC /l 0x809 /d "NDEBUG"
BSC32=bscmake.exe
# ADD BASE BSC32 /nologo
# ADD BSC32 /nologo
LINK32=link.exe
# ADD BASE LINK32 kernel32.lib /nologo /subsystem:console /machine:I386
# ADD LINK32 kernel32.lib, user32.lib, gdi32.lib, winspool.lib, comdgl32.lib, advapi32.lib, shell32.lib, ole32.lib, oleaut32.lib, uuid.lib, odbc32.lib, odbccp32.lib, camxml.lib, camstrings.lib, sprog.lib, cvode.lib, sweep.lib, brushlib.lib /nologo /subsystem:console /machine:I386 /libpath:"C:\Documents and Settings\ecw46\Desktop\Camflow With Limex\camxml\lib" /libpath:"C:\Documents and Settings\ecw46\Desktop\Camflow With Limex\comostrings\lib" /libpath:"C:\Documents and Settings\ecw46\Desktop\Camflow With Limex\sprogc\lib" /libpath:"C:\Documents and Settings\ecw46\Desktop\Camflow With Limex\cvode\lib" /libpath:"C:\Documents and Settings\ecw46\Desktop\Camflow With Limex\sweepc\lib" /libpath:"C:\Documents and Settings\ecw46\Desktop\Camflow With Limex\brush\Release" /libpath:"C:\Program Files\Microsoft Visual Studio\VC98\LIB"

!ELSEIF  "$(CFG)" == "limex - Win32 Debug"

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
# ADD BASE F90 /check:bounds /compile_only /dbglibs /debug:full /nologo /traceback /warn:argument_checking /warn:nofileopt
# ADD F90 /check:bounds /compile_only /dbglibs /debug:full /nologo /traceback /warn:argument_checking /warn:nofileopt
# ADD BASE CPP /nologo /W3 /Gm /GX /ZI /Od /D "WIN32" /D "_DEBUG" /D "_CONSOLE" /D "_MBCS" /YX /FD /GZ /c
# ADD CPP /nologo /W3 /Gm /GX /ZI /Od /D "WIN32" /D "_DEBUG" /D "_CONSOLE" /D "_MBCS" /YX /FD /GZ /c
# ADD BASE RSC /l 0x809 /d "_DEBUG"
# ADD RSC /l 0x809 /d "_DEBUG"
BSC32=bscmake.exe
# ADD BASE BSC32 /nologo
# ADD BSC32 /nologo
LINK32=link.exe
# ADD BASE LINK32 kernel32.lib /nologo /subsystem:console /debug /machine:I386 /pdbtype:sept
# ADD LINK32 kernel32.lib /nologo /subsystem:console /profile /debug /machine:I386

!ENDIF 

# Begin Target

# Name "limex - Win32 Release"
# Name "limex - Win32 Debug"
# Begin Source File

SOURCE=.\limex.f
DEP_F90_LIMEX=\
	".\xlimex_def.h"\
	
# End Source File
# Begin Source File

SOURCE=.\LIMEX4_2_Auxiliaries.f
# End Source File
# Begin Source File

SOURCE=.\LIMEX4_2_Dense.f
# End Source File
# End Target
# End Project
