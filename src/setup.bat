@echo off

setlocal
set LIBPATH=galib\lib
set INCPATH=galib\include
set OPT=-O2 -Wall -Wno-unused-dummy-argument -Wno-unused-function -Wno-uninitialized -Wno-conversion -static -ffree-line-length-0 -std=gnu -fbacktrace -fbounds-check -fopenmp
rem set OPT=-O2 -Wall -static -ffree-line-length-0 -std=gnu -fbacktrace -fbounds-check

if "%1"=="" (
  set EXE_FILE=a.exe
) else (
  set EXE_FILE=%1.exe
)
echo %OS%
if "%OS%"=="Windows_NT" (
  del %EXE_FILE% >NUL 2>&1
) else (
  rm -f %EXE_FILE%
)
rem set OMP_NUM_THREADS=20
rem pushd src
call :main
rem popd
endlocal
exit /b


:compile
setlocal
set cmd=gfortran -c %OPT% -o %~1.o %~1.f08 -I%INCPATH%
call :exec "%cmd%"
endlocal
exit /b


:link
setlocal
set cmd=gfortran %OPT% -o %* -L%LIBPATH% -lfec -llapack -lrefblas
call :exec "%cmd%"
endlocal
exit /b


:mkdll
setlocal
set cmd=gfortran %OPT% -shared -o %* -L%LIBPATH% -lfec -llapack -lrefblas
call :exec "%cmd%"
endlocal
exit /b


:exec
echo %~1
%~1
exit /b


:main
rem call :compile util
rem call :compile matrix
rem call :compile kriging
call :compile dof_kriging
call :compile dof_base
rem call :compile dofdof
call :compile dof_new
rem call :compile dof_old

rem [A] 動的リンクライブラリ作成↓
rem call :mkdll libdof.dll dof_kriging.o dofdof.o

rem [B] 実行ファイル作成↓
call :link %EXE_FILE% dof_kriging.o dof_base.o dof_new.o
rem call :link %EXE_FILE% dof_old.o

copy a.exe ..\dof_new.exe
rem copy a.exe ..\dof_old.exe

rem del *.o *.mod *.smod >NUL 2>&1
exit /b
