@echo off

setlocal
set LIBPATH=galib/lib
set INCPATH=galib/include
set OPT=-O2 -Wall -Wno-unused-dummy-argument -Wno-uninitialized -Wno-conversion -static -ffree-line-length-0 -std=gnu
if "%1"=="" (
  set EXE_FILE=a.exe
) else (
  set EXE_FILE=%1.exe
)
echo %OS%
if "%OS%"=="Windows_NT" (
  del %EXE_FILE% > NUL 2>&1
) else (
  rm -f %EXE_FILE%
)
call :main
endlocal
exit /b

:compile
setlocal
set cmd=gfortran -c %OPT% -o %~1.o %~1.f95 -I%INCPATH%
call :exec "%cmd%"
endlocal
exit /b

:link
setlocal
set cmd=gfortran %OPT% -o %* -L%LIBPATH% -lfec -llapack -lrefblas
call :exec "%cmd%"
endlocal
exit /b

:exec
echo %~1
%~1
exit /b

:main
:: ↓ファイルが1つの場合
call :compile dofdof
call :link %EXE_FILE% dofdof.o

:: ↓ファイルが複数の場合
:: call :compile dof_kriging
:: call :compile dof_rk4
:: call :compile dof_flight
:: call :compile dof_main
:: call :link %EXE_FILE% dof_kriging.o dof_rk4.o dof_flight.o dof_main.o

del *.o *.mod *.smod > NUL 2>&1
exit /b
