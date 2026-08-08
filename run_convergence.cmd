@echo off
setlocal
set "SCRIPT_DIR=%~dp0"
set "SOURCE=%SCRIPT_DIR%convergence_analyzer.cpp"
set "BINARY=%SCRIPT_DIR%convergence_analyzer.exe"

if exist "%BINARY%" goto run

g++ -O3 -std=c++17 -DNDEBUG "%SOURCE%" -o "%BINARY%"
if errorlevel 1 (
  echo error: compilation failed. Ensure a C++17 g++ compiler is available. 1>&2
  exit /b 2
)

:run
"%BINARY%" %*
exit /b %errorlevel%
