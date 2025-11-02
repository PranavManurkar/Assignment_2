@echo off
setlocal enabledelayedexpansion

set LAMMPS_EXE=lmp
set MAX_PARALLEL=9

echo ======================================================================
echo Al-Cr-Co Parallel Simulation Runner
echo Running %MAX_PARALLEL% simulations simultaneously per composition
echo ======================================================================
echo.

REM Check if LAMMPS exists
where %LAMMPS_EXE% >nul 2>&1
if errorlevel 1 (
    echo ERROR: LAMMPS executable '%LAMMPS_EXE%' not found!
    echo Please install LAMMPS or update LAMMPS_EXE in this script
    pause
    exit /b 1
)

set TOTAL_COMP=0
set TOTAL_SIMS=0

REM Loop through each composition folder
for /d %%C in (Comp*) do (
    set /a TOTAL_COMP+=1
    echo.
    echo ======================================================================
    echo [!TOTAL_COMP!] Processing: %%C
    echo ======================================================================
    
    cd %%C
    
    REM Check for required files
    if not exist "library.meam" (
        echo   ERROR: Missing library.meam in %%C
        cd ..
        goto :next_comp
    )
    if not exist "CrCoAl.meam" (
        echo   ERROR: Missing CrCoAl.meam in %%C
        cd ..
        goto :next_comp
    )
    
    echo Starting 9 parallel simulations...
    echo.
    
    REM Start all 9 simulations in parallel
    start /B "" %LAMMPS_EXE% -in in.FCC_100K.lammps -log log.FCC.100K.lammps
    echo   [1/9] Started: FCC_100K
    
    start /B "" %LAMMPS_EXE% -in in.FCC_350K.lammps -log log.FCC.350K.lammps
    echo   [2/9] Started: FCC_350K
    
    start /B "" %LAMMPS_EXE% -in in.FCC_550K.lammps -log log.FCC.550K.lammps
    echo   [3/9] Started: FCC_550K
    
    start /B "" %LAMMPS_EXE% -in in.HCP_100K.lammps -log log.HCP.100K.lammps
    echo   [4/9] Started: HCP_100K
    
    start /B "" %LAMMPS_EXE% -in in.HCP_350K.lammps -log log.HCP.350K.lammps
    echo   [5/9] Started: HCP_350K
    
    start /B "" %LAMMPS_EXE% -in in.HCP_550K.lammps -log log.HCP.550K.lammps
    echo   [6/9] Started: HCP_550K
    
    start /B "" %LAMMPS_EXE% -in in.DHCP_100K.lammps -log log.DHCP.100K.lammps
    echo   [7/9] Started: DHCP_100K
    
    start /B "" %LAMMPS_EXE% -in in.DHCP_350K.lammps -log log.DHCP.350K.lammps
    echo   [8/9] Started: DHCP_350K
    
    start /B "" %LAMMPS_EXE% -in in.DHCP_550K.lammps -log log.DHCP.550K.lammps
    echo   [9/9] Started: DHCP_550K
    
    echo.
    echo Waiting for all 9 simulations to complete...
    
    REM Wait for all LAMMPS processes to finish
    :wait_loop
    tasklist /FI "IMAGENAME eq lmp.exe" 2>NUL | find /I /N "lmp.exe">NUL
    if "%ERRORLEVEL%"=="0" (
        timeout /t 10 /nobreak >nul
        goto wait_loop
    )
    
    echo All simulations in %%C completed!
    set /a TOTAL_SIMS+=9
    
    cd ..
    
    :next_comp
)

echo.
echo ======================================================================
echo All Compositions Completed!
echo ======================================================================
echo Total compositions processed: !TOTAL_COMP!
echo Total simulations run: !TOTAL_SIMS!
echo ======================================================================
pause