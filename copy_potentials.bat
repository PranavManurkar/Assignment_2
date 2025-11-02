@echo off
echo Copying MEAM potential files to all composition folders...
echo.

for /d %%D in (Comp*) do (
    echo Copying to %%D
    copy /Y library.meam "%%D\" >nul
    copy /Y CrCoAl.meam "%%D\" >nul
)

echo.
echo Done! All potential files copied.
pause