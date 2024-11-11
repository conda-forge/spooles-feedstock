make lib CC=%CC% AR=%AR%
if errorlevel 1 exit 1

cd MT\src\
make CC=%CC% AR=%AR%
if errorlevel 1 exit 1

mkdir -p %LIBRARY_LIB%
cp *.a %LIBRARY_LIB%
if errorlevel 1 exit 1

cd ..\..
cp *.a %LIBRARY_LIB%
if errorlevel 1 exit 1

mkdir %LIBRARY_INC%\spooles
if errorlevel 1 exit 1

xcopy /s .\*.h %LIBRARY_INC%\spooles
