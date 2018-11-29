CP %BUILD_PREFIX%\Library\mingw-w64\bin\mingw32-make %BUILD_PREFIX%\Library\mingw-w64\bin\make
mingw32-make lib
if errorlevel 1 exit 1

cd MT\src\
mingw32-make
if errorlevel 1 exit 1

mkdir -p %LIBRARY_PREFIX%\lib
cp *.a %LIBRARY_PREFIX%\lib\
if errorlevel 1 exit 1

cd ..\..
cp *.a %LIBRARY_PREFIX%\lib\
if errorlevel 1 exit 1

mkdir %LIBRARY_PREFIX%\include\spooles
if errorlevel 1 exit 1

find . -type f -name '*.h' -exec cp -p --parents {} %LIBRARY_PREFIX%\include\spooles ";"
