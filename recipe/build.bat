@echo on

set "BOOST_ROOT=%PREFIX%\Library"
set "BOOSTROOT=%PREFIX%\Library"

meson setup build_win ^
    --prefix="%PREFIX%" ^
    --libdir="%PREFIX%\Library\lib" ^
    --includedir="%PREFIX%\Library\include" ^
    --bindir="%PREFIX%\Library\bin" ^
    --buildtype=release ^
    -Dbuild_cli=enabled ^
    -Dbuild_library=enabled ^
    -Dbuild_doc=false ^
    -Dcpp_args="-I%PREFIX%\Library\include" ^
    -Dcpp_link_args="-L%PREFIX%\Library\lib"
if errorlevel 1 exit 1

meson compile -C build_win -v -j 2
if errorlevel 1 exit 1

meson install -C build_win
if errorlevel 1 exit 1
