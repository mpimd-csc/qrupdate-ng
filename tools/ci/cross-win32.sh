#!/bin/bash
export WINEARCH=win32
export WINEPREFIX=$(mktemp -d -p "$(pwd)")
set -x
 (cd tools/; bash install-lapack-mingw-i686.sh)
cmake -S . -B build.win32 -DCMAKE_TOOLCHAIN_FILE=$(pwd)/cmake/mingw-w32-i686.cmake -DBUILD_SHARED_LIBS=ON --fresh
make -C build.win32 all
(cd build.win32; ctest -V )
rm -rf "${WINEPREFIX}"

