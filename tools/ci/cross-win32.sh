#!/bin/bash
set -x
(cd tools/; bash install-lapack-mingw-i686.sh)
cmake -S . -B build.win32 -DCMAKE_TOOLCHAIN_FILE=$(pwd)/cmake/mingw-w32-i686.cmake -DBUILD_SHARED_LIBS=ON
make -C build.win32 all
(cd build.win32; ctest -V )

