#!/bin/bash
set -x
(cd tools/; bash install-lapack-mingw-x86_64.sh)
cmake -S . -B build.win64 -DCMAKE_TOOLCHAIN_FILE=$(pwd)/cmake/mingw-w64-x86_64.cmake -DBUILD_SHARED_LIBS=ON
make -C build.win64 all
(cd build.win64; ctest -V )

