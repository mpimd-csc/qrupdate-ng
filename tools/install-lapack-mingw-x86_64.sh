#!/bin/bash

#!/bin/bash
set -x

[ -d build-lapack-win64 ] && rm -rf build-lapack-win64

[ -e v3.10.0.tar.gz ] || wget -c https://github.com/Reference-LAPACK/lapack/archive/refs/tags/v3.10.0.tar.gz
tar xf v3.10.0.tar.gz
cmake -S lapack-3.10.0 -B build-lapack-win64 \
    -DCMAKE_TOOLCHAIN_FILE=$(pwd)/../cmake/mingw-w64-x86_64.cmake \
    -DBUILD_SHARED_LIBS=ON \
    -DCMAKE_INSTALL_PREFIX=`dirname $(realpath $0)`/x86_64-w64-mingw32
make -C build-lapack-win64 -j 4
make -C build-lapack-win64 install
