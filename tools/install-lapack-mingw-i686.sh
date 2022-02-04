#!/bin/bash
set -x

[ -d build-lapack-win32 ] && rm -rf build-lapack-win32

[ -e v3.10.0.tar.gz ] || wget -c https://github.com/Reference-LAPACK/lapack/archive/refs/tags/v3.10.0.tar.gz
tar xf v3.10.0.tar.gz
cmake -S lapack-3.10.0 -B build-lapack-win32 \
    -DCMAKE_TOOLCHAIN_FILE=$(pwd)/../cmake/mingw-w32-i686.cmake \
    -DBUILD_SHARED_LIBS=ON \
    -DCMAKE_INSTALL_PREFIX=`dirname $(realpath $0)`/i686-w64-mingw32
make -C build-lapack-win32 -j 4
make -C build-lapack-win32 install
