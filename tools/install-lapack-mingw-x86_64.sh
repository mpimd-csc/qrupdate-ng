#!/bin/bash
export LAPACK_VERSION=3.12.1
set -x

[ -d build-lapack-win64 ] && rm -rf build-lapack-win64

[ -e v${LAPACK_VERSION}.tar.gz ] || wget -c https://github.com/Reference-LAPACK/lapack/archive/refs/tags/v${LAPACK_VERSION}.tar.gz
tar xf v${LAPACK_VERSION}.tar.gz
cmake -S lapack-${LAPACK_VERSION} -B build-lapack-win64 \
    -DCMAKE_TOOLCHAIN_FILE=$(pwd)/../cmake/mingw-w64-x86_64.cmake \
    -DBUILD_SHARED_LIBS=ON \
    -DCMAKE_INSTALL_PREFIX=`dirname $(realpath $0)`/x86_64-w64-mingw32
make -C build-lapack-win64 -j 4
make -C build-lapack-win64 install
