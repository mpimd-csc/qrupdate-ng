Installation Notes
==================

The library is written in Fortran. It is Fortran-90/95-compliant. Fortran 77
compilers are no longer supported. This includes the f2c translator as well.

## Requirements

The library requires a Fortran 90/95 compiler, BLAS, and LAPACK. The following
compilers and BLAS/LAPACK implementations are tested under Linux. Other
operating systems, which are supported by CMAKE, should work as well.

Tested Fortran compilers:

 * GNU gfortran 4.8.5, 5.4, 7.3, 8.3, 9.2, 9.3
 * Intel ifort 18.0.1
 * IBM XLF 16.1.1

Tested BLAS/LAPACK implementations:

 * Reference BLAS/LAPACK, http://www.netlib.org/lapack
 * OpenBLAS, https://www.openblas.net
 * FlexiBLAS, https://www.mpi-magdeburg.mpg.de/projects/ flexiblas
 * Intel MKL, https://software.intel.com/en-us/mkl
 * IBM ESSL (Missing LAPACK routines are automatically added.),
   https://www.ibm.com/support/knowledgecenter/en/SSFHY8

The compilation requires at least CMAKE 3.1.


## Configuration and Installation

CMAKE is used to configure the source code. Out-of-source builds are prefered.

The default installation of qrupdate-ng is done by

    mkdir -p build && cd build
    cmake ../ -DCMAKE_INSTALL_PREFIX=/your/installation/path/
    make
    make install

The installation can adjusted with the help of the CMAKE options. Beside the
standard options of CMAKE, the following ones are supported:

| **Option**                     | **Description**                             |
|--------------------------------|---------------------------------------------|
| `-DEBUG=OFF/ON`                | Enable the debug build.                     |
| `-DHOSTOPT=OFF/ON`             | Enable host specific compiler flags.        |
| `-DFORTRAN_BOUND_CHECK=OFF/ON` | Enable the runtime bound checker.           |
| `-DFORTRAN_SANITIZE=OFF/ON`    | Enable the runtime sanitizer.               |
| `-DBUILD_SHARED_LIBS=ON/OFF`   | Enable building of shared libraries.        |

The `FORTRAN_BOUND_CHECK` option is only supported for gfortran and ifort
compilers. The `FORTRAN_SANITIZE` option can only be used with the gfortran
compiler. If the host optimizations are enabled, the build might not be portable
to other computers.

The qrupdate-ng code includes a testsuite which can be executed by `make test`.

### Selecting individual BLAS and LAPACK libraries

By default CMAKE uses the first BLAS and LAPACK library it finds on the system.
If a special selection is required, one can either specifiy the `BLA_VENDOR`
vendor option of CMAKE's `FindBLAS.cmake` module or specify the BLAS and LAPACK
libraries completly on ones own.

Possible values for the `BLA_VENDOR` options can be found in the help of the
`FindBLAS.cmake` module and can be obtained by `cmake --help-module FindBLAS`.

A custom BLAS and LAPACK library can be specified by setting the
`BLAS_LIBRRARIES` and `LAPACK_LIBRARIES` variables. The variable contain a
semicolon separated list of libraries that are required to provide the BLAS and
LAPACK functionality. If a libraries contains both, BLAS and LAPACK, it needs
to be set in both variables. An example to use a self-compiled reference BLAS
and LAPACK could yield the following CMAKE call:

    cmake ../ -DBLAS_LIBRARIES=/home/user/software/libblas.a \
              -DLAPACK_LIBRARIES=/home/user/software/liblapack.a

## Cross Compiling for Windows
The build system supports cross compiling from Unix-like operating systems to
Windows. Therefore, MingW64 (https://www.mingw-w64.org/) and
Wine (https://www.winehq.org/) need to be present on the system. In case of
Ubuntu 20.04 the following packages are required:
* `mingw-w64`, `mingw-w64-i686-dev` , `mingw-w64-x86-64-dev`
* `gfortran-mingw-w64-i686`, `gfortran-mingw-w64-x86-64`
* `wine32`, `wine64`
* `wget`

Furthermore, BLAS and LAPACK are required in the MingW installation. If this is
not the case, the reference BLAS and LAPACK library can be installed via
```shell
 bash ./tools/install-lapack-mingw-i686.sh
```
for the Windows 32-bit environment, and  using
```shelll
bash ./tools/install-lapack-mingw-x86_64.sh
```
for the Windows 64-bit environment.

The Windows 32-bit  library can then be compiled using
```shell
cmake -S . -B build.win32 \
      -DCMAKE_TOOLCHAIN_FILE=$(pwd)/cmake/mingw-w32-i686.cmake \
      -DBUILD_SHARED_LIBS=ON \
      -DCMAKE_INSTALL_PREFIX=$(pwd)/install-win32
make -C build.win32 all
(cd build.win32; ctest -V )
make -C build.win32 install
```

In case of the 64-bit Windows library, the build process looks like
```shell
cmake -S . -B build.win64 \
      -DCMAKE_TOOLCHAIN_FILE=$(pwd)/cmake/mingw-w64-x86_64.cmake \
      -DBUILD_SHARED_LIBS=ON \
      -DCMAKE_INSTALL_PREFIX=$(pwd)/install-win32
make -C build.win64 all
(cd build.win64; ctest -V )
make -C build.win64 install
```


