# Sample toolchain file for building for Windows from an Ubuntu Linux system.
#
# Typical usage:
#    *) install cross compiler: `sudo apt-get install mingw-w64` or `brew install mingw-w64` on macOS
#    *) cmake -DCMAKE_TOOLCHAIN_FILE=~/mingw-w64-x86_64.cmake -G Ninja -B build -S .
#    *) ninja -C build


set(CMAKE_SYSTEM_NAME Windows)
set(TOOLCHAIN_PREFIX i686-w64-mingw32)

# cross compilers to use for C, C++ and Fortran
set(CMAKE_C_COMPILER ${TOOLCHAIN_PREFIX}-gcc)
set(CMAKE_CXX_COMPILER ${TOOLCHAIN_PREFIX}-g++)
set(CMAKE_Fortran_COMPILER ${TOOLCHAIN_PREFIX}-gfortran)
set(CMAKE_RC_COMPILER ${TOOLCHAIN_PREFIX}-windres)

# target environment on the build host system
set(CMAKE_FIND_ROOT_PATH /usr/${TOOLCHAIN_PREFIX} ${CMAKE_SOURCE_DIR}/tools/${TOOLCHAIN_PREFIX})

# modify default behavior of FIND_XXX() commands
set(CMAKE_FIND_ROOT_PATH_MODE_PROGRAM NEVER)
set(CMAKE_FIND_ROOT_PATH_MODE_LIBRARY ONLY)
set(CMAKE_FIND_ROOT_PATH_MODE_INCLUDE ONLY)

set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS}")
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS}")
set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS}")

set(Fortran_LOCAL_FLAGS "-static")
set(C_LOCAL_FLAGS "-static")
set(CXX_LOCAL_FLAGS "-static")

#set(CMAKE_SHARED_LINKER_FLAGS "${CMAKE_SHARED_LINKER_FLAGS} -static-libgcc -Wl,-Bstatic -lgfortran -lquadmath -Wl,-Bdynamic")
#set(CMAKE_EXE_LINKER_FLAGS " -static-libgcc -Wl,-Bstatic -lgfortran -lquadmath -Wl,-Bdynamic ${CMAKE_EXE_LINKER_FLAGS} ")

