CMAKE_MINIMUM_REQUIRED(VERSION 3.1.0)
CMAKE_POLICY(SET CMP0048 NEW)
PROJECT(qrupdate-ng VERSION 1.1.3 LANGUAGES Fortran)
ENABLE_TESTING()

# Options
OPTION(BUILD_SHARED_LIBS "Build shared libraries" ON)
OPTION(DEBUG "Enable Debug Symbol generation" OFF)

SET(CMAKE_POSITION_INDEPENDENT_CODE TRUE)
IF ( DEBUG STREQUAL ON )
	SET (CMAKE_BUILD_TYPE "Debug")
	ADD_DEFINITIONS(-DDEBUG)
ENDIF()

IF(NOT CMAKE_BUILD_TYPE)
	SET (CMAKE_BUILD_TYPE "Release")
ENDIF()

# Add addtional CMAKE Paths.
LIST(APPEND CMAKE_MODULE_PATH "${PROJECT_SOURCE_DIR}/cmake")


INCLUDE(FortranCompilerSettings)
INCLUDE(GNUInstallDirs)

FIND_PACKAGE(BLAS REQUIRED)
FIND_PACKAGE(LAPACK REQUIRED)

ADD_SUBDIRECTORY(src)
ADD_SUBDIRECTORY(test)

INCLUDE(CMakePackageConfigHelpers)
WRITE_BASIC_PACKAGE_VERSION_FILE(
    qrupdateConfigVersion.cmake
    VERSION ${PROJECT_VERSION}
    COMPATIBILITY AnyNewerVersion
 )
INSTALL(EXPORT qrupdatetargets
        FILE qrupdateTargets.cmake
        NAMESPACE qrupdate::
        DESTINATION ${CMAKE_INSTALL_LIBDIR}/cmake/qrupdate
)

CONFIGURE_FILE(qrupdate-config.cmake.in qrupdate-config.cmake @ONLY)
INSTALL(FILES "${CMAKE_CURRENT_BINARY_DIR}/qrupdate-config.cmake"
              "${CMAKE_CURRENT_BINARY_DIR}/qrupdateConfigVersion.cmake"
              DESTINATION ${CMAKE_INSTALL_LIBDIR}/cmake/qrupdate
        )

CONFIGURE_FILE(qrupdate.pc.in qrupdate.pc @ONLY)
INSTALL(FILES "${CMAKE_CURRENT_BINARY_DIR}/qrupdate.pc"
              DESTINATION "${CMAKE_INSTALL_LIBDIR}/pkgconfig")

MESSAGE(STATUS "------------ Build Information --------------")
MESSAGE(STATUS "CMAKE_Fortran_FLAGS:  ${CMAKE_Fortran_FLAGS}")
MESSAGE(STATUS "CMAKE_Fortran_FLAGS_DEBUG:  ${CMAKE_Fortran_FLAGS_DEBUG}")
MESSAGE(STATUS "CMAKE_Fortran_FLAGS_RELEASE:  ${CMAKE_Fortran_FLAGS_RELEASE}")
MESSAGE(STATUS "BLAS Libraries: ${BLAS_LIBRARIES}")
MESSAGE(STATUS "LAPACK Libraries: ${LAPACK_LIBRARIES}")
MESSAGE(STATUS "Build Type: ${CMAKE_BUILD_TYPE}")
MESSAGE(STATUS "Debug: ${DEBUG}")
MESSAGE(STATUS "Shared Libs: ${BUILD_SHARED_LIBS}")
MESSAGE(STATUS "---------------------------------------------")

