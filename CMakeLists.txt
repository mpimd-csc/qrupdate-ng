CMAKE_MINIMUM_REQUIRED(VERSION 3.18.0)
CMAKE_POLICY(SET CMP0048 NEW)
PROJECT(qrupdate-ng VERSION 1.2.0 LANGUAGES Fortran)
ENABLE_TESTING()

# Options
OPTION(BUILD_SHARED_LIBS "Build shared libraries" ON)
OPTION(BUILD_STATIC_LIBS "Build static libraries" ON)
OPTION(DEBUG "Enable Debug Symbol generation" OFF)

if (NOT ( BUILD_SHARED_LIBS OR BUILD_STATIC_LIBS  ) )
    MESSAGE(FATAL_ERROR "Either BUILD_SHARED_LIBS or BUILD_STATIC_LIBS must be set to ON.")
endif()

# SET(CMAKE_POSITION_INDEPENDENT_CODE TRUE)
IF ( DEBUG STREQUAL ON )
	SET (CMAKE_BUILD_TYPE "Debug")
	ADD_DEFINITIONS(-DDEBUG)
ENDIF()

IF(NOT CMAKE_BUILD_TYPE)
	SET (CMAKE_BUILD_TYPE "Release")
ENDIF()

# Output directories
IF(CMAKE_CURRENT_SOURCE_DIR STREQUAL CMAKE_SOURCE_DIR)
   SET(CMAKE_ARCHIVE_OUTPUT_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/lib)
   SET(CMAKE_LIBRARY_OUTPUT_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/lib)
   SET(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/bin)
ENDIF()


# Required to work with Intel MKL
IF (CMAKE_Fortran_COMPILER_ID STREQUAL "Intel" OR BLA_VENDOR MATCHES "Intel.*")
	ENABLE_LANGUAGE(C)
ENDIF()

# Add addtional CMAKE Paths.
LIST(APPEND CMAKE_MODULE_PATH "${PROJECT_SOURCE_DIR}/cmake")

# enable code coverage
find_package(codecov)

SET(Fortran_MODULE_NAME qrupdate)
INCLUDE(FortranCompilerSettings)
INCLUDE(GNUInstallDirs)
INCLUDE(FortranModuleDir)


# IBM XLF Compilation
IF( CMAKE_Fortran_COMPILER_ID STREQUAL "XL")
	SET(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -qfixed -qnosave")
ENDIF()

# Find Doxygen
find_package(Doxygen OPTIONAL_COMPONENTS dot)

if (INTEGER8 STREQUAL ON)
    cmake_minimum_required(VERSION 3.22)
    SET(BLA_SIZEOF_INTEGER 8)
endif()

FIND_PACKAGE(BLAS REQUIRED)
MESSAGE(STATUS "BLAS Vendor: ${BLA_VENDOR}")
IF(NOT BLA_VENDOR STREQUAL "IBMESSL")
	FIND_PACKAGE(LAPACK REQUIRED)
ELSE()
	SET(LAPACK_MISSING_TARGET_NAME lapack_missing)
	INCLUDE(${CMAKE_SOURCE_DIR}/lapack-missing/LapackMissing.cmake)
	SET(LAPACK_FOUND TRUE)
	SET(LAPACK_LIBRARIES ${BLAS_LIBRARIES})
ENDIF()


ADD_SUBDIRECTORY(src)
ADD_SUBDIRECTORY(test)

# Add documentation target
if(DOXYGEN_FOUND)
    # Configure Doxygen to extract all functions including undocumented ones
    # and to browse source files
    set(DOXYGEN_EXTRACT_ALL YES)
    set(DOXYGEN_EXTRACT_STATIC YES)
    set(DOXYGEN_SOURCE_BROWSER YES)
    set(DOXYGEN_GENERATE_TREEVIEW YES)
    # Optimize output for Fortran
    set(DOXYGEN_OPTIMIZE_OUTPUT_FOR_FORTRAN YES)
    set(DOXYGEN_OUTPUT_DIRECTORY "${CMAKE_CURRENT_BINARY_DIR}/doc")
    # Set input directories and files
    set(DOXYGEN_INPUT_DIR
        ${CMAKE_CURRENT_SOURCE_DIR}/README.md
        ${CMAKE_CURRENT_SOURCE_DIR}/INSTALL.md
        ${CMAKE_CURRENT_SOURCE_DIR}/CHANGELOG.md
        ${CMAKE_CURRENT_SOURCE_DIR}/src/categories.doxy
        ${CMAKE_CURRENT_SOURCE_DIR}/LICENSE ${CMAKE_CURRENT_SOURCE_DIR}/src)
    set(DOXYGEN_USE_MDFILE_AS_MAINPAGE ${CMAKE_CURRENT_SOURCE_DIR}/README.md)
    set(DOXYGEN_MARKDOWN_SUPPORT YES)
    set(DOXYGEN_EXCLUDE_PATTERNS "*/build*/*")

    # Configure Doxygen
    doxygen_add_docs(doc
        ${DOXYGEN_INPUT_DIR}
        ALL
        COMMENT "Generate documentation"
    )
    install(DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/doc/html/
            TYPE DOC)
endif()


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

coverage_evaluate()
