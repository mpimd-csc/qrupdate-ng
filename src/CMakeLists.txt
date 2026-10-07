SET(SRC
    qrupdate.f90
    qrupdate_blas.f90
    qrupdate_error.f90
    caxcpy.f90
    cch1dn.f90
    cch1up.f90
    cchdex.f90
    cchinx.f90
    cchshx.f90
    cdotc.f90
    cdotu.f90
    cgqvec.f90
    clu1up.f90
    clup1up.f90
    cqhqr.f90
    cqr1up.f90
    cqrdec.f90
    cqrder.f90
    cqrinc.f90
    cqrinr.f90
    cqrot.f90
    cqrqh.f90
    cqrshc.f90
    cqrtv1.f90
    dch1dn.f90
    dch1up.f90
    dchdex.f90
    dchinx.f90
    dchshx.f90
    dgqvec.f90
    dlu1up.f90
    dlup1up.f90
    dqhqr.f90
    dqr1up.f90
    dqrdec.f90
    dqrder.f90
    dqrinc.f90
    dqrinr.f90
    dqrot.f90
    dqrqh.f90
    dqrshc.f90
    dqrtv1.f90
    sch1dn.f90
    sch1up.f90
    schdex.f90
    schinx.f90
    schshx.f90
    sgqvec.f90
    slu1up.f90
    slup1up.f90
    sqhqr.f90
    sqr1up.f90
    sqrdec.f90
    sqrder.f90
    sqrinc.f90
    sqrinr.f90
    sqrot.f90
    sqrqh.f90
    sqrshc.f90
    sqrtv1.f90
    zaxcpy.f90
    zch1dn.f90
    zch1up.f90
    zchdex.f90
    zchinx.f90
    zchshx.f90
    zdotc.f90
    zdotu.f90
    zgqvec.f90
    zlu1up.f90
    zlup1up.f90
    zqhqr.f90
    zqr1up.f90
    zqrdec.f90
    zqrder.f90
    zqrinc.f90
    zqrinr.f90
    zqrot.f90
    zqrqh.f90
    zqrshc.f90
    zqrtv1.f90
)

SET(ITARGETS)
#
# C Interface
#


# C Interface and Fortran Interop
if(QRUPDATE_C_INTERFACE)
    include(FortranCInterface)
    FortranCInterface_HEADER(
        ${CMAKE_CURRENT_BINARY_DIR}/qrupdate_fortran_mangle.h
        MACRO_NAMESPACE QRUPDATE_FORTRAN_
    )
    FortranCInterface_VERIFY()

    configure_file(qrupdate_config.h.in ${CMAKE_CURRENT_BINARY_DIR}/qrupdate_config.h @ONLY)
    set(SRC ${SRC}
        c_interface/caxcpy.c
        c_interface/cch1dn.c
        c_interface/cch1up.c
        c_interface/cchdex.c
        c_interface/cchinx.c
        c_interface/cchshx.c
        c_interface/cgqvec.c
        c_interface/clu1up.c
        c_interface/clup1up.c
        c_interface/cqhqr.c
        c_interface/cqr1up.c
        c_interface/cqrdec.c
        c_interface/cqrder.c
        c_interface/cqrinc.c
        c_interface/cqrinr.c
        c_interface/cqrot.c
        c_interface/cqrqh.c
        c_interface/cqrshc.c
        c_interface/cqrtv1.c
        c_interface/dch1dn.c
        c_interface/dch1up.c
        c_interface/dchdex.c
        c_interface/dchinx.c
        c_interface/dchshx.c
        c_interface/dgqvec.c
        c_interface/dlu1up.c
        c_interface/dlup1up.c
        c_interface/dqhqr.c
        c_interface/dqr1up.c
        c_interface/dqrdec.c
        c_interface/dqrder.c
        c_interface/dqrinc.c
        c_interface/dqrinr.c
        c_interface/dqrot.c
        c_interface/dqrqh.c
        c_interface/dqrshc.c
        c_interface/dqrtv1.c
        c_interface/sch1dn.c
        c_interface/sch1up.c
        c_interface/schdex.c
        c_interface/schinx.c
        c_interface/schshx.c
        c_interface/sgqvec.c
        c_interface/slu1up.c
        c_interface/slup1up.c
        c_interface/sqhqr.c
        c_interface/sqr1up.c
        c_interface/sqrdec.c
        c_interface/sqrder.c
        c_interface/sqrinc.c
        c_interface/sqrinr.c
        c_interface/sqrot.c
        c_interface/sqrqh.c
        c_interface/sqrshc.c
        c_interface/sqrtv1.c
        c_interface/zaxcpy.c
        c_interface/zch1dn.c
        c_interface/zch1up.c
        c_interface/zchdex.c
        c_interface/zchinx.c
        c_interface/zchshx.c
        c_interface/zgqvec.c
        c_interface/zlu1up.c
        c_interface/zlup1up.c
        c_interface/zqhqr.c
        c_interface/zqr1up.c
        c_interface/zqrdec.c
        c_interface/zqrder.c
        c_interface/zqrinc.c
        c_interface/zqrinr.c
        c_interface/zqrot.c
        c_interface/zqrqh.c
        c_interface/zqrshc.c
        c_interface/zqrtv1.c)


endif()


IF (TARGET lapack_missing)
    if (BUILD_SHARED_LIBS)
        ADD_LIBRARY(qrupdate SHARED ${SRC} $<TARGET_OBJECTS:lapack_missing>)
    endif()
    if (BUILD_STATIC_LIBS)
        ADD_LIBRARY(qrupdate_static STATIC ${SRC} $<TARGET_OBJECTS:lapack_missing>)
    endif()
ELSE()
    if (BUILD_SHARED_LIBS)
        ADD_LIBRARY(qrupdate SHARED ${SRC})
    endif()
    if (BUILD_STATIC_LIBS)
        ADD_LIBRARY(qrupdate_static STATIC ${SRC})
    endif()
ENDIF()

if (QRUPDATE_C_INTERFACE)
    include(GenerateExportHeader)

    if (BUILD_SHARED_LIBS)
        generate_export_header(
            qrupdate
            BASE_NAME QRUPDATE
            EXPORT_FILE_NAME qrupdate_export.h
            STATIC_DEFINE QRUPDATE_STATIC
            INCLUDE_GUARD_NAME QRUPDATE_EXPORT_H
        )
    else()
        generate_export_header(
            qrupdate_static
            BASE_NAME QRUPDATE
            STATIC_DEFINE QRUPDATE_STATIC
            EXPORT_FILE_NAME qrupdate_export.h
            INCLUDE_GUARD_NAME QRUPDATE_EXPORT_H
        )
    endif()


endif()

if (BUILD_SHARED_LIBS)
    TARGET_LINK_LIBRARIES(qrupdate PRIVATE LAPACK::LAPACK BLAS::BLAS)
    add_library(qrupdate::qrupdate ALIAS qrupdate)
    SET_PROPERTY(TARGET qrupdate PROPERTY VERSION ${PROJECT_VERSION})
    SET_PROPERTY(TARGET qrupdate PROPERTY SOVERSION ${PROJECT_VERSION_MAJOR})

    SET_TARGET_PROPERTIES(qrupdate PROPERTIES Fortran_MODULE_DIRECTORY "${CMAKE_CURRENT_BINARY_DIR}/fortran-modules")
    SET_TARGET_PROPERTIES(qrupdate PROPERTIES WINDOWS_EXPORT_ALL_SYMBOLS ON)

    IF (CMAKE_Fortran_COMPILER_ID STREQUAL "LLVMFlang")
        SET_TARGET_PROPERTIES(qrupdate PROPERTIES COMPILE_FLAGS "-fPIC")
    ENDIF()


    TARGET_INCLUDE_DIRECTORIES(qrupdate
        PUBLIC
        $<BUILD_INTERFACE:${CMAKE_CURRENT_BINARY_DIR}/fortran-modules>
        $<BUILD_INTERFACE:${CMAKE_CURRENT_BINARY_DIR}>
        $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}>
        $<INSTALL_INTERFACE:include/qrupdate>
    )

    ADD_COVERAGE(qrupdate)
    list(APPEND ITARGETS qrupdate)

    if (QRUPDATE_C_INTERFACE)
        SET_PROPERTY(TARGET qrupdate PROPERTY
            PUBLIC_HEADER
            "${CMAKE_CURRENT_BINARY_DIR}/qrupdate_config.h"
            "${CMAKE_CURRENT_BINARY_DIR}/qrupdate_export.h"
            "${CMAKE_CURRENT_BINARY_DIR}/qrupdate_fortran_mangle.h"
            "${CMAKE_CURRENT_SOURCE_DIR}/qrupdate_fortran_types.h"
            "${CMAKE_CURRENT_SOURCE_DIR}/qrupdate_error.h"
            "${CMAKE_CURRENT_SOURCE_DIR}/qrupdate.h")
    endif()

    if (CMAKE_SYSTEM_NAME MATCHES "Windows" OR CYGWIN)
	set_target_properties(qrupdate PROPERTIES LINK_FLAGS "-Wl,--export-all-symbols")
    endif()
endif()

if (BUILD_STATIC_LIBS)
    TARGET_LINK_LIBRARIES(qrupdate_static PRIVATE LAPACK::LAPACK BLAS::BLAS)
    add_library(qrupdate::qrupdate_static ALIAS qrupdate_static)

    SET_TARGET_PROPERTIES(qrupdate_static PROPERTIES Fortran_MODULE_DIRECTORY "${CMAKE_CURRENT_BINARY_DIR}/fortran-modules-static")
    SET_TARGET_PROPERTIES(qrupdate_static PROPERTIES OUTPUT_NAME qrupdate)
    TARGET_COMPILE_DEFINITIONS(qrupdate_static PUBLIC QRUPDATE_STATIC)
    TARGET_INCLUDE_DIRECTORIES(qrupdate_static
        PUBLIC
        $<BUILD_INTERFACE:${CMAKE_CURRENT_BINARY_DIR}/fortran-modules-static>
        $<BUILD_INTERFACE:${CMAKE_CURRENT_BINARY_DIR}>
        $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}>
        $<INSTALL_INTERFACE:include/qrupdate>
    )
    if (QRUPDATE_C_INTERFACE AND NOT BUILD_SHARED_LIBS)
        SET_PROPERTY(TARGET qrupdate_static PROPERTY
            PUBLIC_HEADER
            "${CMAKE_CURRENT_BINARY_DIR}/qrupdate_config.h"
            "${CMAKE_CURRENT_BINARY_DIR}/qrupdate_export.h"
            "${CMAKE_CURRENT_BINARY_DIR}/qrupdate_fortran_mangle.h"
            "${CMAKE_CURRENT_SOURCE_DIR}/qrupdate_fortran_types.h"
            "${CMAKE_CURRENT_SOURCE_DIR}/qrupdate_error.h"
            "${CMAKE_CURRENT_SOURCE_DIR}/qrupdate.h")
    endif()

    if (QRUPDATE_STATIC_WITH_PIC)
        SET_TARGET_PROPERTIES(qrupdate_static PROPERTIES POSITION_INDEPENDENT_CODE ON)
        IF (CMAKE_Fortran_COMPILER_ID STREQUAL "LLVMFlang")
            SET_TARGET_PROPERTIES(qrupdate_static PROPERTIES COMPILE_FLAGS "-fPIC")
        ENDIF()
    endif()
    list(APPEND ITARGETS qrupdate_static )
    ADD_COVERAGE(qrupdate_static)
endif()

if (BUILD_SHARED_LIBS)
INSTALL(DIRECTORY "${CMAKE_CURRENT_BINARY_DIR}/fortran-modules/" DESTINATION ${CMAKE_INSTALL_FULL_MODULEDIR} )
else()
INSTALL(DIRECTORY "${CMAKE_CURRENT_BINARY_DIR}/fortran-modules-static/" DESTINATION ${CMAKE_INSTALL_FULL_MODULEDIR} )
endif()

INSTALL(TARGETS ${ITARGETS}
        EXPORT  qrupdatetargets
        LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR}
        ARCHIVE DESTINATION ${CMAKE_INSTALL_LIBDIR}
        RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR}
        PUBLIC_HEADER DESTINATION ${CMAKE_INSTALL_FULL_MODULEDIR}
    )

