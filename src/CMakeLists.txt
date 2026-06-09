SET(SRC
    caxcpy.f90
    cch1dn.f90
    cch1up.f90
    cchdex.f90
    cchinx.f90
    cchshx.f90
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

IF (TARGET lapack_missing)
	ADD_LIBRARY(qrupdate ${SRC} $<TARGET_OBJECTS:lapack_missing>)
ELSE()
	ADD_LIBRARY(qrupdate ${SRC})
ENDIF()
TARGET_LINK_LIBRARIES(qrupdate ${LAPACK_LIBRARIES} ${BLAS_LIBRARIES})
SET_PROPERTY(TARGET qrupdate PROPERTY VERSION ${PROJECT_VERSION})
SET_PROPERTY(TARGET qrupdate PROPERTY SOVERSION ${PROJECT_VERSION_MAJOR})

INSTALL(TARGETS qrupdate
        EXPORT  qrupdatetargets
        LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR}
        ARCHIVE DESTINATION ${CMAKE_INSTALL_LIBDIR}
        RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR}
        )
