SET(SRC
    caxcpy.f
    cch1dn.f
    cch1up.f
    cchdex.f
    cchinx.f
    cchshx.f
    cgqvec.f
    clu1up.f
    clup1up.f
    cqhqr.f
    cqr1up.f
    cqrdec.f
    cqrder.f
    cqrinc.f
    cqrinr.f
    cqrot.f
    cqrqh.f
    cqrshc.f
    cqrtv1.f
    dch1dn.f
    dch1up.f
    dchdex.f
    dchinx.f
    dchshx.f
    dgqvec.f
    dlu1up.f
    dlup1up.f
    dqhqr.f
    dqr1up.f
    dqrdec.f
    dqrder.f
    dqrinc.f
    dqrinr.f
    dqrot.f
    dqrqh.f
    dqrshc.f
    dqrtv1.f
    sch1dn.f
    sch1up.f
    schdex.f
    schinx.f
    schshx.f
    sgqvec.f
    slu1up.f
    slup1up.f
    sqhqr.f
    sqr1up.f
    sqrdec.f
    sqrder.f
    sqrinc.f
    sqrinr.f
    sqrot.f
    sqrqh.f
    sqrshc.f
    sqrtv1.f
    zaxcpy.f
    zch1dn.f
    zch1up.f
    zchdex.f
    zchinx.f
    zchshx.f
    zgqvec.f
    zlu1up.f
    zlup1up.f
    zqhqr.f
    zqr1up.f
    zqrdec.f
    zqrder.f
    zqrinc.f
    zqrinr.f
    zqrot.f
    zqrqh.f
    zqrshc.f
    zqrtv1.f
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
