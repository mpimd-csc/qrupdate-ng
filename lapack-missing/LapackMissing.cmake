SET(CURDIR ${CMAKE_CURRENT_LIST_DIR})

INCLUDE(CheckFortranFunctionExists)

IF(NOT LAPACK_MISSING_TARGET_NAME)
    MESSAGE(FATAL_ERROR "LAPACK_MISSING_TARGET_NAME not set.")
ENDIF()

IF(NOT LAPACK_MISSING_LAPACK_LIBRARIES)
    SET(LAPACK_MISSING_LAPACK_LIBRARIES ${LAPACK_LIBRARIES})
ENDIF()
IF(NOT LAPACK_MISSING_BLAS_LIBRARIES)
    SET(LAPACK_MISSING_BLAS_LIBRARIES ${BLAS_LIBRARIES})
ENDIF()

MACRO(CHECKLAPACKROUTINE ROUTINE VAR_TO_SET)
    SET(_X ${CMAKE_REQUIRED_LIBRARIES})
    SET(CMAKE_REQUIRED_LIBRARIES ${LAPACK_MISSING_LAPACK_LIBRARIES} ${LAPACK_MISSING_BLAS_LIBRARIES})

    CHECK_FORTRAN_FUNCTION_EXISTS(${ROUTINE} ${VAR_TO_SET})
    SET(CMAKE_REQUIRED_LIBRARIES ${_X})
ENDMACRO()

MACRO(ADD_LAPACK_SRC symbol)
CheckLapackRoutine(${symbol} HAVE_${symbol})
IF (NOT HAVE_${symbol})
	MESSAGE(STATUS "Add separate ${symbol} to library.")
	FOREACH(S ${ARGN})
		MESSAGE(STATUS "  --> src/${S}")
        SET(SRC ${SRC} "${CURDIR}/src/${S}")
	ENDFOREACH()
ENDIF()
SET(HAVE_${symbol} ${HAVE_${symbol}} CACHE BOOL "Lapack Symbol ${symbol} exists")
ENDMACRO()


SET (SRC)


ADD_LAPACK_SRC(daxpby   daxpby.f saxpby.f caxpby.f zaxpby.f)
ADD_LAPACK_SRC(dgebak   dgebak.f sgebak.f cgebak.f zggbak.f)
ADD_LAPACK_SRC(dgebal   dgebal.f sgebal.f cgebal.f zgebal.f)
ADD_LAPACK_SRC(dgebd2   dgebd2.f sgebd2.f cgebd2.f zgebd2.f)
ADD_LAPACK_SRC(dgebrd   dgebrd.f sgebrd.f cgebrd.f zgebrd.f)
ADD_LAPACK_SRC(dgecon   dgecon.f )
ADD_LAPACK_SRC(dgees    dgees.f  sgees.f  cgees.f  zgees.f)
ADD_LAPACK_SRC(dgegs    dgegs.f sgegs.f cgegs.f zgegs.f)
ADD_LAPACK_SRC(dgehd2   dgehd2.f sgehd2.f cgehd2.f zgehd2.f)
ADD_LAPACK_SRC(dgehrd   dgehrd.f sgehrd.f cgehrd.f zgehrd.f)
ADD_LAPACK_SRC(dgelq2   dgelq2.f sgelq2.f cgelq2.f zgelq2.f)
ADD_LAPACK_SRC(dgelqf   dgelqf.f sgelqf.f cgelqf.f zgelqf.f)
ADD_LAPACK_SRC(dgelqt   dgelqt.f sgelqt.f cgelqt.f zgelqt.f)
ADD_LAPACK_SRC(dgelqt3  dgelqt3.f sgelqt3.f cgelqt3.f zgelqt3.f)
ADD_LAPACK_SRC(dgemlqt  dgemlqt.f sgemlqt.f cgemlqt.f zgemlqt.f)
ADD_LAPACK_SRC(dgemqr   	dgemqr.f sgemqr.f cgemqr.f zgemqr.f)
ADD_LAPACK_SRC(dgemqrt  dgemqrt.f sgemqrt.f cgemqrt.f zgemqrt.f)
ADD_LAPACK_SRC(dgeqp3   dgeqp3.f sgeqp3.f cgeqp3.f zgeqp3.f)
ADD_LAPACK_SRC(dgeqp3   dgeqp3.f sgeqp3.f cgeqp3.f zgeqp3.f)
ADD_LAPACK_SRC(dgeqpf   dgeqpf.f sgeqpf.f cgeqpf.f zgeqpf.f)
ADD_LAPACK_SRC(dgeqr   	dgeqr.f sgeqr.f cgeqr.f zgeqr.f)
ADD_LAPACK_SRC(dgeqr2   dgeqr2.f sgeqr2.f cgeqr2.f zgeqr2.f)
ADD_LAPACK_SRC(dgeqrf   dgeqrf.f )
ADD_LAPACK_SRC(dgeqrt   dgeqrt.f sgeqrt.f cgeqrt.f zgeqrt.f)
ADD_LAPACK_SRC(dgeqrt2  dgeqrt2.f sgeqrt2.f cgeqrt2.f zgeqrt2.f)
ADD_LAPACK_SRC(dgeqrt3  dgeqrt3.f sgeqrt3.f cgeqrt3.f zgeqrt3.f)
ADD_LAPACK_SRC(dgerq2   dgerq2.f sgerq2.f cgerq2.f zgerq2.f)
ADD_LAPACK_SRC(dgerqf   dgerqf.f)
ADD_LAPACK_SRC(dgesc2   dgesc2.f sgesc2.f cgesc2.f zgesc2.f)
ADD_LAPACK_SRC(dgetc2   dgetc2.f sgetc2.f cgetc2.f zgetc2.f)
ADD_LAPACK_SRC(dggbak   dggbak.f sggbak.f cggbak.f zggbak.f)
ADD_LAPACK_SRC(dggbal   dggbal.f sggbal.f cggbal.f zggbal.f)
ADD_LAPACK_SRC(dgges    dgges.f  sgges.f cgges.f zgges.f)
ADD_LAPACK_SRC(dgges3   dgges3.f)
ADD_LAPACK_SRC(dgghd3   dgghd3.f)
ADD_LAPACK_SRC(dgghrd   dgghrd.f sgghrd.f cgghrd.f zgghrd.f)
ADD_LAPACK_SRC(dhgeqz   dhgeqz.f shgeqz.f chgeqz.f zhgeqz.f)
ADD_LAPACK_SRC(dhseqr   dhseqr.f shseqr.f chseqr.f zhseqr.f)
ADD_LAPACK_SRC(disnan   disnan.f sisnan.f)
ADD_LAPACK_SRC(dlabad   dlabad.f slabad.f)
ADD_LAPACK_SRC(dlabrd   dlabrd.f slabrd.f clabrd.f zlabrd.f)
ADD_LAPACK_SRC(dlacn2   dlacn2.f slacn2.f clacn2.f zlacn2.f)
ADD_LAPACK_SRC(dlacon   dlacon.f slacon.f clacon.f zlacon.f)
ADD_LAPACK_SRC(dlacpy   dlacpy.f slacpy.f clacpy.f zlacpy.f)
ADD_LAPACK_SRC(dladiv   dladiv.f sladiv.f cladiv.f zladiv.f)
ADD_LAPACK_SRC(dlae2    dlae2.f slae2.f)
ADD_LAPACK_SRC(dlaev2   dlaev2.f slaev2.f claev2.f zlaev2.f)
ADD_LAPACK_SRC(dlaexc   dlaexc.f slaexc.f)
ADD_LAPACK_SRC(dlag2    dlag2.f  slag2.f)
ADD_LAPACK_SRC(dlagv2   dlagv2.f slagv2.f )
ADD_LAPACK_SRC(dlahqr   dlahqr.f slahqr.f clahqr.f zlahqr.f)
ADD_LAPACK_SRC(dlahr2   dlahr2.f slahr2.f clahr2.f zlahr2.f)
ADD_LAPACK_SRC(dlaic1   dlaic1.f slaic1.f claic1.f zlaic1.f)
ADD_LAPACK_SRC(dlaisnan dlaisnan.f slaisnan.f)
ADD_LAPACK_SRC(dlaln2 	dlaln2.f slaln2.f)
ADD_LAPACK_SRC(dlamch   dlamch.f slamch.f)
ADD_LAPACK_SRC(dlamtsqr   	dlamtsqr.f slamtsqr.f clamtsqr.f zlamtsqr.f)
ADD_LAPACK_SRC(dlanhs   dlanhs.f slanhs.f clanhs.f zlanhs.f)
ADD_LAPACK_SRC(dlanst   dlanst.f  slanst.f)
ADD_LAPACK_SRC(dlanv2   dlanv2.f slanv2.f)
ADD_LAPACK_SRC(dlapy2   dlapy2.f slapy2.f slapy3.f dlapy3.f)
ADD_LAPACK_SRC(dlaqp2   dlaqp2.f slaqp2.f claqp2.f zlaqp2.f)
ADD_LAPACK_SRC(dlaqps   dlaqps.f slaqps.f claqps.f zlaqps.f)
ADD_LAPACK_SRC(dlaqr0   dlaqr0.f slaqr0.f claqr0.f zlaqr0.f)
ADD_LAPACK_SRC(dlaqr1   dlaqr1.f slaqr1.f claqr1.f zlaqr1.f)
ADD_LAPACK_SRC(dlaqr2   dlaqr2.f slaqr2.f claqr2.f zlaqr2.f)
ADD_LAPACK_SRC(dlaqr3   dlaqr3.f slaqr3.f claqr3.f zlaqr3.f)
ADD_LAPACK_SRC(dlaqr4   dlaqr4.f slaqr4.f claqr4.f zlaqr4.f)
ADD_LAPACK_SRC(dlaqr5   dlaqr5.f slaqr5.f claqr5.f zlaqr5.f)
ADD_LAPACK_SRC(dlarf    dlarf.f slarf.f clarf.f zlarf.f)
ADD_LAPACK_SRC(dlarfb   dlarfb.f slarfb.f clarfb.f zlarfb.f)
ADD_LAPACK_SRC(dlarfg   dlarfg.f slarfg.f clarfg.f zlarfg.f)
ADD_LAPACK_SRC(dlarft   dlarft.f slarft.f clarft.f zlarft.f)
ADD_LAPACK_SRC(dlarfx   dlarfx.f slarfx.f clarfx.f zlarfx.f)
ADD_LAPACK_SRC(dlarnv   dlarnv.f slarnv.f clarnv.f zlarnv.f)
ADD_LAPACK_SRC(dlartg   dlartg.f slartg.f clartg.f zlartg.f  )
ADD_LAPACK_SRC(dlaruv   dlaruv.f slaruv.f)
ADD_LAPACK_SRC(dlascl   dlascl.f slascl.f clascl.f zlascl.f)
ADD_LAPACK_SRC(dlaset   dlaset.f slaset.f claset.f zlaset.f)
ADD_LAPACK_SRC(dlasr    dlasr.f slasr.f clasr.f zlasr.f)
ADD_LAPACK_SRC(dlasrt   dlasrt.f    slasrt.f)
ADD_LAPACK_SRC(dlassq   dlassq.f slassq.f classq.f zlassq.f)
ADD_LAPACK_SRC(dlasv2   dlasv2.f slasv2.f)
ADD_LAPACK_SRC(dlaswp   dlaswp.f slaswp.f claswp.f zlaswp.f)
ADD_LAPACK_SRC(dlasy2   dlasy2.f slasy2.f)
ADD_LAPACK_SRC(dlatdf   dlatdf.f slatdf.f clatdf.f zlatdf.f)
ADD_LAPACK_SRC(dlatrs   dlatrs.f slatrs.f clatrs.f zlatrs.f) 
ADD_LAPACK_SRC(dlatsqr   	dlatsqr.f slatsqr.f clatsqr.f zlatsqr.f)
ADD_LAPACK_SRC(dorg2r   dorg2r.f sorg2r.f)
ADD_LAPACK_SRC(dorghr   dorghr.f sorghr.f)
ADD_LAPACK_SRC(dorgqr   dorgqr.f sorgqr.f)
ADD_LAPACK_SRC(dorgr2   dorgr2.f sorgr2.f)
ADD_LAPACK_SRC(dorgrq   dorgrq.f)
ADD_LAPACK_SRC(dorgtsqr     dorgtsqr.f sorgtsqr.f cungtsqr.f zungtsqr.f)
ADD_LAPACK_SRC(dorm22   dorm22.f)
ADD_LAPACK_SRC(dorm2r   dorm2r.f sorm2r.f)
ADD_LAPACK_SRC(dormbr   dormbr.f sormbr.f)
ADD_LAPACK_SRC(dormhr   dormhr.f sormhr.f)
ADD_LAPACK_SRC(dorml2   dorml2.f sorml2.f)
ADD_LAPACK_SRC(dormlq   dormlq.f sormlq.f)
ADD_LAPACK_SRC(dormqr   dormqr.f sormqr.f)
ADD_LAPACK_SRC(dormr2   dormr2.f sormr2.f)
ADD_LAPACK_SRC(dormrq   dormrq.f)
ADD_LAPACK_SRC(dsteqr   dsteqr.f  ssteqr.f  csteqr.f zsteqr.f)
ADD_LAPACK_SRC(dtgevc   dtgevc.f stgevc.f ctgevc.f ztgevc.f)
ADD_LAPACK_SRC(dtgex2   dtgex2.f stgex2.f ctgex2.f ztgex2.f)
ADD_LAPACK_SRC(dtgexc   dtgexc.f stgexc.f ctgexc.f ztgexc.f )
ADD_LAPACK_SRC(dtgsen   dtgsen.f stgsen.f ctgsen.f ztgsen.f )
ADD_LAPACK_SRC(dtgsy2   dtgsy2.f stgsy2.f ctgsy2.f ztgsy2.f)
ADD_LAPACK_SRC(dtgsyl   dtgsyl.f stgsyl.f ctgsyl.f ztgsyl.f)
ADD_LAPACK_SRC(dtplqt  dtplqt.f stplqt.f ctplqt.f ztplqt.f)
ADD_LAPACK_SRC(dtplqt2 dtplqt2.f stplqt2.f ctplqt2.f ztplqt2.f)
ADD_LAPACK_SRC(dtpmlqt  dtpmlqt.f stpmlqt.f ctpmlqt.f ztpmlqt.f)
ADD_LAPACK_SRC(dtpmqrt  dtpmqrt.f stpmqrt.f ctpmqrt.f ztpmqrt.f)
ADD_LAPACK_SRC(dtpqrt   dtpqrt.f ctpqrt.f stpqrt.f ztpqrt.f)
ADD_LAPACK_SRC(dtpqrt2  dtpqrt2.f ctpqrt2.f stpqrt2.f ztpqrt2.f)
ADD_LAPACK_SRC(dtprfb   dtprfb.f stprfb.f ctprfb.f ztprfb.f)
ADD_LAPACK_SRC(dtrevc   dtrevc.f  strevc.f ctrevc.f ztrevc.f)
ADD_LAPACK_SRC(dtrexc   dtrexc.f strexc.f ctrexc.f ztrexc.f)
ADD_LAPACK_SRC(dtrsen   dtrsen.f strsen.f ctrsen.f ztrsen.f)
ADD_LAPACK_SRC(dtrsyl   dtrsyl.f strsyl.f ctrsyl.f ztrsyl.f)
ADD_LAPACK_SRC(dtrti2   dtrti2.f strti2.f ctrti2.f ztrti2.f)
ADD_LAPACK_SRC(dzsum1   dzsum1.f scsum1.f)
ADD_LAPACK_SRC(ieeeck   ieeeck.f )
ADD_LAPACK_SRC(iladlc   iladlc.f ilaslc.f ilaclc.f ilazlc.f)
ADD_LAPACK_SRC(iladlr   iladlr.f ilaslr.f ilaclr.f ilazlr.f)
ADD_LAPACK_SRC(ilaenv   ilaenv.f )
ADD_LAPACK_SRC(iparam2stage iparam2stage.F)
ADD_LAPACK_SRC(iparmq   iparmq.f )
ADD_LAPACK_SRC(lsame    lsame.f lsamen.f)
ADD_LAPACK_SRC(sdsdot   sdsdot.f )
ADD_LAPACK_SRC(zlacgv   zlacgv.f clacgv.f)
ADD_LAPACK_SRC(zung2r   zung2r.f cung2r.f)
ADD_LAPACK_SRC(zunghr   zunghr.f cunghr.f)
ADD_LAPACK_SRC(zungqr   zungqr.f cungqr.f)
ADD_LAPACK_SRC(zunm2r   zunm2r.f cunm2r.f)
ADD_LAPACK_SRC(zunmbr   zunmbr.f cunmbr.f)
ADD_LAPACK_SRC(zunmqr   zunmqr.f cunmqr.f)
ADD_LAPACK_SRC(zunmqr   zunmqr.f cunmqr.f)
ADD_LAPACK_SRC(zunmr2   zunmr2.f cunmr2.f)
ADD_LAPACK_SRC(zunmhr   zunmhr.f cunmhr.f) 
ADD_LAPACK_SRC(izmax1   izmax1.f icmax1.f) 
ADD_LAPACK_SRC(zunmlq   zunmlq.f cunmlq.f) 
ADD_LAPACK_SRC(zunml2   zunml2.f cunml2.f) 


IF(SRC)
	ADD_LIBRARY(${LAPACK_MISSING_TARGET_NAME} OBJECT ${SRC})
    SET_TARGET_PROPERTIES(${LAPACK_MISSING_TARGET_NAME} PROPERTIES POSITION_INDEPENDENT_CODE TRUE)
    IF( CMAKE_Fortran_COMPILER_ID STREQUAL "XL")
		SET_TARGET_PROPERTIES(${LAPACK_MISSING_TARGET_NAME} PROPERTIES COMPILE_FLAGS "${CMAKE_Fortran_FLAGS} -qfixed -qnosave")
ENDIF()


ENDIF()
UNSET(SRC)
