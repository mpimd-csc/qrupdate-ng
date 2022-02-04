SET(UTILS_SRC utils.f)

ADD_LIBRARY(utils STATIC ${UTILS_SRC})

SET(TEST_PRGM
    tch1dn
    tch1up
    tchdex
    tchinx
    tchshx
    tlu1up
    tlup1up
    tqr1up
    tqrdec
    tqrder
    tqrinc
    tqrinr
    tqrshc
    )


IF ( CMAKE_CROSSCOMPILING AND CMAKE_SYSTEM_NAME MATCHES Windows )
    SET(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} ${Fortran_LOCAL_FLAGS}")
    FIND_PROGRAM(WINE NAMES wine wine32 wine64)
    IF ( NOT WINE )
        MESSAGE(FATAL_ERROR "Wine not found")
    ELSE()
        MESSAGE(STATUS "Wine found: ${WINE}")
    ENDIF()

    EXECUTE_PROCESS(COMMAND ${CMAKE_C_COMPILER} -print-search-dirs
        COMMAND grep programs
        COMMAND cut -d= -f2
        COMMAND sed -e "s/:/;/g"
        COMMAND tr -d "\n" OUTPUT_VARIABLE WINEPATH_TMP)
    SET(WINEPATH)
    SET(ENVPATH_TMP "$ENV{WINEPATH}")
    FOREACH ( P IN LISTS ENVPATH_TMP)
        SET(WINEPATH "${WINEPATH}\\;${P}")
    ENDFOREACH()
    FOREACH ( P IN LISTS CMAKE_FIND_ROOT_PATH )
        SET(WINEPATH "${WINEPATH}\\;${P}/bin")
        SET(WINEPATH "${WINEPATH}\\;${P}/lib")
    ENDFOREACH()


    FOREACH ( P IN LISTS WINEPATH_TMP )
        SET(WINEPATH "${WINEPATH}\\;${P}")
    ENDFOREACH()

    SET(WINEPATH "${CMAKE_BINARY_DIR}/src\\;${WINEPATH}")
    MESSAGE(STATUS "Setting WINEPATH for tests to ${WINEPATH}")
ENDIF()

IF (CMAKE_SYSTEM_NAME MATCHES "CYGWIN")
    MESSAGE(STATUS "CYGWIN build")
    MESSAGE("Set additional PATH for tests to ${CMAKE_BINARY_DIR}/src:/usr/lib/lapack")
    SET(ADDPATH  "${CMAKE_BINARY_DIR}/src:/usr/lib/lapack")
ENDIF()


FOREACH(TEST_CODE ${TEST_PRGM})
    ADD_EXECUTABLE(${TEST_CODE} ${TEST_CODE}.f)
    TARGET_LINK_LIBRARIES(${TEST_CODE} utils qrupdate )
    IF ( CMAKE_CROSSCOMPILING )
        IF ( CMAKE_SYSTEM_NAME MATCHES Windows)
            ADD_TEST(NAME test_${TEST_CODE} COMMAND ${WINE} ${CMAKE_CURRENT_BINARY_DIR}/${TEST_CODE}${CMAKE_EXECUTABLE_SUFFIX})
            SET_TESTS_PROPERTIES(test_${TEST_CODE} PROPERTIES ENVIRONMENT "WINEPATH=${WINEPATH}")
        ELSE()
            MESSAGE(WARNING "Cross-Compile checks are only possible from Linux to Windows.")
        ENDIF()
    ELSE ()
        ADD_TEST(NAME test_${TEST_CODE} COMMAND ${CMAKE_CURRENT_BINARY_DIR}/${TEST_CODE})
        IF ( ADDPATH )
            SET_TESTS_PROPERTIES(test_${TEST_CODE} PROPERTIES ENVIRONMENT "PATH=${ADDPATH}:$ENV{PATH}")
        ENDIF()
    ENDIF()
ENDFOREACH()

