SET(UTILS_SRC utils.f90)

ADD_LIBRARY(utils STATIC ${UTILS_SRC})
target_link_libraries(utils PUBLIC qrupdate)
ADD_COVERAGE(utils)


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
    FIND_PROGRAM(WINE_CMD wine
        HINTS
          /usr/bin
          /usr/local/bin
          ENV WINEROOT
    )
    if (NOT WINE_CMD)
        MESSAGE(FATAL_ERROR "Wine not found. Consider setting WINEROOT")
    endif()
    set(CMAKE_CROSSCOMPILING_EMULATOR "${WINE_CMD}")
    SET(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} ${Fortran_LOCAL_FLAGS}")
    # FIND_PROGRAM(WINE NAMES wine wine32 wine64)
    # IF ( NOT WINE )
        # MESSAGE(FATAL_ERROR "Wine not found")
    # ELSE()
        # MESSAGE(STATUS "Wine found: ${WINE}")
    # ENDIF()

    EXECUTE_PROCESS(COMMAND ${CMAKE_C_COMPILER} -print-search-dirs
        COMMAND grep programs
        COMMAND cut -d= -f2
        COMMAND sed -e "s/:/;/g"
        COMMAND tr -d "\n" OUTPUT_VARIABLE WINEPATH_TMP)
    SET(WINEPATH)
    SET(ENVPATH_TMP "$ENV{WINEPATH}")
    FOREACH ( P IN LISTS ENVPATH_TMP)
        execute_process(COMMAND winepath -0 -w ${P} OUTPUT_VARIABE PX ERROR_QUIET)
        # message(STATUS "PX: ${PX}")
        SET(WINEPATH "${WINEPATH}\\;${PX}")
    ENDFOREACH()
    FOREACH ( P IN LISTS CMAKE_FIND_ROOT_PATH )
        execute_process(COMMAND winepath -0 -w ${P}/bin OUTPUT_VARIABLE PX ERROR_QUIET)
        # message(STATUS "PX: ${PX}")
        SET(WINEPATH "${WINEPATH}\\;${PX}")
        execute_process(COMMAND winepath -0 -w ${P}/lib OUTPUT_VARIABLE PX ERROR_QUIET)
        # message(STATUS "PX: ${PX}")
        SET(WINEPATH "${WINEPATH}\\;${PX}")
    ENDFOREACH()


    FOREACH ( P IN LISTS WINEPATH_TMP )
        execute_process(COMMAND winepath -0 -w ${P} OUTPUT_VARIABLE PX ERROR_QUIET)
        # message(STATUS "PX: ${PX}")
        SET(WINEPATH "${WINEPATH}\\;${PX}")
    ENDFOREACH()

    execute_process(COMMAND winepath -0 -w ${CMAKE_BINARY_DIR}/src OUTPUT_VARIABLE PX ERROR_QUIET)
    SET(WINEPATH "${PX}\\;${WINEPATH}")

    MESSAGE(STATUS "Setting WINEPATH for tests to ${WINEPATH}")
ENDIF()

IF (CMAKE_SYSTEM_NAME MATCHES "CYGWIN")
    MESSAGE(STATUS "CYGWIN build")
    MESSAGE("Set additional PATH for tests to ${CMAKE_BINARY_DIR}/src:/usr/lib/lapack")
    SET(ADDPATH  "${CMAKE_BINARY_DIR}/src:/usr/lib/lapack")
ENDIF()


FOREACH(TEST_CODE ${TEST_PRGM})
    ADD_EXECUTABLE(${TEST_CODE} ${TEST_CODE}.f90)
    ADD_COVERAGE(${TEST_CODE})
    TARGET_LINK_LIBRARIES(${TEST_CODE} utils qrupdate )
    # set_target_properties(${TEST_CODE} PROPERTIES CROSSCOMPILING_EMULATOR "/usr/bin/wine")
    IF ( CMAKE_CROSSCOMPILING )
        IF ( CMAKE_SYSTEM_NAME MATCHES Windows)
            ADD_TEST(NAME test_${TEST_CODE} COMMAND ${TEST_CODE})
            SET_TESTS_PROPERTIES(test_${TEST_CODE} PROPERTIES ENVIRONMENT "WINEPATH=${WINEPATH}")
            # SET_TESTS_PROPERTIES(test_${TEST_CODE} PROPERTIES CROSSCOMPILING_EMULATOR "/usr/bin/wine")
        ELSE()
            MESSAGE(WARNING "Cross-Compile checks are only possible from Linux to Windows.")
        ENDIF()
    ELSE ()
        ADD_TEST(NAME test_${TEST_CODE} COMMAND ${TEST_CODE})
        IF ( ADDPATH )
            SET_TESTS_PROPERTIES(test_${TEST_CODE} PROPERTIES ENVIRONMENT "PATH=${ADDPATH}:$ENV{PATH}")
        ENDIF()
    ENDIF()
ENDFOREACH()

