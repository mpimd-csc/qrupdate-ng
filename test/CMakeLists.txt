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

FOREACH(TEST_CODE ${TEST_PRGM})
    ADD_EXECUTABLE(${TEST_CODE} ${TEST_CODE}.f)
    TARGET_LINK_LIBRARIES(${TEST_CODE} utils qrupdate)
    ADD_TEST(NAME test_${TEST_CODE} COMMAND ${CMAKE_CURRENT_BINARY_DIR}/${TEST_CODE})
ENDFOREACH()

