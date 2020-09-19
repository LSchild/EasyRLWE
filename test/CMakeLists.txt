if(DEBUG)
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -g -O0 -Wall")
    add_compile_definitions(DEBUG)
else()
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -O3")
endif()

set(TEST_SOURCES main_test.cpp common_test.cpp sample_test.cpp )
set(TEST_HEADERS common_test.h parameters.h sample_test.h)

list(TRANSFORM TEST_SOURCES PREPEND ${CMAKE_CURRENT_SOURCE_DIR}/test/src/)
list(TRANSFORM TEST_HEADERS PREPEND ${CMAKE_CURRENT_SOURCE_DIR}/test/include/)

add_executable(lwe_structs ${TEST_SOURCES} ${TEST_HEADERS})
target_include_directories(lwe_structs PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/include)
target_include_directories(lwe_structs PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/test/include)