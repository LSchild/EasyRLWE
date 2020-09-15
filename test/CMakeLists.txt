if(DEBUG)
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -g -O0")
    add_compile_definitions(DEBUG)
else()
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -O3")
endif()

set(TEST_SOURCES "${CMAKE_CURRENT_SOURCE_DIR}/test/main_test.cpp")

add_executable(lwe_structs "${TEST_SOURCES}")
target_include_directories(lwe_structs PRIVATE "${CMAKE_CURRENT_SOURCE_DIR}/include")