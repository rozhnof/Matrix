cmake_minimum_required(VERSION 3.16)

project(matrix CXX)

include(cmake/CompileOptions.cmake)
include(cmake/Sanitize.cmake)
include(cmake/Dependencies.cmake)

set(CMAKE_INCLUDE_CURRENT_DIR ON)
include_directories(${CMAKE_SOURCE_DIR}/include)

file(GLOB_RECURSE HEADER_FILE_LIST CONFIGURE_DEPENDS
    "include/*.h"
)

file(GLOB_RECURSE TEST_FILES CONFIGURE_DEPENDS 
    "tests/*.cc"
)

add_library(${PROJECT_NAME} INTERFACE)
target_sources(${PROJECT_NAME} INTERFACE ${HEADER_FILE_LIST})


if(BUILD_TESTS)
    add_executable(${PROJECT_NAME}-test ${TEST_FILES})
    target_link_libraries(${PROJECT_NAME}-test PUBLIC ${PROJECT_NAME} gtest gtest_main)
endif()
