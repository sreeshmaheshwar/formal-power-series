cmake_minimum_required(VERSION 3.17)
project(FormalPowerSeries)

set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

add_compile_options(-Wall -Wextra -Wpedantic)
add_compile_options(-I${CMAKE_SOURCE_DIR}/../ac-library)

enable_testing()

add_subdirectory(googletest)
include_directories(${gtest_SOURCE_DIR}/include ${gtest_SOURCE_DIR})
include(GoogleTest)

include_directories(.)
include_directories(${CMAKE_SOURCE_DIR}/..)

add_executable(FormalPowerSeriesTest FormalPowerSeriesTest.cpp)
target_link_libraries(FormalPowerSeriesTest gtest gtest_main)
gtest_discover_tests(FormalPowerSeriesTest)
