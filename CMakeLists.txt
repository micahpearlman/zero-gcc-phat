cmake_minimum_required(VERSION 3.0.0)
project(zero-gcc-phat VERSION 0.1.0)
find_package(PkgConfig) # support for Makefile style pkg-config

### Log system info
message("-- system arch:  ${CMAKE_SYSTEM_PROCESSOR}")

### Get the required packages

# FFTW
include(cmake/findFFTW.cmake)
find_package(FFTW REQUIRED)
message(STATUS "Found FFTW ${FFTW_FOUND}")
message(STATUS "    libraries: ${FFTW_LIBRARIES}")


### C++ config
set(CMAKE_CXX_STANDARD 11)

### Create Zero GCC-PHAT library
add_library(zo-gcc-phat src/zo-gcc-phat.cpp)
target_include_directories(zo-gcc-phat PRIVATE ${FFTW_INCLUDES})
target_link_libraries(zo-gcc-phat ${FFTW_LIBRARIES})


### Create Test App
add_executable(zo-gcc-phat-test test/zo-gcc-phat-test.cpp)
target_include_directories(zo-gcc-phat-test PRIVATE ./src)

target_link_libraries(zo-gcc-phat-test zo-gcc-phat)
target_link_libraries(zo-gcc-phat-test ${FFTW_LIBRARIES})
 

