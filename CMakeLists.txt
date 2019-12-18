cmake_minimum_required(VERSION 2.8)

project(ISTAT)

set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11 -pthread -O3 -funroll-loops -g")

set(CMAKE_INSTALL_PREFIX ${PROJECT_SOURCE_DIR})

set(ISTAT_INSTALL_INCLUDE_DIR ${PROJECT_SOURCE_DIR}/include)
set(ISTAT_INSTALL_BIN_DIR ${PROJECT_SOURCE_DIR}/bin)
set(ISTAT_INSTALL_LIB_DIR ${PROJECT_SOURCE_DIR}/lib)

set(ISTAT_HEADERS_DIR ${PROJECT_SOURCE_DIR}/src/lib)

include_directories(${ISTAT_INSTALL_INCLUDE_DIR})
include_directories(${ISTAT_HEADERS_DIR})

add_subdirectory(src)

# To be added
#add_subdirectory(test)


