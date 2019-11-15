CMAKE_MINIMUM_REQUIRED(VERSION 2.8.12)

PROJECT(deal_II_tests)

# load cmake modules (for things like where to find libraries)
set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "${CMAKE_SOURCE_DIR}/cmake/modules/")

# set various compiler flags
add_compile_options(-std=c++11)
set(CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG} -O0 -Wall -Wextra")
set(CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS_RELEASE} -O3 -Wall -Wextra")
set(CMAKE_VERBOSE_MAKEFILE ON)

FIND_PACKAGE(deal.II 8.1.0 REQUIRED HINTS 
  ${DEAL_II_DIR} ../ ../../ $ENV{DEAL_II_DIR})

DEAL_II_INITIALIZE_CACHED_VARIABLES()

add_subdirectory(src)
