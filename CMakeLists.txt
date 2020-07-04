# define the project name
project(sc-clustering)

cmake_minimum_required (VERSION 3.12)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set( CMAKE_C_FLAGS "-Wall -g" )

set(Boost_USE_STATIC_LIBS ON)
set(Boost_USE_MULTITHREADED OFF)

enable_testing()

add_subdirectory(src)
add_subdirectory(simulation)
add_subdirectory(tests)
