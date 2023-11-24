# define the project name
project(PhylEx)

cmake_minimum_required (VERSION 3.12)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_C_FLAGS "-Wall -g" )

#set(Boost_USE_STATIC_LIBS ON)
#set(Boost_USE_MULTITHREADED OFF)

add_compile_definitions(_GLIBCXX_USE_CXX11_ABI=0)

if (APPLE)
    set(OPENMP_CXX_INCLUDE_DIRS /usr/local/opt/libomp/include)
    set(OpenMP_C_FLAGS "-Xpreprocessor -fopenmp -I/usr/local/opt/libomp/include" CACHE INTERNAL "OpenMP flags for Xcode toolchain.")
    set(OpenMP_C_LIB_NAMES "omp" CACHE INTERNAL "OpenMP lib name for Xcode toolchain.")
    set(OpenMP_CXX_FLAGS "-Xpreprocessor -fopenmp -I/usr/local/opt/libomp/include" CACHE INTERNAL "OpenMP flags for Xcode toolchain.")
    set(OpenMP_CXX_LIB_NAMES "omp" CACHE INTERNAL "OpenMP lib name for Xcode toolchain.")
    set(OpenMP_omp_LIBRARY "/usr/local/opt/libomp/lib/libomp.dylib" CACHE INTERNAL "OpenMP lib name for Xcode toolchain.")
endif()

enable_testing()

add_subdirectory(src)
add_subdirectory(simulation)
add_subdirectory(tests)
