# code relating to tree structured stick breaking process
include_directories(/usr/local/include/)

file( GLOB SRCS *.cpp *.hpp)
add_library(tssb STATIC ${SRCS}) # create static library
add_executable(run ${SRCS}) # create the executable file

find_package (PkgConfig REQUIRED)
find_package(Boost 1.61 COMPONENTS regex system filesystem program_options REQUIRED)
find_package (GSL REQUIRED)
if(Boost_FOUND)
    include_directories(${Boost_INCLUDE_DIRS})
endif()
find_package(OpenMP REQUIRED)
if (OPENMP_FOUND)
    set (CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS}")
    set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
    set (CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} ${OpenMP_EXE_LINKER_FLAGS}")
endif()

target_link_libraries(tssb GSL::gsl GSL::gslcblas ${Boost_LIBRARIES} OpenMP::OpenMP_CXX)
target_link_libraries(run GSL::gsl GSL::gslcblas ${Boost_LIBRARIES} OpenMP::OpenMP_CXX)

install(TARGETS tssb DESTINATION bin)
INSTALL(TARGETS run DESTINATION ${CMAKE_SOURCE_DIR}/)
