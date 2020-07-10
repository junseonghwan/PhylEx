# code relating to tree structured stick breaking process
include_directories(/usr/local/include/)

file( GLOB SRCS *.cpp *.hpp)
add_library(tssb STATIC ${SRCS}) # create static library
add_executable(run ${SRCS}) # create the executable file

#find_package (OpenMP REQUIRED)
find_package (PkgConfig REQUIRED)
find_package(Boost 1.61 COMPONENTS regex system filesystem program_options REQUIRED)
find_package (GSL REQUIRED)
#find_package(Eigen3 REQUIRED)
if(Boost_FOUND)
    include_directories(${Boost_INCLUDE_DIRS})
endif()

target_link_libraries(tssb GSL::gsl GSL::gslcblas ${Boost_LIBRARIES})
target_link_libraries(run GSL::gsl GSL::gslcblas ${Boost_LIBRARIES})

install(TARGETS tssb DESTINATION bin)
INSTALL(TARGETS run DESTINATION ${CMAKE_SOURCE_DIR}/)
