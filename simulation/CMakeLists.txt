# code relating to tree structured stick breaking process
include_directories(/usr/local/include/ ${PROJECT_SOURCE_DIR}/src ${PROJECT_SOURCE_DIR}/eigen/)

find_package (PkgConfig REQUIRED)
find_package(Boost COMPONENTS regex system filesystem program_options REQUIRED)
find_package (GSL REQUIRED)
#find_package(Eigen3 REQUIRED)
if(Boost_FOUND)
    include_directories(${Boost_INCLUDE_DIRS})
endif()

file( GLOB SRCS *.cpp *.hpp)

add_executable(simul ${SRCS}) # create the executable file for testing
target_link_libraries(simul GSL::gsl GSL::gslcblas ${Boost_LIBRARIES} tssb)

INSTALL(TARGETS simul DESTINATION ${CMAKE_SOURCE_DIR}/)
