# code relating to tree structured stick breaking process
include_directories(/usr/local/include/ ${PROJECT_SOURCE_DIR}/src ${PROJECT_SOURCE_DIR}/eigen/)

set(Boost_USE_STATIC_LIBS ON)

find_package (PkgConfig REQUIRED)
find_package(Boost 1.61 COMPONENTS unit_test_framework filesystem REQUIRED)
find_package (GSL REQUIRED)
if(Boost_FOUND)
    include_directories(${Boost_INCLUDE_DIRS})
endif()

add_executable(testing main_test.cpp) # create the executable file for testing

target_link_libraries(testing GSL::gsl GSL::gslcblas ${Boost_LIBRARIES} tssb)

add_test(Tester testing) # add test
INSTALL(TARGETS testing DESTINATION ${CMAKE_SOURCE_DIR}/)