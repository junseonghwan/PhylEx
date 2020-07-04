# code relating to tree structured stick breaking process
include_directories(/usr/local/include/ ${PROJECT_SOURCE_DIR}/src ${PROJECT_SOURCE_DIR}/eigen/)

set(Boost_USE_STATIC_LIBS ON)

find_package (PkgConfig REQUIRED)
find_package(Boost 1.61 COMPONENTS unit_test_framework REQUIRED)
find_package (GSL REQUIRED)
if(Boost_FOUND)
    include_directories(${Boost_INCLUDE_DIRS})
endif()

add_executable(tests main_test.cpp) # create the executable file for testing
# excluded tests: 
# test_dynamic_program.cpp

target_link_libraries(tests GSL::gsl GSL::gslcblas ${Boost_LIBRARIES} tssb)

add_test(Tester tests) # add test
