cmake_minimum_required(VERSION 3.16)
project(thermostat)
set(CMAKE_CXX_STANDARD 20)


if (MSVC)
add_compile_options(-bigobj)
if (CMAKE_CXX_COMPILER_VERSION VERSION_GREATER_EQUAL 19.12.25835)
    set(CMAKE_CXX20_STANDARD_COMPILE_OPTION "-std:c++latest")
    set(CMAKE_CXX20_EXTENSION_COMPILE_OPTION "-std:c++latest")
endif()
endif()


add_subdirectory("externals/teqp")
add_executable(THERMOSTAT "src/thermostat.cxx")

if(MSVC)
target_sources(THERMOSTAT PUBLIC "${CMAKE_CURRENT_SOURCE_DIR}/externals/teqp/externals/Eigen/debug/msvc/eigen.natvis")
endif()

include_directories(THERMOSTAT "${CMAKE_CURRENT_SOURCE_DIR}/include")

target_link_libraries(THERMOSTAT PRIVATE autodiff PRIVATE teqpinterface PRIVATE Catch2WithMain PRIVATE teqpcpp)

target_compile_definitions(THERMOSTAT PRIVATE -DTEQP_MULTICOMPLEX_ENABLED)