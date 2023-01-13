cmake_minimum_required(VERSION 3.0.0)
project(medial_constraint_icp VERSION 0.1.0)

include(CTest)
enable_testing()

add_executable(medial_constraint_icp example_usage.cpp)

set(CPACK_PROJECT_NAME ${PROJECT_NAME})
set(CPACK_PROJECT_VERSION ${PROJECT_VERSION})
include(CPack)
