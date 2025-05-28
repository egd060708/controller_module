cmake_minimum_required(VERSION 3.5)
project(Test)
set(CMAKE_CXX_STANDARD 14)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

add_subdirectory(Mpc/Cpp)
# add_subdirectory(Pid/Cpp)

add_executable(Testfile test.cpp)
# target_link_libraries(Testfile MpcCpp PidCpp)
target_link_libraries(Testfile MpcCpp)