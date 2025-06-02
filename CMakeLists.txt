cmake_minimum_required(VERSION 3.5)
project(Test)
set(CMAKE_CXX_STANDARD 14)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

# 设置默认优化等级（影响所有构建类型）
# set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -O2")
# set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -O2")
# 设置默认优化标志（影响所有编译器）
set(CMAKE_BUILD_TYPE "Release")  # 强制Release，会覆盖vscode的配置

# add_compile_options(-Wno-psabi)

add_subdirectory(Mpc/Cpp)
# add_subdirectory(Pid/Cpp)

add_executable(Testfile test.cpp)
# target_link_libraries(Testfile MpcCpp PidCpp)
target_link_libraries(Testfile MpcCpp)