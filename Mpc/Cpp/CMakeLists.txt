# 设置头文件路径
set(INCLUDE_DIR "${CMAKE_CURRENT_SOURCE_DIR}/config"
                "${CMAKE_CURRENT_SOURCE_DIR}/include" 
                "${CMAKE_CURRENT_SOURCE_DIR}/include/osqp100/osqp" 
                "${CMAKE_CURRENT_SOURCE_DIR}/include/OsqpEigen100"
                "${CMAKE_CURRENT_SOURCE_DIR}/include/qpOASES"
                "${CMAKE_CURRENT_SOURCE_DIR}/include/QuadProg"
                "${CMAKE_CURRENT_SOURCE_DIR}/include/TinyMPC")

# 所有源文件形成变量
aux_source_directory(source MPC_interface)
aux_source_directory(source/qpOASES MPC_qpOASES)
aux_source_directory(source/QuadProg MPC_quadprog)
aux_source_directory(source/TinyMPC MPC_tinyMPC)
aux_source_directory(source/OsqpEigen MPC_OSQPE)

# 创建静态库并添加源文件
add_library(MpcCpp STATIC ${MPC_interface} ${MPC_qpOASES} ${MPC_quadprog} ${MPC_tinyMPC} ${MPC_OSQPE})

# 设置包含路径
# include_directories(${EIGEN3_INCLUDE_DIR} libs ${qpOASES_INCLUDE_DIR})

# 设置目标包含路径
target_include_directories(MpcCpp PUBLIC include ${INCLUDE_DIR})

# 链接qpOASES库
# link_directories("${CMAKE_CURRENT_SOURCE_DIR}/libs")

# 查找库文件
# find_library(QPOASES_LIB qpOASES PATHS "${CMAKE_CURRENT_SOURCE_DIR}/libs")

# 指定在链接该库时需要使用的库
if(CMAKE_SYSTEM_NAME STREQUAL "Windows")
    target_link_libraries(MpcCpp "${CMAKE_CURRENT_SOURCE_DIR}/libs/osqp100/osqpstatic.lib")
elseif(CMAKE_SYSTEM_NAME STREQUAL "Linux")
    target_link_libraries(MpcCpp "${CMAKE_CURRENT_SOURCE_DIR}/bin/osqp100/libosqp.so")
    # target_link_libraries(MpcCpp "${CMAKE_CURRENT_SOURCE_DIR}/bin/osqp100/libosqpstatic.a")
endif()