# 设置 Eigen 库的头文件路径
set(EIGEN3_INCLUDE_DIR "D:/Libraries/eigen-3.3.4/eigen-3.3.4")
# set(qpOASES_INCLUDE_DIR "D:/Libraries/qpOASES-stable-3.2/qpOASES-stable-3.2/include")

# 所有源文件形成变量
aux_source_directory(source MPC_DIRS)

# 创建静态库并添加源文件
add_library(MpcCpp STATIC ${MPC_DIRS})

# 设置包含路径
# include_directories(${EIGEN3_INCLUDE_DIR} libs ${qpOASES_INCLUDE_DIR})

# 设置静态库搜索路径
# target_include_directories(MpcCpp PUBLIC include ${EIGEN3_INCLUDE_DIR} ${qpOASES_INCLUDE_DIR})
target_include_directories(MpcCpp PUBLIC include ${EIGEN3_INCLUDE_DIR})

# 链接qpOASES库
# link_directories("${CMAKE_CURRENT_SOURCE_DIR}/libs")

# 查找库文件
# find_library(QPOASES_LIB qpOASES PATHS "${CMAKE_CURRENT_SOURCE_DIR}/libs")

# 指定在链接该库时需要使用的库
# target_link_libraries(MpcCpp ${QPOASES_LIB})
target_link_libraries(MpcCpp "${CMAKE_CURRENT_SOURCE_DIR}/libs/libqpOASES.a")