# 所有源文件形成变量
aux_source_directory(source PID_DIRS)
# 创建静态库并添加源文件
add_library(PidCpp STATIC ${PID_DIRS})
# 设置静态库搜索路径
target_include_directories(PidCpp PUBLIC include ../../Datastruct)
# 指定在链接该库时需要使用的库
# target_link_libraries(PidCpp)