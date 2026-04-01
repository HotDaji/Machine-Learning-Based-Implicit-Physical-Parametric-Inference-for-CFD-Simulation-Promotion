cmake_minimum_required(VERSION 3.20)  # 提高最低版本，3.2 可能过低
project(main)

# 查找源文件
aux_source_directory(. DIR_SRCS)

# 设置编译器
set(CMAKE_CXX_COMPILER "mpicxx")  # 使用 MPI 编译器

# 设置编译选项
set(CMAKE_CXX_STANDARD 17)        # 更新为 C++14 或更高
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -O3")  # 添加优化标志

# 查找 libtorch
find_package(Torch REQUIRED PATHS "/public/home/chaosuan/SHU_24820214/lib/libtorch/share/cmake/Torch")

# 如果找不到 Torch，则给出错误提示
if (NOT Torch_FOUND)
    message(FATAL_ERROR "libtorch not found. Ensure that the path is correct.")
endif()

# 创建可执行文件
add_executable(main ${DIR_SRCS})

# 链接 libtorch 库
target_link_libraries(main "${TORCH_LIBRARIES}")

# 添加 RPATH，确保运行时能找到 libtorch 的共享库
set_target_properties(main PROPERTIES BUILD_RPATH "/public/home/chaosuan/SHU_24820214/lib/libtorch/lib")
