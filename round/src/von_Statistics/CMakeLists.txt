cmake_minimum_required(VERSION 3.24)
project(von_dpmm LANGUAGES CXX CUDA)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CUDA_STANDARD 17)
set(CMAKE_CUDA_STANDARD_REQUIRED ON)

if(NOT DEFINED CMAKE_CUDA_ARCHITECTURES)
  set(CMAKE_CUDA_ARCHITECTURES 75)
endif()

set(CMAKE_RUNTIME_OUTPUT_DIRECTORY            "${CMAKE_SOURCE_DIR}/bin")
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY_RELEASE    "${CMAKE_SOURCE_DIR}/bin")
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY_DEBUG      "${CMAKE_SOURCE_DIR}/bin")

add_executable(von_dpmm
  src/main.cpp
  src/config.cpp
  src/io.cpp
  src/cpu.cpp
  src/gpu.cu
)

target_include_directories(von_dpmm PRIVATE src)

# 中文代码页(936)下 cuda/msvc 头会刷 C4819,统一关掉;源文件按 UTF-8 读
target_compile_options(von_dpmm PRIVATE
  $<$<COMPILE_LANGUAGE:CXX>:/utf-8 /wd4819 /W3 /O2>
  $<$<COMPILE_LANGUAGE:CUDA>:-Xcompiler=/utf-8 -Xcompiler=/wd4819 -Xcompiler=/O2>
)

find_package(Threads REQUIRED)
target_link_libraries(von_dpmm PRIVATE Threads::Threads)
