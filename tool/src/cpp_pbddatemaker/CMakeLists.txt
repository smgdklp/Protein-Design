cmake_minimum_required(VERSION 3.25)

# .pdb / .pdb.gz  ->  定长二面角二进制(.bin)
#   struct List { int _len; (char _kind, float32 ca, sa, cb, sb)[N] }  // 每单元 17 字节
# 配置在工程根目录 config.json, 构建时复制到编译产物同目录, 运行时由主程序从该目录读取。
project(pbddatemaker LANGUAGES C CXX)

# ---------------------------------------------------------------- 全局设置
set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)
set(CMAKE_C_STANDARD 11)

# 静态 CRT: 产物不依赖 VC 运行库分发
set(CMAKE_MSVC_RUNTIME_LIBRARY "MultiThreaded$<$<CONFIG:Debug>:Debug>")

if(NOT CMAKE_CONFIGURATION_TYPES AND NOT CMAKE_BUILD_TYPE)
    set(CMAKE_BUILD_TYPE Release CACHE STRING "构建类型" FORCE)
endif()

# ---------------------------------------------------------------- 产物目录
# 可执行文件与 config.json 同目录: tool/result/pbd_datamaker
get_filename_component(PBD_OUT_DIR "${CMAKE_CURRENT_SOURCE_DIR}/../../result/pbd_datamaker" ABSOLUTE)
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY                 "${PBD_OUT_DIR}")
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY_DEBUG           "${PBD_OUT_DIR}")
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY_RELEASE         "${PBD_OUT_DIR}")
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY_RELWITHDEBINFO  "${PBD_OUT_DIR}")
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY_MINSIZEREL      "${PBD_OUT_DIR}")
set(CMAKE_PDB_OUTPUT_DIRECTORY "${CMAKE_CURRENT_SOURCE_DIR}/build/symbols")

# ---------------------------------------------------------------- zlib 1.3.1
# 只取 gz 读取所需的最小源码集(inflate + gz 层), 静态链接
set(ZLIB_SRC_DIR "D:/Cpp/Packages/zlib-1.3.1" CACHE PATH "zlib 1.3.1 源码目录")
if(NOT EXISTS "${ZLIB_SRC_DIR}/zlib.h")
    message(FATAL_ERROR "找不到 zlib 源码目录: ${ZLIB_SRC_DIR}")
endif()

add_library(pbd_zlib STATIC
    "${ZLIB_SRC_DIR}/adler32.c"
    "${ZLIB_SRC_DIR}/crc32.c"
    "${ZLIB_SRC_DIR}/deflate.c"
    "${ZLIB_SRC_DIR}/gzclose.c"
    "${ZLIB_SRC_DIR}/gzlib.c"
    "${ZLIB_SRC_DIR}/gzread.c"
    "${ZLIB_SRC_DIR}/gzwrite.c"
    "${ZLIB_SRC_DIR}/inffast.c"
    "${ZLIB_SRC_DIR}/inflate.c"
    "${ZLIB_SRC_DIR}/inftrees.c"
    "${ZLIB_SRC_DIR}/trees.c"
    "${ZLIB_SRC_DIR}/zutil.c")
target_include_directories(pbd_zlib PUBLIC "${ZLIB_SRC_DIR}")
if(MSVC)
    target_compile_options(pbd_zlib PRIVATE /W0)
endif()

# ---------------------------------------------------------------- 可执行文件
add_executable(pbddatemaker
    src/main.cpp
    src/config.cpp
    src/reader.cpp)

target_include_directories(pbddatemaker PRIVATE "${CMAKE_CURRENT_SOURCE_DIR}/src")
target_link_libraries(pbddatemaker PRIVATE pbd_zlib)

if(MSVC)
    target_compile_options(pbddatemaker PRIVATE /utf-8 /W4 /permissive- /EHsc /Zc:__cplusplus)
    target_compile_definitions(pbddatemaker PRIVATE _CRT_SECURE_NO_WARNINGS NOMINMAX WIN32_LEAN_AND_MEAN)
endif()

# 配置必须躺在产物旁边(主程序按 exe 所在目录找 config.json)
add_custom_command(TARGET pbddatemaker POST_BUILD
    COMMAND ${CMAKE_COMMAND} -E make_directory "$<TARGET_FILE_DIR:pbddatemaker>"
    COMMAND ${CMAKE_COMMAND} -E copy_if_different
            "${CMAKE_CURRENT_SOURCE_DIR}/config.json"
            "$<TARGET_FILE_DIR:pbddatemaker>/config.json"
    COMMENT "复制 config.json 到编译产物目录")

set_property(DIRECTORY "${CMAKE_CURRENT_SOURCE_DIR}" PROPERTY VS_STARTUP_PROJECT pbddatemaker)
