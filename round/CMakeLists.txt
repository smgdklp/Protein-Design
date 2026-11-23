cmake_minimum_required(VERSION 3.25)
project(Round LANGUAGES CXX)

# ---------------------------------------------------------------------------
# 工作目录主 CMake：只负责全局设置、统一产物目录、挂载各子项目。
# 每个子项目在自己的目录里有独立 CMakeLists.txt。
# ---------------------------------------------------------------------------

set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

if(MSVC)
    add_compile_options(/utf-8 /W4 /permissive- /Zc:__cplusplus)
    add_compile_definitions(_CRT_SECURE_NO_WARNINGS NOMINMAX WIN32_LEAN_AND_MEAN UNICODE _UNICODE)
endif()

# 统一编译产物目录：<源码根>/bin/<Config>/
set(ROUND_OUTPUT_ROOT "${CMAKE_SOURCE_DIR}/bin" CACHE PATH "Unified runtime output root")
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY "${ROUND_OUTPUT_ROOT}/$<CONFIG>")

# 子项目
add_subdirectory(src/date)
add_subdirectory(src/process_Statistics)
add_subdirectory(src/peak_Statistics)
