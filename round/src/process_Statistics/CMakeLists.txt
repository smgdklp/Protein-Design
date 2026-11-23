# src/process_Statistics —— 二面角条件概率统计（20x20x20 张 72x72 矩阵）
# 只负责编译本项目的可执行文件；全局设置与产物目录由根 CMakeLists 控制。

add_executable(statgen
    main.cpp
    config.cpp
    npy.cpp
)

target_include_directories(statgen PRIVATE "${CMAKE_CURRENT_SOURCE_DIR}")
target_compile_features(statgen PRIVATE cxx_std_20)

# 运行期配置随 exe 一起落到统一产物目录，保证"exe 读同目录 config.json"
add_custom_command(TARGET statgen POST_BUILD
    COMMAND ${CMAKE_COMMAND} -E copy_if_different
            "${CMAKE_CURRENT_SOURCE_DIR}/config.json"
            "$<TARGET_FILE_DIR:statgen>/config.json"
    COMMENT "Copy config.json next to statgen.exe"
)
