# src/date —— PDB 片段提取器（第一个子项目）
# 只负责编译本项目的可执行文件；全局设置与产物目录由根 CMakeLists 控制。

add_executable(pdbextract
    main.cpp
    config.cpp
    pdb.cpp
    segment.cpp
    util.cpp
)

target_include_directories(pdbextract PRIVATE "${CMAKE_CURRENT_SOURCE_DIR}")
target_compile_features(pdbextract PRIVATE cxx_std_20)

# 运行期配置随 exe 一起落到统一产物目录，保证"exe 读同目录 pbd_config.json"
add_custom_command(TARGET pdbextract POST_BUILD
    COMMAND ${CMAKE_COMMAND} -E copy_if_different
            "${CMAKE_CURRENT_SOURCE_DIR}/pbd_config.json"
            "$<TARGET_FILE_DIR:pdbextract>/pbd_config.json"
    COMMENT "Copy pbd_config.json next to pdbextract.exe"
)
