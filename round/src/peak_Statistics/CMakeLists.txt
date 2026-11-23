# src/peak_Statistics —— 仅峰成分分析（Python 流水线，无编译产物）

set(PEAK_PYTHON "${CMAKE_SOURCE_DIR}/src/process_Statistics/.venv/Scripts/python.exe"
    CACHE FILEPATH "Python interpreter for the peak pipeline")

add_custom_target(peak_analysis
    COMMAND "${PEAK_PYTHON}" "${CMAKE_CURRENT_SOURCE_DIR}/peak_pipeline.py"
    WORKING_DIRECTORY "${CMAKE_SOURCE_DIR}"
    COMMENT "Run peak-only analysis -> result/peak"
    VERBATIM
)
