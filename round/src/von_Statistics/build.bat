@echo off
setlocal
rem CUDA 12.8 + MSVC v143 (VS2026 default v145 is incompatible with CUDA 12.8).
rem On this machine nvcc lives in D:\Cpp\CUDA\12.8, not the VS-installed CUDA 13.3.
call "D:\Program Files\vs\Common7\Tools\VsDevCmd.bat" -arch=x64 -host_arch=x64 -vcvars_ver=14.40.33807 >nul
cd /d "%~dp0"
if "%1"=="clean" (
  if exist build rmdir /s /q build
)
cmake -S . -B build -G "NMake Makefiles" -DCMAKE_BUILD_TYPE=Release -DCMAKE_CUDA_COMPILER="D:/Cpp/CUDA/12.8/bin/nvcc.exe" -DCMAKE_CUDA_ARCHITECTURES=75
if errorlevel 1 exit /b 1
cmake --build build --parallel
if errorlevel 1 exit /b 1
echo BUILD_OK
