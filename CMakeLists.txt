cmake_minimum_required(VERSION 3.16)
project(PurineCFD_R2 LANGUAGES CXX)

set(CMAKE_EXPORT_COMPILE_COMMANDS ON)

find_package(OpenMP REQUIRED)
add_executable(purinecfd
  src/main.cpp src/config.cpp src/classconfig.cpp src/readmesh.cpp src/geometry.cpp
  src/initialize.cpp src/physic.cpp src/boundary.cpp src/interpolate.cpp src/grad.cpp
  src/convect.cpp src/dissipation.cpp src/SA.cpp src/timarch.cpp src/residual.cpp
  src/io.cpp src/parallel.cpp
)
target_compile_features(purinecfd PRIVATE cxx_std_17)
target_include_directories(purinecfd PRIVATE include)
target_include_directories(purinecfd SYSTEM PRIVATE third_party)
target_link_libraries(purinecfd PRIVATE OpenMP::OpenMP_CXX)
# Flow validity is checked explicitly; the solver does not inspect libm errno.
target_compile_options(purinecfd PRIVATE -Wall -Wextra -Wpedantic -Werror -fno-math-errno)

option(PURINE_NATIVE "Optimize for the build machine" OFF)
option(PURINE_SANITIZE "Enable AddressSanitizer and UBSan" OFF)
if(PURINE_NATIVE)
  target_compile_options(purinecfd PRIVATE -march=native)
endif()
if(PURINE_SANITIZE)
  target_compile_options(purinecfd PRIVATE -fsanitize=address,undefined -fno-omit-frame-pointer)
  target_link_options(purinecfd PRIVATE -fsanitize=address,undefined)
endif()
