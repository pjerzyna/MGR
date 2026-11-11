cmake_minimum_required(VERSION 3.16)
project(evimo2_pipeline LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

if(NOT CMAKE_BUILD_TYPE)
    set(CMAKE_BUILD_TYPE Release)
endif()

add_compile_options(-Wall -Wextra)

# Aplikacja: evimo2_viewer
add_executable(evimo2_viewer apps/viewer/main.cpp)
target_include_directories(evimo2_viewer PRIVATE include)

# Aplikacja: imo_detect (główny detektor IMO C++)
add_executable(imo_detect apps/imo_detect/main.cpp)
target_include_directories(imo_detect PRIVATE include)

# Testy numeryczne
add_executable(test_numerics tests/test_numerics.cpp)
target_include_directories(test_numerics PRIVATE include)

