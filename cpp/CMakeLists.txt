cmake_minimum_required(VERSION 3.16)
project(evimo2_pipeline LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

include_directories(include)

# Aplikacja: evimo2_viewer
add_executable(evimo2_viewer
    apps/viewer/main.cpp
)
target_include_directories(evimo2_viewer PRIVATE include)
