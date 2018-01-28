cmake_minimum_required(VERSION 3.9)
project(CpuFreeLunch)

set(CMAKE_CXX_FLAGS  "${CMAKE_CXX_FLAGS} -std=c++1z -ffast-math -O2 -Wno-deprecated-declarations")


add_executable(cpufreelunch cpufreelunch.cpp)
target_link_libraries(cpufreelunch)

