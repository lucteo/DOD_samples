cmake_minimum_required(VERSION 2.8)
project(MyTest)

find_package(OpenGL REQUIRED)
find_package(GLUT REQUIRED)

set(CMAKE_CXX_FLAGS  "${CMAKE_CXX_FLAGS} -std=gnu++0x -Wno-deprecated-declarations -O2")


add_executable(gameobject gameobject.cpp)
target_link_libraries(gameobject ${OPENGL_LIBRARIES} ${GLUT_LIBRARY})

