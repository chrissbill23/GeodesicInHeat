cmake_minimum_required(VERSION 3.1)

set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} ${CMAKE_CURRENT_SOURCE_DIR}/cmake)

# libigl
option(LIBIGL_WITH_OPENGL            "Use OpenGL"         ON)
option(LIBIGL_WITH_OPENGL_GLFW       "Use GLFW"           ON)

find_package(LIBIGL REQUIRED QUIET)

file(GLOB SRCFILES src/*.cpp)
file(GLOB HEADERSFILES include/*.h)
add_executable(main ${HEADERSFILES} ${SRCFILES})

target_link_libraries(main igl::core igl::opengl_glfw )
include_directories(include/ src/)

