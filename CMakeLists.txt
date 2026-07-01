cmake_minimum_required(VERSION 3.20)
project(TreeSim)

set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

find_package(OpenGL REQUIRED)
find_package(glm CONFIG REQUIRED)
find_package(glfw3 CONFIG QUIET)

if (NOT glfw3_FOUND)
    find_package(PkgConfig REQUIRED)
    pkg_check_modules(GLFW REQUIRED glfw3)
endif()

add_library(glad
    external/glad/src/glad.c
)

target_include_directories(glad PUBLIC
    external/glad/include
)

add_library(TreeSimCore
    src/simulation/world.cpp
    src/simulation/organism.cpp
    src/simulation/limb.cpp
)
target_include_directories(TreeSimCore PUBLIC
    src/simulation
)

add_library(TreeSimRenderer
    src/renderer/renderer.cpp
    src/renderer/shader.cpp
    src/renderer/camera.cpp
)
target_include_directories(TreeSimRenderer PUBLIC
    src/renderer
)

if (TARGET glfw)
    target_link_libraries(TreeSimRenderer PUBLIC glfw)
else()
    target_link_libraries(TreeSimRenderer PUBLIC ${GLFW_LIBRARIES})
endif()
target_link_libraries(TreeSimRenderer PUBLIC
    glad
    OpenGL::GL
    TreeSimCore
)

# GLFW linking
if (TARGET glfw)
    set(GLFW_LIB glfw)
else()
    set(GLFW_LIB ${GLFW_LIBRARIES})
endif()

add_executable(TreeSim src/main.cpp)

target_link_libraries(TreeSim PRIVATE
    OpenGL::GL
    glm::glm
    TreeSimCore
    TreeSimRenderer
)