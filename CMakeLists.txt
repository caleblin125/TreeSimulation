cmake_minimum_required(VERSION 3.20)
project(TreeSim)

set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

include(FetchContent)

find_package(OpenGL REQUIRED)
find_package(Threads REQUIRED)
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

FetchContent_Declare(
    imgui
    URL https://github.com/ocornut/imgui/archive/refs/tags/v1.91.9b.zip
    DOWNLOAD_EXTRACT_TIMESTAMP TRUE
)
FetchContent_MakeAvailable(imgui)

FetchContent_Declare(
    implot
    URL https://github.com/epezent/implot/archive/refs/tags/v0.16.zip
    DOWNLOAD_EXTRACT_TIMESTAMP TRUE
)
FetchContent_MakeAvailable(implot)

add_library(imgui STATIC
    ${imgui_SOURCE_DIR}/imgui.cpp
    ${imgui_SOURCE_DIR}/imgui_draw.cpp
    ${imgui_SOURCE_DIR}/imgui_tables.cpp
    ${imgui_SOURCE_DIR}/imgui_widgets.cpp
    ${imgui_SOURCE_DIR}/backends/imgui_impl_glfw.cpp
    ${imgui_SOURCE_DIR}/backends/imgui_impl_opengl3.cpp
)
target_include_directories(imgui PUBLIC
    ${imgui_SOURCE_DIR}
    ${imgui_SOURCE_DIR}/backends
)

if (TARGET glfw)
    target_link_libraries(imgui PUBLIC glfw)
else()
    target_include_directories(imgui PUBLIC ${GLFW_INCLUDE_DIRS})
    target_link_libraries(imgui PUBLIC ${GLFW_LIBRARIES})
endif()
target_link_libraries(imgui PUBLIC glad OpenGL::GL)

add_library(implot STATIC
    ${implot_SOURCE_DIR}/implot.cpp
    ${implot_SOURCE_DIR}/implot_items.cpp
)
target_include_directories(implot PUBLIC
    ${implot_SOURCE_DIR}
)
target_link_libraries(implot PUBLIC imgui)

add_library(TreeSimCore
    src/simulation/world.cpp
    src/simulation/organism.cpp
    src/simulation/limb.cpp
    src/simulation/geometry.cpp
    src/simulation/light.cpp
    src/simulation/physics.cpp
    src/simulation/distributed.cpp
)
target_include_directories(TreeSimCore PUBLIC
    include
)

option(TREESIM_USE_MPI "Split the world across MPI ranks when launched under mpirun" OFF)
if (TREESIM_USE_MPI)
    find_package(MPI REQUIRED COMPONENTS CXX)
    target_compile_definitions(TreeSimCore PRIVATE TREESIM_USE_MPI)
    target_link_libraries(TreeSimCore PUBLIC MPI::MPI_CXX)
endif()

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
    imgui
    implot
    OpenGL::GL
    Threads::Threads
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

enable_testing()

add_executable(TreeSimTests
    tests/simulation_tests.cpp
)

target_link_libraries(TreeSimTests PRIVATE
    glm::glm
    TreeSimCore
)

add_test(NAME TreeSimTests COMMAND TreeSimTests)
