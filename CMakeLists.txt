cmake_minimum_required(VERSION 3.10)

project(TP_projekt3)

include_directories (${AudioFile_SOURCE_DIR})

add_subdirectory(pybind11)
add_subdirectory(AudioFile)
add_subdirectory(matplotplusplus)
add_definitions (-DPROJECT_BINARY_DIR="${CMAKE_CURRENT_BINARY_DIR}")
file (COPY test-audio.wav DESTINATION ${CMAKE_CURRENT_BINARY_DIR})


pybind11_add_module(TP_projekt3 main.cpp)

target_link_libraries(TP_projekt3 PUBLIC matplot)
target_link_libraries (TP_projekt3 PUBLIC AudioFile)