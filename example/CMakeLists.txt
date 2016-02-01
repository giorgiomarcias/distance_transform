cmake_minimum_required(VERSION 2.8)

project(distance_transform)

# source files
set(hdr_dir ${CMAKE_CURRENT_SOURCE_DIR}/..)
set(src
    ${hdr_dir}/multiple_array.hpp
    ${hdr_dir}/distance_transform.hpp
    ${CMAKE_CURRENT_SOURCE_DIR}/main.cpp
)

# create the executable
add_executable(${PROJECT_NAME} ${src})

set_target_properties(${PROJECT_NAME} PROPERTIES
    CXX_STANDARD 11
    CXX_STANDARD_REQUIRED ON
)

target_include_directories(${PROJECT_NAME} PUBLIC ${hdr_dir})
