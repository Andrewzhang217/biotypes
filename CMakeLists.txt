cmake_minimum_required(VERSION 3.19)
include(cmake/Warnings.cmake)

project(biotypes)

set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY ${PROJECT_BINARY_DIR}/lib)
set(CMAKE_LIBRARY_OUTPUT_DIRECTORY ${PROJECT_BINARY_DIR}/lib)
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${PROJECT_BINARY_DIR}/bin)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

add_library(biotypes 
  ${CMAKE_CURRENT_LIST_DIR}/src/lbcb/biotypes/bioseq.cc)
target_include_directories(biotypes 
  PUBLIC
    "$<BUILD_INTERFACE:${CMAKE_CURRENT_LIST_DIR}/include>")

target_compile_options(biotypes PUBLIC
  $<$<CONFIG:Debug,RelWithDebInfo>:-fsanitize=address>
  $<$<CONFIG:Debug,RelWithDebInfo>:-fno-omit-frame-pointer>)
target_link_options(biotypes PUBLIC
  $<$<CONFIG:Debug,RelWithDebInfo>:-fsanitize=address>)

include(${CMAKE_CURRENT_LIST_DIR}/test/BiotypesTest.cmake)
