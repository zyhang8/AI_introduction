cmake_minimum_required(VERSION 3.15)
project(AI_introduction)

set(CMAKE_CXX_STANDARD 14)

set(SOURCE_FILES astar.cpp)
#set(SOURCE_FILES genetic.cpp)
add_executable(ruadhan ${SOURCE_FILES})
