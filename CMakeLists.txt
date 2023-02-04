cmake_minimum_required(VERSION 3.23)
project(matrix)

set(CMAKE_CXX_STANDARD 11)

add_executable(matrix main.cpp matrix.h matrix.cpp test.h parser_evaluator.h parser_evaluator.cpp)
