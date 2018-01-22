cmake_minimum_required(VERSION 3.7)
project(5G_LDPC)

set(CMAKE_CXX_STANDARD 11)
include_directories(./common/)


set(SOURCE_FILES main.cpp common/LDPC_helper.cpp common/LDPC_helper.h common/LDPC_5G.cpp common/LDPC_5G.h common/LDPC.h common/Channel.h common/GF.h)
add_executable(5G_LDPC ${SOURCE_FILES})