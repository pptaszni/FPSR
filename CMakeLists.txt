cmake_minimum_required( VERSION 2.8 )
project(FPSR)

include_directories(inc inc/RS232Connector inc/CUDASample)

include_directories(/usr/local/cuda-6.0/samples/common/inc)

add_subdirectory(lib)

include_directories(lib/gmock/include lib/gmock/gtest/include tests/mocks)

add_subdirectory(src)

add_subdirectory(tests)
