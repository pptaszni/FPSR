cmake_minimum_required( VERSION 2.8 )
project(FPSR)

include_directories(inc inc/RS232Connector)

add_subdirectory(lib)

include_directories(lib/gmock/include lib/gmock/gtest/include)

add_subdirectory(src)

add_subdirectory(tests)
