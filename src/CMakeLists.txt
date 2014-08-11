set(CMAKE_CXX_FLAGS "-Wall -Werror -Wextra")

set(CMAKE_EXE_LINKER_FLAGS "-lm")

add_subdirectory(RS232Connector)

add_executable(main main.cpp)

target_link_libraries(main RS232Connector)
