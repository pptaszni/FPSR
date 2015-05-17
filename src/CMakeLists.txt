set(CMAKE_CXX_FLAGS "-Wall -Wextra")

set(CMAKE_EXE_LINKER_FLAGS "-lm")

add_subdirectory(RS232Connector)
add_subdirectory(CUDASample)

add_executable(main main.cpp)

target_link_libraries(main RS232Connector CUDASample)
