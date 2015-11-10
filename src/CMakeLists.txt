set(CMAKE_CXX_FLAGS "-Wall -Wextra")

set(OPENCV_LINKER_FLAGS "-lopencv_core -lopencv_imgproc -lopencv_highgui")

set(CMAKE_EXE_LINKER_FLAGS "-lm ${OPENCV_LINKER_FLAGS}")

add_subdirectory(RS232Connector)
add_subdirectory(CUDASample)

add_executable(main main.cpp)

target_link_libraries(main RS232Connector CUDASample)
