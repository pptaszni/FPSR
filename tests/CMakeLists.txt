file(GLOB TEST_FILES
    ${CMAKE_CURRENT_SOURCE_DIR}/*.cpp)

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11")

add_executable(ut ${TEST_FILES})

target_link_libraries(ut gmock gtest pthread RS232Connector Endogenous)
