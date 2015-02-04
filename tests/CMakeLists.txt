file(GLOB TEST_FILES
    ${CMAKE_CURRENT_SOURCE_DIR}/*.cpp)

add_executable(ut ${TEST_FILES})

target_link_libraries(ut gmock gtest pthread)
