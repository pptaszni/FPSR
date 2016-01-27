file(GLOB Endogenous_SRC
    ${CMAKE_CURRENT_SOURCE_DIR}/*.cpp)

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11 -g")

add_library(Endogenous STATIC ${Endogenous_SRC})