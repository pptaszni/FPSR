file(GLOB Endogenous_SRC
    ${CMAKE_CURRENT_SOURCE_DIR}/*.cpp)

add_library(Endogenous STATIC ${Endogenous_SRC})