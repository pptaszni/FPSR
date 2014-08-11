file(GLOB RS232Connector_SRC
    ${CMAKE_CURRENT_SOURCE_DIR}/*.cpp)

add_library(RS232Connector STATIC ${RS232Connector_SRC})
