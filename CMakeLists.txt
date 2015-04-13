cmake_minimum_required(VERSION 2.6)
project(Infotaxis)
set(BOOST_ROOT "/usr/include/boost")
find_package(Boost REQUIRED)
find_package(PNG REQUIRED)

include_directories(${Boost_INCLUDE_DIRS} ${PNG_INCLUDE_DIRS})

SET(EXECUTABLE_OUTPUT_PATH ${PROJECT_BINARY_DIR}/bin)
SET(LIBRARY_OUTPUT_PATH ${PROJECT_BINARY_DIR}/lib)

add_library(infotaxis STATIC "src/infotaxis.cpp")
add_executable (main "src/main.cpp")

target_link_libraries(infotaxis ${Boost_LIBRARIES})
target_link_libraries(main infotaxis ${PNG_LIBRARIES})