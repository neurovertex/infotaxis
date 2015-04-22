cmake_minimum_required(VERSION 2.6)
project(Infotaxis C CXX)
#set(BOOST_ROOT "/usr/include/boost")
find_package(Boost REQUIRED)
find_package(PNG REQUIRED)
find_package (Threads REQUIRED)

include_directories(${PNG_INCLUDE_DIRS} ${Boost_INCLUDE_DIRS})

SET(CMAKE_CXX_FLAGS "-pedantic -Wall")
SET(CMAKE_CXX_STANDARD 11)
SET(CMAKEC_XX_STANDARD_REQUIRED true)
SET(EXECUTABLE_OUTPUT_PATH ${PROJECT_BINARY_DIR}/bin)
SET(LIBRARY_OUTPUT_PATH ${PROJECT_BINARY_DIR}/lib)

add_library(infotaxis STATIC "src/infotaxis.cpp")
add_executable (main "src/main.cpp")

target_link_libraries(main infotaxis ${PNG_LIBRARIES} ${Boost_LIBRARIES} ${CMAKE_THREAD_LIBS_INIT})