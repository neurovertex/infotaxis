cmake_minimum_required(VERSION 2.6)
project(Infotaxis C CXX)

find_package(Boost REQUIRED)
find_package(PNG)
find_package (Threads REQUIRED)

if (PNG_FOUND)
	include_directories(${PNG_INCLUDE_DIRS} ${Boost_INCLUDE_DIRS})
elseif(NOT INFOTAXIS_NOPNG)
	message('Can't find Png++, ')
	SET(INFOTAXIS_NOPNG 'true')
endif(PNG_FOUND)


if (INFOTAXIS_NOPNG)
	add_definitions(-DINFOTAXIS_NOPNG)
endif (INFOTAXIS_NOPNG)

SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -pedantic -Wall")
SET(CMAKE_CXX_STANDARD 11)
SET(CMAKEC_XX_STANDARD_REQUIRED true)
SET(EXECUTABLE_OUTPUT_PATH ${PROJECT_BINARY_DIR}/bin)
SET(LIBRARY_OUTPUT_PATH ${PROJECT_BINARY_DIR}/lib)

add_library(infotaxis STATIC "src/infotaxis.cpp")

if (NOT INFOTAXIS_NOPNG)
	add_executable (main "src/main.cpp")
	target_link_libraries(main infotaxis ${PNG_LIBRARIES} ${Boost_LIBRARIES} ${CMAKE_THREAD_LIBS_INIT})
endif (NOT INFOTAXIS_NOPNG)


