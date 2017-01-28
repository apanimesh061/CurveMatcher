cmake_minimum_required(VERSION 3.5)
project(CurveMatcher)

message(STATUS "start running cmake...")

find_package(Boost 1.61.0 COMPONENTS system filesystem REQUIRED)

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11")

set(SOURCE_FILES Graph.cpp Graph.h filter.h)
add_executable(CurveMatcher ${SOURCE_FILES})

if (Boost_FOUND)
    message(STATUS "Boost_INCLUDE_DIRS: ${Boost_INCLUDE_DIRS}")
    message(STATUS "Boost_LIBRARIES: ${Boost_LIBRARIES}")
    message(STATUS "Boost_VERSION: ${Boost_VERSION}")
    include_directories(CurveMatcher ${Boost_INCLUDE_DIRS})
    target_link_libraries(CurveMatcher ${Boost_LIBRARIES})
endif ()
