cmake_minimum_required(VERSION 3.15...3.26)
project(${SKBUILD_PROJECT_NAME} LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 17 CACHE STRING "C++ version selection")
set(CMAKE_CXX_STANDARD_REQUIRED ON)

# Force macOS toolchains to use C++17
set(CMAKE_OSX_DEPLOYMENT_TARGET 10.15)

set(PYBIND11_FINDPYTHON ON)
find_package(pybind11 CONFIG REQUIRED)

set(Boost_DEBUG ON)
find_package(Boost 1.53 REQUIRED)
if(Boost_FOUND)
  include_directories(${Boost_INCLUDE_DIRS})
endif()

file(GLOB SRC_FILES ${SKBUILD_PROJECT_NAME}/extern/*.cpp)
file(GLOB HEADER_FILES ${SKBUILD_PROJECT_NAME}/extern/*.h)

pybind11_add_module(_core MODULE ${SRC_FILES} ${HEADER_FILES})
install(TARGETS _core DESTINATION ${SKBUILD_PROJECT_NAME})
