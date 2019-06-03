cmake_minimum_required (VERSION 3.8)
project (cfd3d)
list(APPEND CMAKE_MODULE_PATH ${CMAKE_CURRENT_SOURCE_DIR}/CMake)
#MESSAGE(STATUS "CMAKE_MODULE_PATH: " ${CMAKE_MODULE_PATH})

set (CMAKE_CXX_STANDARD 11)

file(GLOB_RECURSE SOURCES src/*.cpp src/*.c)
add_executable(cfd3d ${SOURCES})
include_directories(src)

#make VERBOSE=1

cmake_policy(SET CMP0012 NEW)
find_package(OpenMP REQUIRED)
find_package(Boost COMPONENTS system filesystem REQUIRED)
find_package(NetCDF REQUIRED)
find_package(glm REQUIRED)

include_directories(${NETCDF_INCLUDES})

target_link_libraries(cfd3d ${Boost_FILESYSTEM_LIBRARY} ${Boost_SYSTEM_LIBRARY} ${NETCDF_LIBRARIES})

if(OPENMP_FOUND)
set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS}")
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} ${OpenMP_EXE_LINKER_FLAGS}")
endif()
