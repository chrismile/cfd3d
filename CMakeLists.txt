cmake_minimum_required (VERSION 3.8)
project (cfd3d)
set (CMAKE_CXX_STANDARD 11)
option(USE_MPI "USE_MPI" OFF)
option(USE_CUDA "USE_CUDA" OFF)
option(USE_SYCL "USE_SYCL" OFF)
list(APPEND CMAKE_MODULE_PATH ${CMAKE_CURRENT_SOURCE_DIR}/CMake)
#MESSAGE(STATUS "CMAKE_MODULE_PATH: " ${CMAKE_MODULE_PATH})

file(GLOB_RECURSE SOURCES src/*.cpp src/*.c)
file(GLOB_RECURSE CUDA_SOURCES src/*.cu)
include_directories(src)

if (USE_MPI)
    find_package(MPI)
    include_directories(SYSTEM ${MPI_INCLUDE_PATH})
    add_definitions(-DUSE_MPI)
else()
    list(REMOVE_ITEM SOURCES ${CMAKE_CURRENT_SOURCE_DIR}/src/CfdSolver/Mpi/BoundaryValuesMpi.cpp)
    list(REMOVE_ITEM SOURCES ${CMAKE_CURRENT_SOURCE_DIR}/src/CfdSolver/Mpi/CfdSolverMpi.cpp)
    list(REMOVE_ITEM SOURCES ${CMAKE_CURRENT_SOURCE_DIR}/src/CfdSolver/Mpi/MpiHelpers.cpp)
    list(REMOVE_ITEM SOURCES ${CMAKE_CURRENT_SOURCE_DIR}/src/CfdSolver/Mpi/SorSolverMpi.cpp)
    list(REMOVE_ITEM SOURCES ${CMAKE_CURRENT_SOURCE_DIR}/src/CfdSolver/Mpi/UvwMpi.cpp)
endif()

if (USE_SYCL)
    #add_compile_definitions(USE_SYCL) # Doesn't work on older versions of CMake
    add_definitions(-DUSE_SYCL)
else()
    list(REMOVE_ITEM SOURCES ${CMAKE_CURRENT_SOURCE_DIR}/src/CfdSolver/Sycl/BoundaryValuesSycl.cpp)
    list(REMOVE_ITEM SOURCES ${CMAKE_CURRENT_SOURCE_DIR}/src/CfdSolver/Sycl/CfdSolverSycl.cpp)
    list(REMOVE_ITEM SOURCES ${CMAKE_CURRENT_SOURCE_DIR}/src/CfdSolver/Sycl/SorSolverSycl.cpp)
    list(REMOVE_ITEM SOURCES ${CMAKE_CURRENT_SOURCE_DIR}/src/CfdSolver/Sycl/UvwSycl.cpp)
endif()

if (USE_CUDA)
    # For more details see https://devblogs.nvidia.com/building-cuda-applications-cmake/
    enable_language(CUDA)
    add_executable(cfd3d ${SOURCES} ${CUDA_SOURCES})
    set_target_properties(cfd3d PROPERTIES CUDA_SEPARABLE_COMPILATION ON)
    #add_compile_definitions(USE_CUDA) # Doesn't work on older versions of CMake
    add_definitions(-DUSE_CUDA)

else()
    add_executable(cfd3d ${SOURCES})
endif()

#make VERBOSE=1

cmake_policy(SET CMP0012 NEW)
find_package(Boost COMPONENTS system filesystem REQUIRED)
find_package(NetCDF REQUIRED)
find_package(glm REQUIRED)

if (CMAKE_CXX_COMPILER MATCHES "clang\\+\\+$")
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fopenmp=libiomp5 -std=c++11")
else()
    set (CMAKE_CXX_STANDARD 11)
    find_package(OpenMP REQUIRED)
    if(OPENMP_FOUND)
        set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS}")
        set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
        set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} ${OpenMP_EXE_LINKER_FLAGS}")
    endif()
endif()

include_directories(${NETCDF_INCLUDE_DIR})
include_directories(${Boost_INCLUDE_DIR})
include_directories(${GLM_INCLUDE_DIR})

if (USE_MPI)
    target_link_libraries(cfd3d ${Boost_FILESYSTEM_LIBRARY} ${Boost_SYSTEM_LIBRARY} ${NETCDF_LIBRARIES} ${MPI_LIBRARIES})
else()
    target_link_libraries(cfd3d ${Boost_FILESYSTEM_LIBRARY} ${Boost_SYSTEM_LIBRARY} ${NETCDF_LIBRARIES})
endif()

