#---------------------------------------------------
# OpenMM Dynamo Plugin
#----------------------------------------------------
PROJECT(OPENMM-DYNAMO)

CMAKE_MINIMUM_REQUIRED(VERSION 3.0.0)

IF (NOT CMAKE_BUILD_TYPE)
SET(CMAKE_BUILD_TYPE "Release" CACHE STRING "Select build type" FORCE)
ENDIF()

SET(DYNAMO_SINGLE_PRECISION OFF CACHE BOOL "Compile with this floating point precision")

SET(CMAKE_C_FLAGS_RELEASE "-Ofast -DNDEBUG -D_FLOAT=float -std=c++11" CACHE STRING "C Optimization flags for release builds"  )
SET(CMAKE_CXX_FLAGS_RELEASE "-Ofast -DNDEBUG -D_FLOAT=float -std=c++11" CACHE STRING "C++ Optimization flags for release builds" )
SET(CUDA_NVCC_FLAGS_RELEASE "-Xptxas -O3" CACHE STRING "NVCC Optimization flags for release builds" FORCE)

# We need to know where OpenMM is installed so we can access the headers and libraries.
SET(OPENMM_DIR "/usr/local/openmm" CACHE PATH "Where OpenMM is installed")
INCLUDE_DIRECTORIES("${OPENMM_DIR}/include")
LINK_DIRECTORIES("${OPENMM_DIR}/lib" "${OPENMM_DIR}/lib/plugins")

# We also need to know where Dynamo is installed.
SET(DYNAMO_LIBRARY_DIR "$ENV{DYNAMOPLUGIN}" CACHE PATH "Where the Dynamo libraries are installed")
SET(DYNAMO_INCLUDE_DIR "$ENV{DYNAMOPLUGIN}" CACHE PATH "Where the Dynamo header files are installed")
INCLUDE_DIRECTORIES("${DYNAMO_INCLUDE_DIR}")
LINK_DIRECTORIES("${DYNAMO_LIBRARY_DIR}")

IF(DYNAMO_SINGLE_PRECISION)
    SET(EXTRA_COMPILE_FLAGS "${EXTRA_COMPILE_FLAGS} -D_FLOAT=float")
ELSE(DYNAMO_SINGLE_PRECISION)
    SET(EXTRA_COMPILE_FLAGS "${EXTRA_COMPILE_FLAGS} -D_FLOAT=double")
ENDIF(DYNAMO_SINGLE_PRECISION)

SET(EXTRA_COMPILE_FLAGS "${EXTRA_COMPILE_FLAGS} -Ofast -std=c++11" )
# set flags for linking on mac
IF(APPLE)
    SET(CMAKE_INSTALL_NAME_DIR "@rpath")
    SET(EXTRA_COMPILE_FLAGS "${EXTRA_COMPILE_FLAGS} -msse2 -stdlib=libc++")
ENDIF(APPLE)

# Select where to install
IF(${CMAKE_INSTALL_PREFIX_INITIALIZED_TO_DEFAULT})
    IF(WIN32)
        SET(CMAKE_INSTALL_PREFIX "$ENV{ProgramFiles}/OpenMM" CACHE PATH "Where to install the plugin" FORCE)
    ELSE(WIN32)
        SET(CMAKE_INSTALL_PREFIX "/usr/local/openmm" CACHE PATH "Where to install the plugin" FORCE)
    ENDIF(WIN32)
ENDIF(${CMAKE_INSTALL_PREFIX_INITIALIZED_TO_DEFAULT})

# The source is organized into subdirectories, but we handle them all from
# this CMakeLists file rather than letting CMake visit them as SUBDIRS.
SET(DYNAMO_PLUGIN_SOURCE_SUBDIRS openmmapi serialization)

# Set the library name
SET(DYNAMO_LIBRARY_NAME OpenMMDynamo)
SET(SHARED_DYNAMO_TARGET ${DYNAMO_LIBRARY_NAME})

# These are all the places to search for header files which are to be part of the API.
SET(API_INCLUDE_DIRS "openmmapi/include" "openmmapi/include/internal")

# Locate header files.
SET(API_INCLUDE_FILES)
FOREACH(dir ${API_INCLUDE_DIRS})
    FILE(GLOB fullpaths ${dir}/*.h)
    SET(API_INCLUDE_FILES ${API_INCLUDE_FILES} ${fullpaths})
ENDFOREACH(dir)

# Collect up source files
SET(SOURCE_FILES) # empty
SET(SOURCE_INCLUDE_FILES)
FOREACH(subdir ${DYNAMO_PLUGIN_SOURCE_SUBDIRS})
    FILE(GLOB src_files  ${CMAKE_CURRENT_SOURCE_DIR}/${subdir}/src/*.cpp)
    FILE(GLOB incl_files ${CMAKE_CURRENT_SOURCE_DIR}/${subdir}/src/*.h)
    SET(SOURCE_FILES         ${SOURCE_FILES}         ${src_files})   #append
    SET(SOURCE_INCLUDE_FILES ${SOURCE_INCLUDE_FILES} ${incl_files})

    ## Make sure we find these locally before looking in OpenMM/include if
    ## OpenMM was previously installed there.
    INCLUDE_DIRECTORIES(BEFORE ${CMAKE_CURRENT_SOURCE_DIR}/${subdir}/include)
ENDFOREACH(subdir)

SET(OPENMM_BUILD_C_AND_FORTRAN_WRAPPERS ON CACHE BOOL "Build wrappers for C and Fortran")
IF(OPENMM_BUILD_C_AND_FORTRAN_WRAPPERS)
    ADD_SUBDIRECTORY(wrappers)
    SET(SOURCE_FILES ${SOURCE_FILES} wrappers/DynamoOpenMMCWrapper.cpp  wrappers/DynamoOpenMMFortranWrapper.cpp)
    SET_SOURCE_FILES_PROPERTIES(wrappers/DynamoOpenMMCWrapper.cpp wrappers/DynamoOpenMMFortranWrapper.cpp PROPERTIES GENERATED TRUE)
ENDIF(OPENMM_BUILD_C_AND_FORTRAN_WRAPPERS)

# Create the library.

ADD_LIBRARY(${SHARED_DYNAMO_TARGET} SHARED ${SOURCE_FILES} ${SOURCE_INCLUDE_FILES} ${API_INCLUDE_FILES})
SET_TARGET_PROPERTIES(${SHARED_DYNAMO_TARGET}
    PROPERTIES COMPILE_FLAGS "-DDYNAMO_BUILDING_SHARED_LIBRARY ${EXTRA_COMPILE_FLAGS}"
    LINK_FLAGS "${EXTRA_COMPILE_FLAGS}")
TARGET_LINK_LIBRARIES(${SHARED_DYNAMO_TARGET} OpenMM dynamo)
INSTALL_TARGETS(/lib RUNTIME_DIRECTORY /lib ${SHARED_DYNAMO_TARGET})

# install headers
FILE(GLOB API_ONLY_INCLUDE_FILES "openmmapi/include/*.h")
INSTALL (FILES ${API_ONLY_INCLUDE_FILES} DESTINATION include)
FILE(GLOB API_ONLY_INCLUDE_FILES_INTERNAL "openmmapi/include/internal/*.h")
INSTALL (FILES ${API_ONLY_INCLUDE_FILES_INTERNAL} DESTINATION include/internal)

IF(OPENMM_BUILD_C_AND_FORTRAN_WRAPPERS)
    ADD_DEPENDENCIES(${SHARED_DYNAMO_TARGET} DynamoApiWrappers)
ENDIF(OPENMM_BUILD_C_AND_FORTRAN_WRAPPERS)

# Enable testing

ENABLE_TESTING()
ADD_SUBDIRECTORY(serialization/tests)

# Build the implementations for different platforms

ADD_SUBDIRECTORY(platforms/reference)

SET(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "${CMAKE_SOURCE_DIR}")
FIND_PACKAGE(OPENCL QUIET)
IF(OPENCL_FOUND)
    SET(DYNAMO_BUILD_OPENCL_LIB ON CACHE BOOL "Build implementation for OpenCL")
ELSE(OPENCL_FOUND)
    SET(DYNAMO_BUILD_OPENCL_LIB OFF CACHE BOOL "Build implementation for OpenCL")
ENDIF(OPENCL_FOUND)
IF(DYNAMO_BUILD_OPENCL_LIB)
    ADD_SUBDIRECTORY(platforms/opencl)
ENDIF(DYNAMO_BUILD_OPENCL_LIB)

FIND_PACKAGE(CUDA QUIET)
IF(CUDA_FOUND)
    SET(DYNAMO_BUILD_CUDA_LIB ON CACHE BOOL "Build implementation for CUDA")
ELSE(CUDA_FOUND)
    SET(DYNAMO_BUILD_CUDA_LIB OFF CACHE BOOL "Build implementation for CUDA")
ENDIF(CUDA_FOUND)
IF(DYNAMO_BUILD_CUDA_LIB)
    ADD_SUBDIRECTORY(platforms/cuda)
ENDIF(DYNAMO_BUILD_CUDA_LIB)

# Build the Python API

FIND_PROGRAM(PYTHON_EXECUTABLE python)
FIND_PROGRAM(SWIG_EXECUTABLE swig)
IF(PYTHON_EXECUTABLE AND SWIG_EXECUTABLE)
    SET(DYNAMO_BUILD_PYTHON_WRAPPERS ON CACHE BOOL "Build wrappers for Python")
ELSE(PYTHON_EXECUTABLE AND SWIG_EXECUTABLE)
    SET(DYNAMO_BUILD_PYTHON_WRAPPERS OFF CACHE BOOL "Build wrappers for Python")
ENDIF(PYTHON_EXECUTABLE AND SWIG_EXECUTABLE)
IF(DYNAMO_BUILD_PYTHON_WRAPPERS)
    ADD_SUBDIRECTORY(python)
ENDIF(DYNAMO_BUILD_PYTHON_WRAPPERS)
