#---------------------------------------------------
# OpenMM Struna Plugin
#----------------------------------------------------

CMAKE_MINIMUM_REQUIRED(VERSION 2.8)

IF (NOT CMAKE_BUILD_TYPE)
SET(CMAKE_BUILD_TYPE "Release" CACHE STRING "Select build type" FORCE)
ENDIF()

SET(CMAKE_C_FLAGS_RELEASE "-Ofast -DNDEBUG -std=c++11" CACHE STRING "C Optimization flags for release builds" FORCE)
SET(CMAKE_CXX_FLAGS_RELEASE "-Ofast -DNDEBUG -std=c++11" CACHE STRING "C++ Optimization flags for release builds" FORCE)
SET(CUDA_NVCC_FLAGS_RELEASE "-Xptxas -O3" CACHE STRING "NVCC Optimization flags for release builds" FORCE)

# We need to know where OpenMM is installed so we can access the headers and libraries.
SET(OPENMM_DIR "/usr/local/openmm" CACHE PATH "Where OpenMM is installed")
INCLUDE_DIRECTORIES("${OPENMM_DIR}/include")
LINK_DIRECTORIES("${OPENMM_DIR}/lib" "${OPENMM_DIR}/lib/plugins")

# We also need to know where Struna is installed.
SET(STRUNA_LIBRARY_DIR "$ENV{STRUNAPLUGIN}" CACHE PATH "Where the Struna libraries are installed")
SET(STRUNA_INCLUDE_DIR "$ENV{STRUNAPLUGIN}" CACHE PATH "Where the Struna header files are installed")
INCLUDE_DIRECTORIES("${STRUNA_INCLUDE_DIR}")
LINK_DIRECTORIES("${STRUNA_LIBRARY_DIR}")

# set flags for linking on mac
IF(APPLE)
    SET (CMAKE_INSTALL_NAME_DIR "@rpath")
    SET(EXTRA_COMPILE_FLAGS "-msse2 -stdlib=libc++")
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
SET(STRUNA_PLUGIN_SOURCE_SUBDIRS openmmapi serialization)

# Set the library name
SET(STRUNA_LIBRARY_NAME OpenMMStruna)
SET(SHARED_STRUNA_TARGET ${STRUNA_LIBRARY_NAME})

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
FOREACH(subdir ${STRUNA_PLUGIN_SOURCE_SUBDIRS})
    FILE(GLOB src_files  ${CMAKE_CURRENT_SOURCE_DIR}/${subdir}/src/*.cpp)
    FILE(GLOB incl_files ${CMAKE_CURRENT_SOURCE_DIR}/${subdir}/src/*.h)
    SET(SOURCE_FILES         ${SOURCE_FILES}         ${src_files})   #append
    SET(SOURCE_INCLUDE_FILES ${SOURCE_INCLUDE_FILES} ${incl_files})

    ## Make sure we find these locally before looking in OpenMM/include if
    ## OpenMM was previously installed there.
    INCLUDE_DIRECTORIES(BEFORE ${CMAKE_CURRENT_SOURCE_DIR}/${subdir}/include)
ENDFOREACH(subdir)

# Create the library.

ADD_LIBRARY(${SHARED_STRUNA_TARGET} SHARED ${SOURCE_FILES} ${SOURCE_INCLUDE_FILES} ${API_INCLUDE_FILES})
SET_TARGET_PROPERTIES(${SHARED_STRUNA_TARGET}
    PROPERTIES COMPILE_FLAGS "-DSTRUNA_BUILDING_SHARED_LIBRARY ${EXTRA_COMPILE_FLAGS}"
    LINK_FLAGS "${EXTRA_COMPILE_FLAGS}")
TARGET_LINK_LIBRARIES(${SHARED_STRUNA_TARGET} OpenMM struna)
INSTALL_TARGETS(/lib RUNTIME_DIRECTORY /lib ${SHARED_STRUNA_TARGET})

# install headers
FILE(GLOB API_ONLY_INCLUDE_FILES "openmmapi/include/*.h")
INSTALL (FILES ${API_ONLY_INCLUDE_FILES} DESTINATION include)
FILE(GLOB API_ONLY_INCLUDE_FILES_INTERNAL "openmmapi/include/internal/*.h")
INSTALL (FILES ${API_ONLY_INCLUDE_FILES_INTERNAL} DESTINATION include/internal)

# Enable testing

ENABLE_TESTING()
ADD_SUBDIRECTORY(serialization/tests)

# Build the implementations for different platforms

ADD_SUBDIRECTORY(platforms/reference)

SET(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "${CMAKE_SOURCE_DIR}")
FIND_PACKAGE(OpenCL QUIET)
IF(OPENCL_FOUND)
    SET(STRUNA_BUILD_OPENCL_LIB ON CACHE BOOL "Build implementation for OpenCL")
ELSE(OPENCL_FOUND)
    SET(STRUNA_BUILD_OPENCL_LIB OFF CACHE BOOL "Build implementation for OpenCL")
ENDIF(OPENCL_FOUND)
IF(STRUNA_BUILD_OPENCL_LIB)
    ADD_SUBDIRECTORY(platforms/opencl)
ENDIF(STRUNA_BUILD_OPENCL_LIB)

FIND_PACKAGE(CUDA QUIET)
IF(CUDA_FOUND)
    SET(STRUNA_BUILD_CUDA_LIB ON CACHE BOOL "Build implementation for CUDA")
ELSE(CUDA_FOUND)
    SET(STRUNA_BUILD_CUDA_LIB OFF CACHE BOOL "Build implementation for CUDA")
ENDIF(CUDA_FOUND)
IF(STRUNA_BUILD_CUDA_LIB)
    ADD_SUBDIRECTORY(platforms/cuda)
ENDIF(STRUNA_BUILD_CUDA_LIB)

# Build the Python API

FIND_PROGRAM(PYTHON_EXECUTABLE python)
FIND_PROGRAM(SWIG_EXECUTABLE swig)
IF(PYTHON_EXECUTABLE AND SWIG_EXECUTABLE)
    SET(STRUNA_BUILD_PYTHON_WRAPPERS ON CACHE BOOL "Build wrappers for Python")
ELSE(PYTHON_EXECUTABLE AND SWIG_EXECUTABLE)
    SET(STRUNA_BUILD_PYTHON_WRAPPERS OFF CACHE BOOL "Build wrappers for Python")
ENDIF(PYTHON_EXECUTABLE AND SWIG_EXECUTABLE)
IF(STRUNA_BUILD_PYTHON_WRAPPERS)
    ADD_SUBDIRECTORY(python)
ENDIF(STRUNA_BUILD_PYTHON_WRAPPERS)
