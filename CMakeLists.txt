cmake_minimum_required(VERSION 3.10)

project(fCWT
    VERSION 2.0
    DESCRIPTION "Highly optimized implementation of the Continuous Wavelet Transform"
    LANGUAGES C CXX)

include(GNUInstallDirs)

if(NOT (UNIX OR WIN32))
    message(FATAL_ERROR "unsupported target platform")
endif()

if(NOT CMAKE_BUILD_TYPE AND NOT CMAKE_CONFIGURATION_TYPES)
    set(CMAKE_BUILD_TYPE Release CACHE STRING
        "Build type: Debug Release RelWithDebInfo MinSizeRel." FORCE)
endif()
message(STATUS "Build type is ${CMAKE_BUILD_TYPE}")

option(BUILD_SHARED_LIBS "Build shared library?" ON)
option(USE_OWN_FFTW "Use FFTW from the system or CMAKE_PREFIX_PATH instead of bundled libraries?" ON)
option(BUILD_EXAMPLE "Build example executable when src/main.cpp is available?" OFF)
option(BUILD_BENCHMARK "Build benchmark executable when benchmark sources are available?" OFF)
option(BUILD_MATLAB "Build Matlab mex?" OFF)

set(FCWT_MATLAB_DIR "../MATLAB" CACHE STRING
    "Name of the build directory for the MATLAB MEX-files")
set(FCWT_CMAKE_INSTALL_DIR "${CMAKE_INSTALL_LIBDIR}/cmake/fCWT" CACHE STRING
    "Installation directory for CMake package files relative to CMAKE_INSTALL_PREFIX")
set(FCWT_LIB_INSTALL_DIR "${CMAKE_INSTALL_LIBDIR}" CACHE STRING
    "Installation directory for libraries relative to CMAKE_INSTALL_PREFIX")
set(FCWT_INCLUDE_INSTALL_DIR "${CMAKE_INSTALL_INCLUDEDIR}" CACHE STRING
    "Installation directory for includes relative to CMAKE_INSTALL_PREFIX")

set(fcwt_library_target fCWT)

set(library_headers
    include/api.h
    include/morlet.h
    include/scales.h
    include/wavelet.h)

set(library_sources
    src/api.cpp
    src/morlet.cpp
    src/scales.cpp)

add_library(${fcwt_library_target} ${library_headers} ${library_sources})
add_library(fCWT::fCWT ALIAS ${fcwt_library_target})

target_compile_features(${fcwt_library_target} PUBLIC cxx_std_23)
set_target_properties(${fcwt_library_target}
    PROPERTIES
        VERSION ${PROJECT_VERSION}
        SOVERSION ${PROJECT_VERSION_MAJOR}
        CXX_STANDARD_REQUIRED ON)

target_include_directories(${fcwt_library_target}
    PUBLIC
        $<BUILD_INTERFACE:${PROJECT_SOURCE_DIR}/include>
        $<BUILD_INTERFACE:${PROJECT_SOURCE_DIR}/libs>
        $<INSTALL_INTERFACE:${FCWT_INCLUDE_INSTALL_DIR}>)

if(BUILD_SHARED_LIBS)
    target_compile_definitions(${fcwt_library_target} PRIVATE FCWT_LIBRARY_DLL_BUILDING)
endif()

if(MSVC)
    target_compile_options(${fcwt_library_target} PRIVATE /arch:AVX /W1)
else()
    target_compile_options(${fcwt_library_target} PRIVATE -mavx)
endif()

find_package(OpenMP REQUIRED)
target_link_libraries(${fcwt_library_target} PRIVATE OpenMP::OpenMP_CXX)

if(APPLE AND NOT USE_OWN_FFTW)
    set(_fcwt_fftw_names fftw3fmac)
    set(_fcwt_fftw_omp_names fftw3f_ompmac)
elseif(UNIX AND NOT APPLE AND NOT USE_OWN_FFTW)
    set(_fcwt_fftw_names fftw3fl)
    set(_fcwt_fftw_omp_names fftw3f_ompl)
else()
    set(_fcwt_fftw_names fftw3f libfftw3f libfftw3f.so.3)
    set(_fcwt_fftw_omp_names fftw3f_omp libfftw3f_omp libfftw3f_omp.so.3)
endif()

if(USE_OWN_FFTW)
    find_library(FFTW NAMES ${_fcwt_fftw_names} REQUIRED)
    if(UNIX)
        find_library(FFTW_OMP NAMES ${_fcwt_fftw_omp_names} REQUIRED)
    endif()
else()
    find_library(FFTW
        NAMES ${_fcwt_fftw_names}
        HINTS "${PROJECT_SOURCE_DIR}/libs"
        NO_DEFAULT_PATH
        REQUIRED)
    if(UNIX)
        find_library(FFTW_OMP
            NAMES ${_fcwt_fftw_omp_names}
            HINTS "${PROJECT_SOURCE_DIR}/libs"
            NO_DEFAULT_PATH
            REQUIRED)
    endif()
endif()

message(STATUS "FFTW: ${FFTW}")
target_link_libraries(${fcwt_library_target} PRIVATE ${FFTW})

if(UNIX)
    message(STATUS "FFTW_OMP: ${FFTW_OMP}")
    target_link_libraries(${fcwt_library_target} PRIVATE ${FFTW_OMP})
endif()

if(UNIX AND NOT APPLE)
    target_link_libraries(${fcwt_library_target} PRIVATE m)
endif()

install(FILES ${library_headers}
    DESTINATION ${FCWT_INCLUDE_INSTALL_DIR})

install(FILES libs/fftw3.h
    DESTINATION ${FCWT_INCLUDE_INSTALL_DIR})

install(TARGETS ${fcwt_library_target}
    EXPORT ${fcwt_library_target}_export
    ARCHIVE DESTINATION ${FCWT_LIB_INSTALL_DIR}
    LIBRARY DESTINATION ${FCWT_LIB_INSTALL_DIR}
    RUNTIME DESTINATION ${FCWT_LIB_INSTALL_DIR}
    INCLUDES DESTINATION ${FCWT_INCLUDE_INSTALL_DIR})

install(EXPORT ${fcwt_library_target}_export
    NAMESPACE fCWT::
    DESTINATION ${FCWT_CMAKE_INSTALL_DIR}
    FILE fCWTConfig.cmake)
