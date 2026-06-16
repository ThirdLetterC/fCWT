cmake_minimum_required(VERSION 3.16)

project(fCWT
    VERSION 2.0
    DESCRIPTION "Highly optimized implementation of the Continuous Wavelet Transform"
    LANGUAGES CXX)

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
option(FCWT_FETCH_KFR "Download KFR when an installed package is not found?" ON)

set(FCWT_KFR_GIT_TAG "7.0.1" CACHE STRING "KFR git tag used by FetchContent")
set(FCWT_CMAKE_INSTALL_DIR "${CMAKE_INSTALL_LIBDIR}/cmake/fCWT" CACHE STRING
    "Installation directory for CMake package files relative to CMAKE_INSTALL_PREFIX")

if(BUILD_SHARED_LIBS)
    set(CMAKE_POSITION_INDEPENDENT_CODE ON)
endif()

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
        $<INSTALL_INTERFACE:${CMAKE_INSTALL_INCLUDEDIR}>)

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

find_package(KFR CONFIG QUIET)
if(NOT KFR_FOUND)
    if(NOT FCWT_FETCH_KFR)
        message(FATAL_ERROR "KFR was not found. Set CMAKE_PREFIX_PATH to a KFR install or enable FCWT_FETCH_KFR.")
    endif()

    include(FetchContent)
    set(KFR_ENABLE_DFT ON CACHE BOOL "" FORCE)
    set(KFR_ENABLE_DSP OFF CACHE BOOL "" FORCE)
    set(KFR_ENABLE_IO OFF CACHE BOOL "" FORCE)
    set(KFR_ENABLE_AUDIO OFF CACHE BOOL "" FORCE)
    set(KFR_ENABLE_CAPI_BUILD OFF CACHE BOOL "" FORCE)
    FetchContent_Declare(kfr
        GIT_REPOSITORY https://github.com/kfrlib/kfr.git
        GIT_TAG ${FCWT_KFR_GIT_TAG})
    FetchContent_MakeAvailable(kfr)
endif()

target_link_libraries(${fcwt_library_target} PRIVATE kfr kfr_dft)

if(UNIX AND NOT APPLE)
    target_link_libraries(${fcwt_library_target} PRIVATE m)
endif()

install(FILES ${library_headers}
    DESTINATION ${CMAKE_INSTALL_INCLUDEDIR})

install(TARGETS ${fcwt_library_target}
    EXPORT ${fcwt_library_target}_export
    ARCHIVE DESTINATION ${CMAKE_INSTALL_LIBDIR}
    LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR}
    RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR}
    INCLUDES DESTINATION ${CMAKE_INSTALL_INCLUDEDIR})

install(EXPORT ${fcwt_library_target}_export
    NAMESPACE fCWT::
    DESTINATION ${FCWT_CMAKE_INSTALL_DIR}
    FILE fCWTConfig.cmake)
