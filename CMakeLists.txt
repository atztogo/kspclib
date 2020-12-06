cmake_minimum_required(VERSION 2.4)

project(kspclib C)
set(CMAKE_MACOSX_RPATH 1)
set(CMAKE_C_FLAGS_RELEASE "-Wall -O2")
set(CMAKE_C_FLAGS_DEBUG "-g -DTHMWARNING -DMATWARNING -DKPTWARNING")
if(NOT CMAKE_BUILD_TYPE)
  set(CMAKE_BUILD_TYPE Release)
endif(NOT CMAKE_BUILD_TYPE)
message(STATUS "Build type: ${CMAKE_BUILD_TYPE}")

include(GNUInstallDirs)

set(CMAKE_POSITION_INDEPENDENT_CODE ON)


# Version numbers
file(READ ${PROJECT_SOURCE_DIR}/src/version.h version_file)
string(REGEX MATCH "KSPCLIB_MAJOR_VERSION ([0-9]+)" kspclib_major_version ${version_file})
set(kspclib_major_version ${CMAKE_MATCH_1})
string(REGEX MATCH "KSPCLIB_MINOR_VERSION ([0-9]+)" kspclib_minor_version ${version_file})
set(kspclib_minor_version ${CMAKE_MATCH_1})
string(REGEX MATCH "KSPCLIB_MICRO_VERSION ([0-9]+)" kspclib_micro_version ${version_file})
set(kspclib_micro_version ${CMAKE_MATCH_1})
set(serial "${kspclib_major_version}.${kspclib_minor_version}.${kspclib_micro_version}")
set(soserial "1")

# Source code
include_directories("${PROJECT_SOURCE_DIR}/src")
set(SOURCES
  ${PROJECT_SOURCE_DIR}/src/kgengrid.c
  ${PROJECT_SOURCE_DIR}/src/kgrid.c
  ${PROJECT_SOURCE_DIR}/src/kpoint.c
  ${PROJECT_SOURCE_DIR}/src/kspclib.c
  ${PROJECT_SOURCE_DIR}/src/mathfunc.c
  ${PROJECT_SOURCE_DIR}/src/tetrahedron_method.c)

# Shared library
add_library(kspc SHARED ${SOURCES})

if(NOT MSVC)
  target_link_libraries(kspc m)
endif()
set_property(TARGET kspc PROPERTY VERSION ${serial})
set_property(TARGET kspc PROPERTY SOVERSION ${soserial})
install(TARGETS kspc LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR})

# Static link library
add_library(kspc_static STATIC ${SOURCES})
set_property(TARGET kspc_static PROPERTY VERSION ${serial})
set_property(TARGET kspc_static PROPERTY SOVERSION ${soserial})
set_property(TARGET kspc_static PROPERTY OUTPUT_NAME kspc)
install(TARGETS kspc_static ARCHIVE DESTINATION ${CMAKE_INSTALL_LIBDIR})

# Header file
install(FILES ${PROJECT_SOURCE_DIR}/src/kspclib.h DESTINATION ${CMAKE_INSTALL_INCLUDEDIR})

add_executable(kspclibtest EXCLUDE_FROM_ALL ${PROJECT_SOURCE_DIR}/test/test_kgrid.c ${SOURCES})
if(NOT MSVC)
  target_link_libraries(kspclibtest m)
endif()
enable_testing()
add_test(kspclibtest kspclibtest)

# make check

# cf. https://stackoverflow.com/questions/5709914/using-cmake-how-do-i-get-verbose-output-from-ctest
add_custom_target(check
        ${CMAKE_COMMAND} -E echo CWD=${CMAKE_BINARY_DIR}
        COMMAND ${CMAKE_COMMAND} -E echo CMD=${CMAKE_CTEST_COMMAND} -C $<CONFIG>
        COMMAND ${CMAKE_COMMAND} -E echo ----------------------------------
        COMMAND ${CMAKE_COMMAND} -E env CTEST_OUTPUT_ON_FAILURE=1
            ${CMAKE_CTEST_COMMAND} -C $<CONFIG> --verbose
    WORKING_DIRECTORY ${CMAKE_BINARY_DIR}
    DEPENDS kspclibtest
    )
