cmake_minimum_required(VERSION 3.10)
set(project_name "GenerateAndTest")
set(CMAKE_MODULE_PATH ${CMAKE_SOURCE_DIR}/cmake_modules ${CMAKE_MODULE_PATH})
set(CMAKE_CXX_STANDARD 17)

# Load modules
include(SuperColliderServerPlugin RESULT_VARIABLE server_plugin_found)
if(NOT server_plugin_found)
    message(FATAL_ERROR "Could not find server plugin functions module")
endif()

include(SuperColliderCompilerConfig RESULT_VARIABLE compiler_config_found)
if(NOT compiler_config_found)
    message(FATAL_ERROR "Could not find compiler config module")
endif()

# Windows - puts redistributable DLLs in install directory
include(InstallRequiredSystemLibraries)

sc_check_sc_path("${SC_PATH}")
message(STATUS "Found SuperCollider: ${SC_PATH}")
set(SC_PATH "${SC_PATH}" CACHE PATH
    "Path to SuperCollider source. Relative paths are treated as relative to this script" FORCE)

include("${SC_PATH}/SCVersion.txt")
message(STATUS "Building plugins for SuperCollider version: ${SC_VERSION}")

project(${project_name})
sc_do_initial_compiler_config()

if(CMAKE_INSTALL_PREFIX_INITIALIZED_TO_DEFAULT OR NOT CMAKE_INSTALL_PREFIX)
    message(WARNING "No install prefix provided, defaulting to $BUILD_DIR/install")
    set(CMAKE_INSTALL_PREFIX "${CMAKE_BINARY_DIR}/install" CACHE PATH "Install prefix" FORCE)
endif()

message(STATUS "Install directory set to: ${CMAKE_INSTALL_PREFIX}")

# Options
option(SUPERNOVA "Build plugins for supernova" ON)
option(SCSYNTH "Build plugins for scsynth" ON)
option(NATIVE "Optimize for native architecture" OFF)
option(STRICT "Use strict warning flags" OFF)
option(NOVA_SIMD "Build plugins with nova-simd support." ON)

# Include directories for shared headers
include_directories(
    ${CMAKE_SOURCE_DIR}/external
    ${CMAKE_SOURCE_DIR}/plugins/Utils    
)

# Nova SIMD support
if (NOVA_SIMD)
    add_definitions(-DNOVA_SIMD)
    include_directories(${SC_PATH}/external_libraries/nova-simd)
endif()

# FilterBank Plugin
set(FilterBank_cpp_files
    plugins/FilterBank/FilterBank.cpp
)
set(FilterBank_sc_files
    plugins/FilterBank/Classes/FilterBank.sc
)
set(FilterBank_schelp_files
    plugins/FilterBank/HelpSource/FilterBank.schelp
)
sc_add_server_plugin(
    "GenerateAndTest/FilterBank" # destination directory
    "FilterBank" # target name
    "${FilterBank_cpp_files}"
    "${FilterBank_sc_files}"
    "${FilterBank_schelp_files}"
    ""
)

# SuperSaw Plugin
set(SuperSawOS_cpp_files
    plugins/SuperSawOS/SuperSawOS.cpp
)
set(SuperSawOS_sc_files
    plugins/SuperSawOS/Classes/SuperSawOS.sc
)
set(SuperSawOS_schelp_files
    plugins/SuperSawOS/HelpSource/SuperSawOS.schelp
)
sc_add_server_plugin(
    "GenerateAndTest/SuperSawOS" # destination directory
    "SuperSawOS" # target name
    "${SuperSawOS_cpp_files}"
    "${SuperSawOS_sc_files}"
    "${SuperSawOS_schelp_files}"
    ""
)

# Disperser Plugin
set(Disperser_cpp_files
    plugins/Disperser/Disperser.cpp
)
set(Disperser_sc_files
    plugins/Disperser/Classes/Disperser.sc
)
set(Disperser_schelp_files
    plugins/Disperser/HelpSource/Disperser.schelp
)
sc_add_server_plugin(
    "GenerateAndTest/Disperser" # destination directory
    "Disperser" # target name
    "${Disperser_cpp_files}"
    "${Disperser_sc_files}"
    "${Disperser_schelp_files}"
    ""
)

message(STATUS "Generating plugin targets done")