cmake_minimum_required (VERSION 3.1.1)

project(RUNTIME_ANALYSIS CXX)

# Gathering source files
file(GLOB_RECURSE source_files ../../source/*)
file(GLOB_RECURSE test_files ../../script/test_cases/*)
file(GLOB_RECURSE runtime_files ../../script/runtime_analysis/*)

# Include

# DACE
include_directories(../../../lib/dace/interfaces/cxx/include)
include_directories(../../../lib/dace/build/core/include)
include_directories(../../../lib/dace/core/include)

# Local files
include_directories(../../include)
include_directories(../../source)
include_directories(../../script/test_cases)
include_directories(../../script/runtime_analysis)
include_directories(../../script)

# Linking librairies (DACE)
link_directories(../../../lib/dace/build/interfaces/cxx)

# Main executable
add_executable(
	runtime_analysis
	${source_files} 
	${test_files}
	${test_files}
	../../script/runtime_analysis/run_runtime_analysis.cpp)

# Properties
set_target_properties(runtime_analysis PROPERTIES
    CXX_STANDARD 14)
	
# Link
target_link_libraries(
	runtime_analysis
    dace)
