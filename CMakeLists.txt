cmake_minimum_required(VERSION 3.1.0)
project(bff)

option(BFF_BUILD_CLI  "Build BFF command line" ON)
option(BFF_BUILD_GUI  "Build BFF GUI"          ON)
option(BFF_USE_CGAL   "Use CGAL library"       ON)

set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

set(CGAL_USE_BOOST 1)

# dependencies
set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "${CMAKE_CURRENT_SOURCE_DIR}/cmake" "${CMAKE_CURRENT_SOURCE_DIR}/deps/CGAL-6.1/cmake/Modules")

# define a variable for the path to the USD directory
set(USD_DIR "" CACHE PATH "Path to USD directory")
if(NOT USD_DIR STREQUAL "")
	set(CMAKE_CXX_STANDARD 17)
	include("${USD_DIR}/pxrConfig.cmake")
	MESSAGE(STATUS "PXR_INCLUDE_DIRS: ${PXR_INCLUDE_DIRS}")
	MESSAGE(STATUS ": ${PXR_LIBRARIES}")
endif()

#CGAL-6.1

if(BFF_USE_CGAL)
    # detect CGAL from ./deps
    set(CGAL_DIR "${CMAKE_CURRENT_SOURCE_DIR}/deps/CGAL-6.1")
    find_package(CGAL 6.1 REQUIRED)
    if(CGAL_FOUND)
        MESSAGE(STATUS "CGAL found: ${CGAL_DIR}")
        MESSAGE(STATUS "CGAL libraries: ${CGAL_LIBRARIES}")
    else()
        MESSAGE(FATAL_ERROR "CGAL-6.1 not found in ${CGAL_DIR}. Please check the path.")
    endif()
endif()

# suitesparse
# ------------------------------------------------------------------
# Detect SuiteSparse libraries:
# If not found automatically, set SuiteSparse_DIR in CMake to the
# directory where SuiteSparse was built.
# ------------------------------------------------------------------
set(SuiteSparse_USE_LAPACK_BLAS ON)
find_package(SuiteSparse QUIET NO_MODULE)  # 1st: Try to locate the *config.cmake file.
if(NOT SuiteSparse_FOUND)
	set(SuiteSparse_VERBOSE ON)
	find_package(SuiteSparse REQUIRED) # 2nd: Use FindSuiteSparse.cmake module
endif()
MESSAGE(STATUS "SuiteSparse_LIBS: ${SuiteSparse_LIBRARIES}")

# bff lib
list(APPEND BFF_INCLUDE_DIRS ${SuiteSparse_INCLUDE_DIRS})
list(APPEND BFF_INCLUDE_DIRS "deps/rectangle-bin-pack")
list(APPEND BFF_INCLUDE_DIRS "include")

file(GLOB BFF_SOURCES "deps/rectangle-bin-pack/GuillotineBinPack.cpp"
					  "deps/rectangle-bin-pack/SkylineBinPack.cpp"
					  "deps/rectangle-bin-pack/Rect.cpp"
					  "src/linear-algebra/*.cpp"
					  "src/mesh/*.cpp"
					  "src/project/*.cpp")

# create bff static library
add_library(bff ${BFF_SOURCES})
target_link_libraries(bff ${SuiteSparse_LIBRARIES})
target_include_directories(bff PUBLIC ${BFF_INCLUDE_DIRS})
if(NOT USD_DIR STREQUAL "")
	target_link_libraries(bff ${PXR_LIBRARIES})
	target_include_directories(bff PRIVATE ${PXR_INCLUDE_DIRS})
	target_compile_definitions(bff PRIVATE -DUSE_USD)
endif()

if (BFF_BUILD_CLI)
	# build command line app
	file(GLOB BFF_COMMAND_LINE_SOURCES "apps/command-line/src/*.cpp")

	add_executable(bff-command-line ${BFF_COMMAND_LINE_SOURCES})
	target_link_libraries(bff-command-line bff)
endif()

if (BFF_BUILD_GUI)
	# nanogui
	set(NANOGUI_BUILD_EXAMPLE OFF CACHE BOOL " " FORCE)
	set(NANOGUI_BUILD_PYTHON  OFF CACHE BOOL " " FORCE)
	set(NANOGUI_INSTALL       OFF CACHE BOOL " " FORCE)
	add_subdirectory(deps/nanogui)
	add_definitions(${NANOGUI_EXTRA_DEFS})

	# glm
	add_subdirectory(deps/glm)
	
	#CGAL
	if(BFF_USE_CGAL)
        		list(APPEND BFF_VIEWER_INCLUDE_DIRS ${CGAL_INCLUDE_DIRS})
    	endif()
	# viewer
	list(APPEND BFF_VIEWER_INCLUDE_DIRS ${NANOGUI_EXTRA_INCS})
	list(APPEND BFF_VIEWER_INCLUDE_DIRS "deps/nanogui/include")
	list(APPEND BFF_VIEWER_INCLUDE_DIRS "apps/viewer/include")

	# build viewer app
	file(GLOB BFF_VIEWER_SOURCES "apps/viewer/src/*.cpp")

	add_executable(bff-viewer ${BFF_VIEWER_SOURCES})
	target_include_directories(bff-viewer PUBLIC ${BFF_VIEWER_INCLUDE_DIRS})
	target_link_libraries(bff-viewer bff glm nanogui ${NANOGUI_EXTRA_LIBS})
	if(BFF_USE_CGAL)
        		target_link_libraries(bff-viewer CGAL::CGAL)
    	endif()
endif()
