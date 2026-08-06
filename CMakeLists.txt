cmake_minimum_required(VERSION 3.14)
project(H5ReaderTest CXX)

# Force C++17
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

# Force to Release build type if not specified
if(NOT CMAKE_BUILD_TYPE)
	set(CMAKE_BUILD_TYPE Release)
endif()

# Find the HDF5 library...
find_package(HDF5 REQUIRED COMPONENTS C)

# Define the source files for the test_h5_reader executable
set(SOURCE_FILES
	test_h5_reader.cxx
	src/ndlar_reader_app.cpp
	src/geometry.cpp
	src/interface.cpp
	src/hdf5/types.cpp
	src/hdf5/readers.cpp
	src/hdf5/readers/base.cpp
	src/hdf5/readers/packet_fraction.cpp
	src/hdf5/readers/region.cpp
	src/hdf5/readers/pairs.cpp
	src/hdf5/readers/segments.cpp
	src/hdf5/readers/trajectory.cpp
	src/hdf5/readers/interaction.cpp
	src/hdf5/readers/prompt_hit.cpp
	src/hdf5/collector.cpp
)
set(TARGET_LIBRARIES
	hdf5::hdf5
	HepEVD
	HighFive::HighFive
)

# Is ROOT installed? If so, find it and link against it.
find_package(ROOT QUIET)

if (ROOT_FOUND)
	message(STATUS "ROOT found: ${ROOT_VERSION}")
	include(${ROOT_USE_FILE})

	# Add a preprocessor definition to indicate that ROOT is available.
	add_definitions(-DROOT_FOUND)

	# We can add the ROOT-dependent source files here.
	list(APPEND SOURCE_FILES src/to_root.cpp)
	list(APPEND TARGET_LIBRARIES ${ROOT_LIBRARIES})
else()
	message(WARNING "ROOT not found. ROOT-dependent features will be disabled.")
endif()

# Add our executable
add_executable(test_h5_reader ${SOURCE_FILES})

# Add the HighFive library
add_subdirectory(highfive)

# Add the HepEVD library
add_subdirectory(HepEVD)

# Run the HepEVD setup script to configure the library
# This is just a case of running HepEVD/get_extern_deps.sh to download the
# required runtime dependencies for HepEVD. The script will download the
# dependencies and place them in the HepEVD/extern directory.
set(HEPEVD_SETUP_SCRIPT ${CMAKE_CURRENT_SOURCE_DIR}/HepEVD/get_extern_deps.sh)

add_custom_target(hep_evd_setup
	COMMAND bash ${HEPEVD_SETUP_SCRIPT}
	WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/HepEVD
	COMMENT "Fetching HepEVD external dependencies"
)

if (NOT EXISTS ${CMAKE_CURRENT_SOURCE_DIR}/HepEVD/include/extern/httplib.h)
	message(WARNING
		"HepEVD extern dependencies do not appear to be present. "
		"Run: cmake --build <build_dir> --target hep_evd_setup"
	)
endif()


# Include headers
target_include_directories(test_h5_reader PRIVATE
	${CMAKE_CURRENT_SOURCE_DIR}/include
)

# Link the HDF5 library and its include directories
target_link_libraries(test_h5_reader PRIVATE ${TARGET_LIBRARIES})

# Library Export Code...
set(LIB_SOURCES ${SOURCE_FILES})
list(REMOVE_ITEM LIB_SOURCES test_h5_reader.cxx src/ndlar_reader_app.cpp src/geometry.cpp src/to_root.cpp)

set(LIB_LIBRARIES ${TARGET_LIBRARIES})
list(REMOVE_ITEM LIB_LIBRARIES HepEVD ${ROOT_LIBRARIES})

add_library(ndlar_hdf5_reader SHARED ${LIB_SOURCES})
target_link_libraries(ndlar_hdf5_reader PUBLIC ${LIB_LIBRARIES})
target_include_directories(ndlar_hdf5_reader PUBLIC
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>
    $<INSTALL_INTERFACE:include>
)

include(CMakePackageConfigHelpers)
include(GNUInstallDirs)

install(TARGETS ndlar_hdf5_reader
    EXPORT ndlar_hdf5_readerTargets
    LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR}
    ARCHIVE DESTINATION ${CMAKE_INSTALL_LIBDIR}
    RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR}
    INCLUDES DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}
)

install(DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/include/ndlar
    DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}
)

install(EXPORT ndlar_hdf5_readerTargets
    FILE ndlar_hdf5_readerTargets.cmake
    NAMESPACE ndlar_hdf5_reader::
    DESTINATION ${CMAKE_INSTALL_LIBDIR}/cmake/ndlar_hdf5_reader
)

file(WRITE "${CMAKE_CURRENT_BINARY_DIR}/ndlar_hdf5_readerConfig.cmake"
    "include(CMakeFindDependencyMacro)\n"
    "find_dependency(HDF5 COMPONENTS C)\n"
    "include(\"\${CMAKE_CURRENT_LIST_DIR}/ndlar_hdf5_readerTargets.cmake\")\n"
)

install(FILES "${CMAKE_CURRENT_BINARY_DIR}/ndlar_hdf5_readerConfig.cmake"
    DESTINATION ${CMAKE_INSTALL_LIBDIR}/cmake/ndlar_hdf5_reader
)
