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

# Add our executable
add_executable(test_h5_reader
	test_h5_reader.cxx
	src/ndlar_reader_app.cpp
	src/geometry.cpp
	src/hdf5/types.cpp
	src/hdf5/readers.cpp
	src/hdf5/readers/packet_fraction.cpp
	src/hdf5/readers/region.cpp
	src/hdf5/readers/pairs.cpp
	src/hdf5/readers/segments.cpp
	src/hdf5/readers/trajectory.cpp
	src/hdf5/readers/interaction.cpp
	src/hdf5/readers/prompt_hit.cpp
	src/hdf5/collector.cpp
)

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
target_link_libraries(test_h5_reader PRIVATE
	hdf5::hdf5
	HepEVD
	HighFive::HighFive
)
