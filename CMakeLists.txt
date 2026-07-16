cmake_minimum_required(VERSION 3.14)
project(H5FlowReaderTest CXX)

# Force C++17
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

# Find the HDF5 library installed by Brew
find_package(HDF5 REQUIRED COMPONENTS C)

# Add our executable
add_executable(test_h5_reader test_h5_reader.cxx)

# Include the HighFive headers
target_include_directories(test_h5_reader PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/HighFive/include)

# Link the HDF5 library and its include directories
target_link_libraries(test_h5_reader PRIVATE hdf5::hdf5)
