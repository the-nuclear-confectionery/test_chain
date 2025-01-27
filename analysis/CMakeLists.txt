cmake_minimum_required(VERSION 3.12)
project(uspscape_converters)

set(CMAKE_CXX_STANDARD 20)

# Use CONDA_PREFIX as the base prefix
if(DEFINED ENV{CONDA_PREFIX})
    set(PREFIX $ENV{CONDA_PREFIX})
else()
    message(FATAL_ERROR "CONDA_PREFIX environment variable is not set. Please activate a Conda environment.")
endif()

# Debug output for the PREFIX variable
message(STATUS "Using CONDA_PREFIX: ${PREFIX}")

# Locate ROOT
find_package(ROOT 6 REQUIRED)

# Locate HepMC3
find_package(HepMC3 REQUIRED)

# Specify SMASH directories
set(SMASH_INCLUDE_DIR "${PREFIX}/include/smash-3.1")
set(SMASH_LIB_DIR "${PREFIX}/lib/smash-3.1")

# Add CUBA directory
set(CUBA_INCLUDE_DIR "${PREFIX}/include/smash-3.1/cuba")

# Include directories
include_directories(SYSTEM ${ROOT_INCLUDE_DIRS} ${SMASH_INCLUDE_DIR} ${CUBA_INCLUDE_DIR})

# Add the executable
add_executable(spectra_hepmc3.exe ${CMAKE_CURRENT_SOURCE_DIR}/spectra_hepmc3.cxx)

# Link SMASH and other required libraries
target_link_libraries(spectra_hepmc3.exe
    ${ROOT_LIBRARIES}
    ${HEPMC3_ROOTIO_LIB}
    -L${SMASH_LIB_DIR} -lsmash
    -lHepMC3 -lHepMC3rootIO
)


# Endianness detection
include(TestBigEndian)
TEST_BIG_ENDIAN(IS_BIG_ENDIAN)
if(IS_BIG_ENDIAN)
    message(STATUS "Big endian architecture detected.")
    add_definitions("-DBIG_ENDIAN_ARCHITECTURE")
else()
    message(STATUS "Little endian architecture detected.")
    add_definitions("-DLITTLE_ENDIAN_ARCHITECTURE")
endif()
