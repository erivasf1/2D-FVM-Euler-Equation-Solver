cmake_minimum_required(VERSION 3.10)

# Project name
project(1DQuasiSteadyNozzle)

# Set C++ standard
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED True)

# Set the default build type to Debug if not specified
if(NOT CMAKE_BUILD_TYPE)
    set(CMAKE_BUILD_TYPE Debug)
endif()

# Enable debug symbols
set(CMAKE_CXX_FLAGS_DEBUG "-g")
set(CMAKE_C_FLAGS_DEBUG "-g")

# Define source files
set(SOURCES
    DataManager.cpp
    EulerOperator.cpp
    ExactNozzle.cpp
    MeshGen.cpp
    Main.cpp
    Output.cpp
    Test.cpp
)

# Define header files
set(HEADERS
    DataManager.h
    EulerOperator.h
    ExactNozzle.h
    MeshGen.h
    Output.h
    Test.h
)

# Create the executable
add_executable(nozzle ${SOURCES} ${HEADERS})

# Include directories (if needed)
#target_include_directories(NozzleSimulation PRIVATE ${CMAKE_SOURCE_DIR})

# Link libraries (if any)
# target_link_libraries(NozzleSimulation <LIBRARY_NAME>)
