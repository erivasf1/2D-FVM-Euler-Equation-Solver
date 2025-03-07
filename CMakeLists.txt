cmake_minimum_required(VERSION 3.10)

# Project name
project(1DQuasiSteadyNozzle)

# Set C++ standard
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED True)

#set(CMAKE_PREFIX_PATH "/apps/easybuild/software/owl-genoa/Catch2/2.13.10-GCCcore-13.3.0/lib")
# Fetch Catch2 in the system
find_package(Catch2 REQUIRED)


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
    TimeIntegrator.cpp
)
# Define unit testing source files
set(TESTSOURCES
    DataManager_TEST.cpp
    EulerOperator_TEST.cpp
    ExactNozzle.cpp
    #Main_TEST.cpp #for system level testing
    MeshGen.cpp
    Output.cpp
    TimeIntegrator_TEST.cpp
    UnitTest.cpp #for unit testing
)

# Define header files
set(HEADERS
    DataManager.h
    EulerOperator.h
    ExactNozzle.h
    MeshGen.h
    Output.h
    Test.h
    TimeIntegrator.h
)

set(TESTHEADERS
    DataManager_TEST.h
    EulerOperator_TEST.h
    ExactNozzle.h
    MeshGen.h
    Output.h
    TimeIntegrator_TEST.h
)

# Create the executable
add_executable(nozzle ${SOURCES} ${HEADERS})
#add_executable(testnozzle UnitTest.cpp)
add_executable(testnozzle ${TESTSOURCES} ${TESTHEADERS})

# Include directories (if needed)
#target_include_directories(NozzleSimulation PRIVATE ${CMAKE_SOURCE_DIR})

# Link libraries (if any)
#target_link_libraries(nozzle PRIVATE ${CMAKE_PREFIX_PATH})
#target_link_libraries(NozzleSimulation <LIBRARY_NAME>)
