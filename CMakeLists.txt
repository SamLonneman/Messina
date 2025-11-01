# ==============================================================================
# Messina Audio Processing Project
# ==============================================================================

# Minimum CMake version and project definition
cmake_minimum_required(VERSION 3.16)
project(Messina 
    VERSION 1.0.0
    DESCRIPTION "Audio processing application"
    LANGUAGES CXX
)

# ==============================================================================
# Project Configuration
# ==============================================================================

# Set C++ standard
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

# Set build type if not specified
if(NOT CMAKE_BUILD_TYPE)
    set(CMAKE_BUILD_TYPE Release)
endif()

# ==============================================================================
# Dependencies
# ==============================================================================

# Include PortAudio
add_subdirectory(external/portaudio)

# Force static build of RtMidi and include it
set(BUILD_SHARED_LIBS OFF CACHE BOOL "Build RtMidi as static library" FORCE)
add_subdirectory(external/rtmidi)

# ==============================================================================
# Target Definition
# ==============================================================================

# Create the main executable
add_executable(messina 
    messina.cpp
)

# Set target properties
set_target_properties(messina PROPERTIES
    OUTPUT_NAME "messina"
    CXX_STANDARD 17
)

# Link libraries
target_link_libraries(messina 
    PRIVATE 
        portaudio
        rtmidi
)

# ==============================================================================
# Installation (Optional)
# ==============================================================================

# Install the executable
install(TARGETS messina
    RUNTIME DESTINATION bin
)