# ==============================================================================
# Messina Audio Processing Project
# ==============================================================================

cmake_minimum_required(VERSION 3.16)
project(Messina 
    VERSION 1.0.0
    DESCRIPTION "Audio processing application using AudioFile library"
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

# Add the AudioFile library
message(STATUS "Adding AudioFile library from external/AudioFile")
add_subdirectory(external/AudioFile)

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
        AudioFile
)

# ==============================================================================
# Custom Targets
# ==============================================================================

# Custom target to run the executable
add_custom_target(run
    COMMAND messina
    DEPENDS messina
    WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}
    COMMENT "Running Messina application"
)

# ==============================================================================
# Installation (Optional)
# ==============================================================================

# Install the executable
install(TARGETS messina
    RUNTIME DESTINATION bin
)