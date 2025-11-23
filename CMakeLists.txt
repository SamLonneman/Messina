cmake_minimum_required(VERSION 3.16)

project(Messina 
    VERSION 1.0.0
    DESCRIPTION "Audio processing application"
    LANGUAGES CXX
)

# Default to Release build
if(NOT CMAKE_BUILD_TYPE)
    set(CMAKE_BUILD_TYPE Release)
endif()

# ==============================================================================
# Dependencies
# ==============================================================================

if(NOT EXISTS "${CMAKE_SOURCE_DIR}/external/portaudio/CMakeLists.txt")
    message(FATAL_ERROR "PortAudio submodule missing. Run: git submodule update --init --recursive")
endif()

if(NOT EXISTS "${CMAKE_SOURCE_DIR}/external/rtmidi/CMakeLists.txt")
    message(FATAL_ERROR "RtMidi submodule missing. Run: git submodule update --init --recursive")
endif()

add_subdirectory(external/portaudio)

set(BUILD_SHARED_LIBS OFF CACHE BOOL "Build RtMidi as static library" FORCE)
add_subdirectory(external/rtmidi)

# ==============================================================================
# Target
# ==============================================================================

add_executable(messina 
    src/messina.cpp
)

set_target_properties(messina PROPERTIES
    CXX_STANDARD 17
    CXX_STANDARD_REQUIRED YES
    CXX_EXTENSIONS NO
)

if (MSVC)
    target_compile_options(messina PRIVATE /W4)
else()
    target_compile_options(messina PRIVATE -Wall -Wextra -Wpedantic)
endif()

target_link_libraries(messina 
    PRIVATE 
        portaudio
        rtmidi
)

install(TARGETS messina
    RUNTIME DESTINATION bin
)
