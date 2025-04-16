cmake_minimum_required(VERSION 3.25.2)

# Include FetchContent to download and configure Google Test
include(FetchContent)

set(TestHungarian testhungarian)

# Download and configure Google Test
FetchContent_Declare(
  googletest
  URL https://github.com/google/googletest/archive/refs/tags/release-1.12.1.zip
)

# For Windows: Prevent overriding the parent project's compiler/linker settings
set(gtest_force_shared_crt ON CACHE BOOL "" FORCE)

# Make Google Test available
FetchContent_MakeAvailable(googletest)

# Include CTest for test discovery
include(CTest)
include(GoogleTest)

add_subdirectory(unit_testing)