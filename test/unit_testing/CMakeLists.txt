cmake_minimum_required(VERSION 3.25.2)

find_package(GTest REQUIRED QUIET)

set(HUNGARIAN_TEST hungarian_test)

FetchContent_MakeAvailable(googletest)

# Include CTest for test discovery
include(CTest)
include(GoogleTest)

# Add the test executable
add_executable(hungarian_test TestHungarian.cpp ${CMAKE_SOURCE_DIR}/src/Hungarian.cpp)

# Add the include directory for Hungarian.h
target_include_directories(hungarian_test PRIVATE ${CMAKE_SOURCE_DIR}/include)

# Link Google Test libraries
target_link_libraries(hungarian_test PRIVATE GTest::GTest GTest::Main)

# Enable test discovery
gtest_discover_tests(hungarian_test)