cmake_minimum_required(VERSION 3.25.2)

# Set C++ standard
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

# Define executable name
set(HUNGARIAN_MIGRATION hungarian_migration)

# Add executable
add_executable(
    ${HUNGARIAN_MIGRATION}
    src/Hungarian.cpp
    src/main.cpp
)

# Include directories (if headers are in 'include/')
target_include_directories(${HUNGARIAN_MIGRATION} PRIVATE include)