cmake_minimum_required(VERSION 3.25.2)
project(hungarian_migration)

# Set C++ standard
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

# Define executable name
set(HUNGARIAN_MIGRATION hungarian_migration)

# Find ament_cmake (required for colcon builds)
find_package(ament_cmake REQUIRED)

# Add executable
add_executable(
    ${HUNGARIAN_MIGRATION}
    src/Hungarian.cpp
    src/main.cpp
)

# Include directories (if headers are in 'include/')
target_include_directories(${HUNGARIAN_MIGRATION} PRIVATE include)

# Link libraries (if needed, add libraries here)
# target_link_libraries(${HUNGARIAN_MIGRATION} ...)

# Install targets
install(TARGETS ${HUNGARIAN_MIGRATION}
        RUNTIME DESTINATION bin
        LIBRARY DESTINATION lib
        ARCHIVE DESTINATION lib
)

# Install include directory
install(DIRECTORY include/
        DESTINATION include
)

# Option to build tests
option(BUILD_TESTS "Build tests" ON)

# Conditionally build tests if the option is enabled
if(BUILD_TESTS)
  message(STATUS "[HUNGARIAN] - Building Tests")
  enable_testing()
  add_subdirectory(test)
endif()

# Export package for ament
ament_package()