cmake_minimum_required (VERSION 3.4)
project (try_idx C CXX)
set(default_build_type "Release")
set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wall -Wextra -std=c++17 -fPIC -fopenmp -W -Wall -pedantic -O3 -Ofast -march=native")

# ----------------------------------------------------------------------------
# kProcessor Setup
# ----------------------------------------------------------------------------

add_subdirectory(kProcessor)

set(kProcessor_INCLUDE_PATH "${PROJECT_SOURCE_DIR}/kProcessor/include")

# Add include directories.
include_directories("${PROJECT_SOURCE_DIR}/kProcessor/ThirdParty/MQF/include")
include_directories (${kProcessor_INCLUDE_PATH})


# ----------------------------------------------------------------------------
# END kProcessor Setup
# ----------------------------------------------------------------------------


add_executable (try_idx main.cpp)


target_link_libraries (try_idx kProcessor)
target_include_directories(try_idx INTERFACE ${kProcessor_INCLUDE_PATH})