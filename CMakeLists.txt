cmake_minimum_required(VERSION 3.15...3.31)
project(Tree2Sexaquark LANGUAGES CXX)

# support folders in IDEs #

set_property(GLOBAL PROPERTY USE_FOLDERS ON)

# prevent in-source builds #

if(PROJECT_SOURCE_DIR STREQUAL PROJECT_BINARY_DIR)
  message(FATAL_ERROR "You don't want to configure in the source directory!")
endif()

# C++ standard #

set(CMAKE_CXX_EXTENSIONS OFF) # ensure -std=c++xx instead of -std=g++xx
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED True)

# find dependencies #

message(STATUS "Looking for ROOT")
find_package(ROOT 6.32 CONFIG REQUIRED COMPONENTS Core RIO Tree EG)
if(ROOT_FOUND)
  message(STATUS "Looking for ROOT -- found")
endif()

message(STATUS "Looking for KFParticle")
find_package(KFParticle CONFIG REQUIRED)
if(KFParticle_FOUND)
  message(STATUS "Looking for KFParticle -- found")
endif()

# main #

add_subdirectory(src)
add_subdirectory(apps)
