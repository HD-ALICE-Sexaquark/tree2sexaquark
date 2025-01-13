cmake_minimum_required(VERSION 3.15...3.31)

project(
    Tree2Sexaquark
    VERSION 0.1
    DESCRIPTION "Tralala"
    LANGUAGES CXX)

set(CMAKE_CXX_EXTENSIONS OFF) # ensure -std=c++xx instead of -std=g++xx
set_property(GLOBAL PROPERTY USE_FOLDERS ON) # support folders in IDEs

# dependencies
find_package(ROOT 6.32 CONFIG REQUIRED COMPONENTS Core RIO Tree)

add_subdirectory(src) # compiled library code
add_subdirectory(apps) # executable code
