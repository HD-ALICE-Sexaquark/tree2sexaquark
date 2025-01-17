cmake_minimum_required(VERSION 3.15...3.31)

project(
    Tree2Sexaquark
    VERSION 0.1
    DESCRIPTION "Tralala"
    LANGUAGES CXX)

set(CMAKE_CXX_EXTENSIONS OFF) # ensure -std=c++xx instead of -std=g++xx
set_property(GLOBAL PROPERTY USE_FOLDERS ON) # support folders in IDEs

# dependencies
find_package(ROOT 6.32 CONFIG REQUIRED COMPONENTS Core RIO Tree EG)
find_package(KFParticle CONFIG REQUIRED)

include(${ROOT_USE_FILE})
set(CMAKE_INCLUDE_DIRECTORIES_BEFORE ON)
include_directories(${CMAKE_SOURCE_DIR} ${KFParticle_INCLUDE_DIR} ${INCLUDE_DIRECTORIES})
set(CMAKE_INCLUDE_DIRECTORIES_BEFORE OFF)
set(LINK_DIRECTORIES ${ROOT_LIBRARY_DIR})
link_directories(${LINK_DIRECTORIES})

ROOT_GENERATE_DICTIONARY(G__KFParticle ${HEADERS} LINKDEF KFLinkDef.h OPTIONS "-DDO_TPCCATRACKER_EFF_PERFORMANCE" "-DHomogeneousField" "-DUSE_TIMERS")
add_library(KFParticle SHARED ${SOURCES} G__KFParticle.cxx)

add_subdirectory(src) # compiled library code
add_subdirectory(apps) # executable code
