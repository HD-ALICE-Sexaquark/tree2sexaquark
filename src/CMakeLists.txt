add_library(Logger STATIC Logger.cxx)
target_include_directories(Logger
  PUBLIC
    ${ROOT_INCLUDE_DIRS}
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/../include>
)

add_library(Parser STATIC Parser.cxx)
target_include_directories(Parser
  PUBLIC
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/../include>
)

add_library(Reader STATIC Reader.cxx)
target_include_directories(Reader
  PUBLIC
    ${ROOT_INCLUDE_DIRS}
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/../include>
)

add_library(KalmanFilter STATIC KalmanFilter.cxx)
target_include_directories(KalmanFilter
  PUBLIC
    ${ROOT_INCLUDE_DIRS}
    ${KFParticle_INCLUDE_DIR}
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/../include>
    )

add_library(Manager STATIC Manager.cxx)
target_include_directories(Manager
  PUBLIC
    ${ROOT_INCLUDE_DIRS}
    ${KFParticle_INCLUDE_DIR}
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/../include>
)

source_group(
  TREE ${PROJECT_SOURCE_DIR}/include
  FILES
    ${PROJECT_SOURCE_DIR}/include/Utilities/Logger.hxx
    ${PROJECT_SOURCE_DIR}/include/Trees/Structure.hxx
    ${PROJECT_SOURCE_DIR}/include/Trees/Reader.hxx
    ${PROJECT_SOURCE_DIR}/include/Analysis/Settings.hxx
    ${PROJECT_SOURCE_DIR}/include/Utilities/CL1.hpp
    ${PROJECT_SOURCE_DIR}/include/Utilities/Parser.hxx
    ${PROJECT_SOURCE_DIR}/include/Math/KalmanFilter.hxx
    ${PROJECT_SOURCE_DIR}/include/Analysis/Manager.hxx
)
