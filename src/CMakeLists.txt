add_library(Manager STATIC Manager.cxx)
target_include_directories(Manager
  PUBLIC
    ${ROOT_INCLUDE_DIRS}
    ${KFParticle_INCLUDE_DIR}
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/../include>)

add_library(KalmanFilter STATIC KalmanFilter.cxx)
target_include_directories(KalmanFilter
  PUBLIC
    ${ROOT_INCLUDE_DIRS}
    ${KFParticle_INCLUDE_DIR}
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/../include>)

add_library(Parser STATIC Parser.cxx)
target_include_directories(Parser
  PUBLIC
    ${ROOT_INCLUDE_DIRS}
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/../include>)

add_library(Settings STATIC Settings.cxx)
target_include_directories(Settings
  PUBLIC
    ${ROOT_INCLUDE_DIRS}
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/../include>)

add_library(Common STATIC Common.cxx)
target_include_directories(Common
  PUBLIC
    ${ROOT_INCLUDE_DIRS}
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/../include>)
