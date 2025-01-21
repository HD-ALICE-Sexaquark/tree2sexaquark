# Make an automatic library - will be static or dynamic based on user setting
add_library(Parser Parser.cxx)
add_library(Logger Logger.cxx)
add_library(Reader Reader.cxx)
add_library(Manager Manager.cxx)

target_include_directories(Parser PUBLIC
  $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/../include>)
target_include_directories(Logger PUBLIC
  $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/../include>)
target_include_directories(Reader PUBLIC
  $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/../include>)
target_include_directories(Manager PUBLIC
  $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/../include>)

target_link_libraries(Parser PUBLIC)
target_link_libraries(Logger PUBLIC ROOT::Core)
target_link_libraries(Reader PUBLIC ROOT::Core ROOT::RIO ROOT::Tree Logger)
target_link_libraries(Manager PUBLIC ROOT::EG Reader)

source_group(
  TREE ${PROJECT_SOURCE_DIR}/include
  FILES
    ${PROJECT_SOURCE_DIR}/include/Utilities/Logger.hxx
    ${PROJECT_SOURCE_DIR}/include/Trees/Structure.hxx
    ${PROJECT_SOURCE_DIR}/include/Trees/Reader.hxx
    ${PROJECT_SOURCE_DIR}/include/Analysis/Settings.hxx
    ${PROJECT_SOURCE_DIR}/include/Utilities/CL1.hpp
    ${PROJECT_SOURCE_DIR}/include/Utilities/Parser.hxx
    ${PROJECT_SOURCE_DIR}/include/Analysis/Manager.hxx)
