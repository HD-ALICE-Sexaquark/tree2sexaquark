# Make an automatic library - will be static or dynamic based on user setting
add_library(Parser Parser.cxx "${Tree2Sexaquark_SOURCE_DIR}/include/Utilities/Parser.hxx")
add_library(Logger Logger.cxx "${Tree2Sexaquark_SOURCE_DIR}/include/Utilities/Logger.hxx")
add_library(Reader Reader.cxx "${Tree2Sexaquark_SOURCE_DIR}/include/Trees/Reader.hxx")
add_library(Manager Manager.cxx "${Tree2Sexaquark_SOURCE_DIR}/include/Analysis/Manager.hxx")

# We need this directory, and users of our library will need it too
target_include_directories(Parser PUBLIC ../include)
target_include_directories(Logger PUBLIC ../include)
target_include_directories(Reader PUBLIC ../include)
target_include_directories(Manager PUBLIC ../include)

target_link_libraries(Parser PUBLIC)
target_link_libraries(Logger PUBLIC ROOT::Core)
target_link_libraries(Reader PUBLIC ROOT::Core ROOT::RIO ROOT::Tree Logger)
target_link_libraries(Manager PUBLIC ROOT::EG KFParticle Reader)

target_compile_features(Parser PUBLIC cxx_std_11)
target_compile_features(Logger PUBLIC cxx_std_11)
target_compile_features(Reader PUBLIC cxx_std_11)
target_compile_features(Manager PUBLIC cxx_std_11)

# IDEs should put the headers in a nice place
source_group(
  TREE "${PROJECT_SOURCE_DIR}/include"
  FILES
  "${PROJECT_SOURCE_DIR}/include/Utilities/Logger.hxx"
  "${PROJECT_SOURCE_DIR}/include/Trees/Structure.hxx"
  "${PROJECT_SOURCE_DIR}/include/Trees/Reader.hxx"
  "${PROJECT_SOURCE_DIR}/include/Analysis/Settings.hxx"
  "${PROJECT_SOURCE_DIR}/include/Utilities/CL1.hpp"
  "${PROJECT_SOURCE_DIR}/include/Utilities/Parser.hxx"
  "${PROJECT_SOURCE_DIR}/include/Analysis/Manager.hxx")
