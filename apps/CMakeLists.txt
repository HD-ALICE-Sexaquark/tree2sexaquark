add_executable(App App.cxx)
target_link_libraries(App PUBLIC ROOT::Core KFParticle::KFParticle Logger Parser Manager)

message(STATUS "ROOT_INCLUDE_DIRS: ${ROOT_INCLUDE_DIRS}")
message(STATUS "KFParticle_INCLUDE_DIRXX: ${KFParticle_INCLUDE_DIRXX}")
message(STATUS "KFParticle_COMPILE_FLAGS: ${KFParticle_COMPILE_FLAGS}")
message(STATUS "KFParticle_INCLUDE_DIR: ${KFParticle_INCLUDE_DIR}")
message(STATUS "PROJECT_SOURCE_DIR: ${PROJECT_SOURCE_DIR}")

target_include_directories(App PUBLIC
  ${ROOT_INCLUDE_DIRS}
  ${KFParticle_INCLUDE_DIR}
  $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}>
)
set_target_properties(App PROPERTIES
  COMPILE_FLAGS "${KFParticle_COMPILE_FLAGS}"
)

