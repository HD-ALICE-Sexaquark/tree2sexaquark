add_executable(App App.cxx)

target_link_libraries(App PRIVATE Parser Manager Writer Reader KalmanFilter Logger KFParticle::KFParticle ROOT::EG ROOT::Tree ROOT::RIO ROOT::Core)
target_include_directories(App PUBLIC
  ${ROOT_INCLUDE_DIRS}
  ${KFParticle_INCLUDE_DIR}
  $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}>
)
