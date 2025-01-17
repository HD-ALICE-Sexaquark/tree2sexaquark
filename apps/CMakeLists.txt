add_executable(App App.cxx)

target_compile_features(App PRIVATE cxx_std_17)

include_directories(App ${CMAKE_SOURCE_DIR} ${KFParticle_INCLUDE_DIR} ${INCLUDE_DIRECTORIES})

target_link_libraries(App PRIVATE ROOT::Core ROOT::RIO ROOT::Tree KFParticle Logger Reader Parser Manager)
add_target_property(App COMPILE_FLAGS "-DDO_TPCCATRACKER_EFF_PERFORMANCE -DHomogeneousField -DUSE_TIMERS")
