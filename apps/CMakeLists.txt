add_executable(App App.cxx)

target_compile_features(App PRIVATE cxx_std_17)

target_link_libraries(App PRIVATE ROOT::Physics ROOT::TIOFeatures Logger Manager)
