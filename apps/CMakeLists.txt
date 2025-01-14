add_executable(App App.cxx)

target_compile_features(App PRIVATE cxx_std_17)

target_link_libraries(App PRIVATE ROOT::Core ROOT::RIO ROOT::Tree Logger Reader Parser Manager)
