set(main2_target "TwoDimensionalDiffusion")

set(compile_options_tmp ${base_compile_options})
string(REPLACE " " ";" compile_options_tmp ${compile_options_tmp})

add_executable(${main2_target} main.cpp)
target_sources(${main2_target}
  PRIVATE
    main.cpp
)



target_include_directories(${main2_target}
  PUBLIC ${PROJECT_SOURCE_DIR}/include
)

target_link_libraries(${main2_target} PRIVATE
  Base
  ${MKL_LIBRARIES}
  iomp5
)

target_compile_options(${main2_target} PUBLIC ${compile_options_tmp})

install(TARGETS ${main2_target} DESTINATION ${CMAKE_INSTALL_PREFIX})