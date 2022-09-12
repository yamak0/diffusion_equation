set(main2_target "TwoDimensionalDiffusion")

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
)

install(TARGETS ${main2_target} DESTINATION ${CMAKE_INSTALL_PREFIX})