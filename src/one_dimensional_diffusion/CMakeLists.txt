set(main_target "OneDimensionalDiffusion")

add_executable(${main_target} main.cpp)
target_sources(${main_target}
  PRIVATE
    main.cpp
)

target_include_directories(${main_target}
  PUBLIC ${PROJECT_SOURCE_DIR}/include
)

target_link_libraries(${main_target} PRIVATE
  Base
)

install(TARGETS ${main_target} DESTINATION ${CMAKE_INSTALL_PREFIX})