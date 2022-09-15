add_library(Base STATIC)
set(compile_options_tmp ${base_compile_options})
string(REPLACE " " ";" compile_options_tmp ${compile_options_tmp})

target_sources(Base
  PRIVATE
    matrix.cpp
    two_dimensional.cpp
    fileIO.cpp
)

target_include_directories(Base PUBLIC 
  ${PROJECT_SOURCE_DIR}/include
)

target_link_libraries(Base PRIVATE
  ${TP_LIBRARY}
  iomp5
)

target_compile_options(Base PUBLIC ${compile_options_tmp})