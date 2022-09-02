add_library(Base STATIC)
target_sources(Base
  PRIVATE
    matrix.cpp
)

target_include_directories(Base
  PUBLIC ${PROJECT_SOURCE_DIR}/include
)