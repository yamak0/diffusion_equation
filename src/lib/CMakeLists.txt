add_library(Base STATIC)
target_sources(Base
  PRIVATE
    matrix.cpp
    two_dimensional.cpp
)

target_include_directories(Base PUBLIC 
  ${PROJECT_SOURCE_DIR}/include
)

target_link_libraries(Base PRIVATE
  ${TP_LIBRARY}
)