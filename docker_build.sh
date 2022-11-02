#!/bin/sh
mkdir -p build
cd build
cmake -D compiler=intel \
    -D TP_DIR=/opt/TextParser \
    -D CMAKE_INSTALL_PREFIX=/work \
    ..
make && make install