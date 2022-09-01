#!/bin/sh
mkdir build
cd build
cmake -D compiler=intel \
-D CMAKE_INSTALL_PREFIX=/mnt/d/diffusion \
..

make && make install