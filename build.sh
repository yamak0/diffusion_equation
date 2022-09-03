#!/bin/sh
mkdir build
cd build
cmake -D compiler=intel \
-D TP_DIR=/mnt/c/share/lib/TextParser \
-D CMAKE_INSTALL_PREFIX=/home/syusaku625/work/diffusion_equation \
..

make && make install