#!/bin/bash
source modules_load.sh
rm -rf build
mkdir build
cd build
cM ffC
make -j $(nproc)
cd ..
