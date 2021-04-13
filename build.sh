#!/bin/bash

module load buildenv-gcc/2018a-eb
wait

module load CMake/3.12.1
wait

module load Boost/1.61.0-nsc1
wait

module load GSL/2.4-nsc1
wait

mkdir build
cd build || exit
cmake -DCMAKE_INSTALL_PREFIX:PATH=/home/x_vitza/opt ..
wait

make install

