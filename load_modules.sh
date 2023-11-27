#!/bin/bash

module load cmake
module load boost
module load gsl
module load openmpi
module load gcc/11.2.0

export CC=/software/gcc/11.2.0/b1/bin/gcc
export CXX=/software/gcc/11.2.0/b1/bin/g++
echo $CC
echo $CXX

