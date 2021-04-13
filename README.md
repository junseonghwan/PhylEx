# PhylEx: Integrated analysis of bulk and single cell RNA-seq

[![Build Status](https://travis-ci.com/junseonghwan/PhylEx.svg?token=wxZvzvzdwz1aU7zpr7vw&branch=master)](https://travis-ci.com/junseonghwan/PhylEx)

This software performs integrated analysis of bulk DNA-seq and single cell RNA-seq data. 

## Requirements

+ `gsl`
+ `boost`
+ `cmake`
+ `libomp`

## Instructions

To run the main software, just build the `src` code

```bash
mkdir build
cd build
cmake ../src
make install
```

then execute 

```bash
./run -c <config_file>
```

An example of the configuration file is provided in a separate repository for analysis [here](https://github.com/junseonghwan/PhylExAnalysis).

## Simulation and tests

To compile the simulation code, it is required first to 
initialize [eigen](https://gitlab.com/libeigen/eigen) submodule

```bash
git submodule init
git submodule update
```

then one can build all components with

```bash
cd PhylEx
mkdir build
cd build
cmake ..
make install
```

This will create `run`, `simul` and `testing` executables.