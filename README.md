# PhylEx: Integrated analysis of bulk and single cell RNA-seq

[![Build Status](https://travis-ci.com/junseonghwan/PhylEx.svg?token=wxZvzvzdwz1aU7zpr7vw&branch=master)](https://travis-ci.com/junseonghwan/PhylEx)

This software performs integrated analysis of bulk DNA-seq and single cell RNA-seq data. 

## Requirements

+ `gsl`
+ `boost`
+ `cmake`
+ `libomp`
+ `pkg-config`

Recommended approach is to install [homebrew](https://brew.sh/) package manager for OSX.

```
brew install cmake
brew install pkg-config
brew install gsl
brew install boost
brew install libomp
```

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

_Note_: this will install `libtssb.a` library file in your `local/bin` directory
unless you specify it with the `cmake` command option
```bash
cmake -DCMAKE_INSTALL_PREFIX=$HOME/opt ../src
```

## Simulation and tests

To compile the simulation code, it is first required to 
initialize [eigen](https://gitlab.com/libeigen/eigen) submodule

```bash
git submodule init
git submodule update
```

then one can build all components with

```bash
mkdir build
cd build
cmake ..
make install
```

This will create `run`, `simul` and `testing` executables.
