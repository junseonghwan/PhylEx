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

On most compute clouds requirements can be satisfied by loading the modules. For example

```bash
module load gsl
module load cmake
module load boost
module load openmpi
module load gcc
```

## Instructions

To compile the code, it is first required to 
initialize [eigen](https://gitlab.com/libeigen/eigen) submodule:

```bash
git submodule init
git submodule update
```

To build:

```bash
mkdir build
cd build
cmake ..
make install
```

This will create `run`, `simul` and `testing` executables.
To run:

```bash
./run -c <config_file>
```

An example of the configuration file is provided in a separate repository for analysis [here](https://github.com/junseonghwan/PhylExAnalysis).

Note 1: this will install `libtssb.a` library file in your `local/bin` directory
unless you specify it with the `cmake` command option

```bash
cmake -DCMAKE_INSTALL_PREFIX=$HOME/opt ../src
```

Note 2: there seems to be an issue with ABI, especially with regards to the compilation of the `boost` library on a computer cloud.

Adding the following line 

```
add_compile_definitions(_GLIBCXX_USE_CXX11_ABI=0)
```

to `CMakeLists.txt` seemed to do the trick. Set to 1 or remove the line if any problem is encountered.

Alternatively, installing the latest version of the `boost` library may solve the problem as per [this issue](https://github.com/junseonghwan/PhylExAnalysis/issues/8).
