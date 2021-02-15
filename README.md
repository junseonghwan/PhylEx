# PhylEx: Integrated analysis of bulk and single cell RNA-seq

[![Build Status](https://travis-ci.com/junseonghwan/PhylEx.svg?token=wxZvzvzdwz1aU7zpr7vw&branch=master)](https://travis-ci.com/junseonghwan/PhylEx)

This software performs integrated analysis of bulk DNA-seq and single cell RNA-seq data. 

## Requirements to build the software
+ gsl
+ boost
+ cmake
+ libomp

```c++
make build
cd build
cmake ..
make install
```

This should create `run`, `simul`, and `testing`.

```
./run.sh --config_file main.config
```

An example of main.config is provided in a separate repository for analysis [here](https://github.com/junseonghwan/PhylExAnalysis).
