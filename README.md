# PhylEx: Integrated analysis of bulk and single cell RNA-seq

This software performs integrated analysis of bulk DNA-seq and single cell RNA-seq data. 

## Requirements to build the software
+ gsl
+ boost
+ cmake

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