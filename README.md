# GCMFaces.jl


[![Travis Build Status](https://api.travis-ci.org/gaelforget/GCMFaces_jl.svg?branch=master)](https://travis-ci.org/gaelforget/GCMFaces_jl)
[![codecov](https://codecov.io/gh/gaelforget/GCMFaces_jl/branch/master/graph/badge.svg)](https://codecov.io/gh/gaelforget/GCMFaces_jl)
[![Coverage Status](https://coveralls.io/repos/github/gaelforget/GCMFaces_jl/badge.svg?branch=master)](https://coveralls.io/github/gaelforget/GCMFaces_jl?branch=master)

This repository contains the `GCMFaces.jl` package introduced at the [JuliaCon-2018](http://juliacon.org/2018/) conference by [this presentation](https://youtu.be/RDxAy_zSUvg). The provided code has been tested with `julia v0.7 to  v1.0` but is still regarded as a **preliminary implementation**.

### Installation And Usage

To activate `GCMFaces`  (e.g., by `using GCMFaces` at `julia`'s prompt), it first needs to be added to `julia` by executing the following command at the `pkg>` prompt: `add https://github.com/gaelforget/GCMFaces_jl`


### Notebooks And Grids

The [JuliaCon-2018 presentation](https://youtu.be/RDxAy_zSUvg) relied on two `Jupyter notebooks` available in [this repository](https://github.com/gaelforget/JuliaCon2018Notebooks.git) (`demo_type.ipynb` and `demo_exch.ipynb`) and on pre-defined grids available from [this site](http://mit.ecco-group.org/opendap/gforget/grid_examples/contents.html) that can, e.g., be downloaded as follows:

```
setenv DemoGrids 'ftp://mit.ecco-group.org/gforget/grid_examples/'
wget --recursive {$DemoGrids}/GRID_CS32/
wget --recursive {$DemoGrids}/GRID_LLC90/
mv mit.ecco-group.org/gforget/grid_examples/GRID_* .
```






