# MeshArrays.jl


[![Travis Build Status](https://travis-ci.org/gaelforget/MeshArrays.jl.svg?branch=master)](https://travis-ci.org/gaelforget/MeshArrays.jl)
[![codecov](https://codecov.io/gh/gaelforget/GCMFaces.jl/branch/master/graph/badge.svg)](https://codecov.io/gh/gaelforget/GCMFaces.jl)
[![Coverage Status](https://coveralls.io/repos/github/gaelforget/GCMFaces_jl/badge.svg?branch=master)](https://coveralls.io/github/gaelforget/GCMFaces_jl?branch=master)
[![DOI](https://zenodo.org/badge/143987632.svg)](https://zenodo.org/badge/latestdoi/143987632)

This repository contains the `MeshArrays.jl` package introduced at the [JuliaCon-2018](http://juliacon.org/2018/) conference by [this presentation](https://youtu.be/RDxAy_zSUvg). The code has passed a full test suite with `julia v0.7 and v1.0` but is still regarded as a **preliminary implementation**.

### Installation And Usage

To install this `Julia` package first execute `add MeshArrays` at the `pkg>` prompt. To use it then execute `using MeshArrays` at the main `>` REPL prompt or include this command in your modules or `startup.jl` file. `Julia`'s package manager, `Pkg`, is currently documented [here within docs.julialang.org](https://docs.julialang.org/en/stable/stdlib/Pkg/). Notebooks that illustrate the use `MeshArrays.jl` in practice are linked below.

### Main Package Features

`MeshArrays.jl` primarily defines composite types that embed inter-connected array collections within a `struct` and provides an `exchange` function that effectively transfers data between connected arrays. 

The composite types specify how each array collection forms a mesh and provide information to allow `exchange` to dispatch to the appropriate method. Various configurations that are commonly used in `Earth System Models` are implemented using the concrete type called `gcmfaces`. Embedded arrays, or meshes, each represent a subdomain within, e.g., an Earth System Model grid.

The `gcmfaces` name derives from a [previous Matlab / Octave package](https://gcmfaces.readthedocs.io/en/latest/) that was introduced in [Forget et al., 2015](http://www.geosci-model-dev.net/8/3071/2015/), `doi:10.5194/gmd-8-3071-2015`, and inspired this `Julia` package. Here, `GCM` is an acronym for General Circulation Model, or Global Climate Model as described in [this wikipedia entry](https://en.wikipedia.org/wiki/General_circulation_model), and `faces` is just another name for meshes, arrays, facets, or subdomains.


### Notebooks And Grids

The [JuliaCon-2018 presentation](https://youtu.be/RDxAy_zSUvg) relied on two `Jupyter notebooks` available in [this repository](https://github.com/gaelforget/JuliaCon2018Notebooks.git) (`demo_type.ipynb` and `demo_exch.ipynb`) and pre-defined grids available from [this website](http://mit.ecco-group.org/opendap/gforget/grid_examples/contents.html). These grids can, e.g., be downloaded as follows:

```
setenv DemoGrids 'ftp://mit.ecco-group.org/gforget/grid_examples/'
wget --recursive {$DemoGrids}/GRID_CS32/
wget --recursive {$DemoGrids}/GRID_LLC90/
mv mit.ecco-group.org/gforget/grid_examples/GRID_* .
```

But the [JuliaCon2018Notebooks repo](https://github.com/gaelforget/JuliaCon2018Notebooks.git) also contains `demo_smooth.ipynb` which illustrates how the `smooth` function is used for CI testing purposes and does not require downloading any of the predefined grids. 




