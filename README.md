# GCMFaces.jl


[![Travis Build Status](https://travis-ci.org/gaelforget/MeshArrays.jl.svg?branch=master)](https://travis-ci.org/gaelforget/MeshArrays.jl)
[![codecov](https://codecov.io/gh/gaelforget/GCMFaces.jl/branch/master/graph/badge.svg)](https://codecov.io/gh/gaelforget/GCMFaces.jl)
[![Coverage Status](https://coveralls.io/repos/github/gaelforget/GCMFaces_jl/badge.svg?branch=master)](https://coveralls.io/github/gaelforget/GCMFaces_jl?branch=master)

This repository contains the `MeshArrays.jl` package introduced at the [JuliaCon-2018](http://juliacon.org/2018/) conference by [this presentation](https://youtu.be/RDxAy_zSUvg). The code has passed a full test suite with `julia v0.7 to  v1.0` but is still regarded as a **preliminary implementation**.

### Installation And Usage

To install the package first execute `add MeshArrays` command at `julia`'s `pkg>` prompt. To use it then execute `using MeshArrays` at `julia`'s `>` REPL prompt or include this command in your modules. `Julia`'s package manager, `Pkg`, is currently documented [here within docs.julialang.org](https://docs.julialang.org/en/stable/stdlib/Pkg/).

### Notebooks And Grids

The [JuliaCon-2018 presentation](https://youtu.be/RDxAy_zSUvg) relied on two `Jupyter notebooks` available in [this repository](https://github.com/gaelforget/JuliaCon2018Notebooks.git) (`demo_type.ipynb` and `demo_exch.ipynb`) and on pre-defined grids available from [this site](http://mit.ecco-group.org/opendap/gforget/grid_examples/contents.html) that can, e.g., be downloaded as follows:

```
setenv DemoGrids 'ftp://mit.ecco-group.org/gforget/grid_examples/'
wget --recursive {$DemoGrids}/GRID_CS32/
wget --recursive {$DemoGrids}/GRID_LLC90/
mv mit.ecco-group.org/gforget/grid_examples/GRID_* .
```

The [JuliaCon2018Notebooks repo](https://github.com/gaelforget/JuliaCon2018Notebooks.git) also contains `demo_smooth.ipynb` which illustrates how the `smooth` function is used for CI testing purposes (without the `GRID_*` input listed above). 

<!--- ### Package Name --->

<!--- GCM is an acronym that stands for General Circulation Model, or Global Climate Model, as discussed in [this wikipedia entry](https://en.wikipedia.org/wiki/General_circulation_model). The name for `GCMFaces.jl` derives from the `Matlab / Octave` package introduced in [Forget et al., 2015](http://www.geosci-model-dev.net/8/3071/2015/), `doi:10.5194/gmd-8-3071-2015` which inspired this `Julia` package. --->
