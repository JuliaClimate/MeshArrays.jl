# MeshArrays.jl


[![Travis Build Status](https://travis-ci.org/juliaclimate/MeshArrays.jl.svg?branch=master)](https://travis-ci.org/juliaclimate/MeshArrays.jl)
[![codecov](https://codecov.io/gh/juliaclimate/MeshArrays.jl/branch/master/graph/badge.svg)](https://codecov.io/gh/juliaclimate/MeshArrays.jl)
[![Coverage Status](https://coveralls.io/repos/github/juliaclimate/MeshArrays.jl/badge.svg?branch=master)](https://coveralls.io/github/juliaclimate/MeshArrays.jl?branch=master)

[![DOI](https://zenodo.org/badge/143987632.svg)](https://zenodo.org/badge/latestdoi/143987632)
[![](https://img.shields.io/badge/docs-stable-blue.svg)](https://juliaclimate.github.io/MeshArrays.jl/stable)
[![](https://img.shields.io/badge/docs-dev-blue.svg)](https://juliaclimate.github.io/MeshArrays.jl/dev)

`MeshArrays.jl` mainly defines an array type that can contain / organize / distribute collections of inter-connected arrays as done in climate models (see below & docs). It provides a simple but efficient and general solution to analyze and simulate climate system variables ([JuliaCon-2018 presentation](https://youtu.be/RDxAy_zSUvg)).

### Installation & Test

```
using Pkg
Pkg.add("MeshArrays")
Pkg.test("MeshArrays")
```

### Use Examples

Example below (1) generates a grid decomposition, (2) seeds random everywhere, (3) smoothes out the global map, and (4) plot subdomain arrays. The diffusion-based smoother illustrates how `MeshArray` computes partial derivatives over the global domain by transfering data between neighboring subdomains. 

```
using MeshArrays; p=dirname(pathof(MeshArrays))
GridVariables=GridOfOnes("dpdo",16,20)

include(joinpath(p,"../examples/Demos.jl"))
(in,out,_,_)=demo2(GridVariables);
show(out)

using Plots; plotlyjs()
include(joinpath(p,"../examples/Plots.jl"))
heatmap(out,clims=(-0.25,0.25))
```

Initial noise           |  Smoothed noise 
:------------------------------:|:---------------------------------:
![](docs/images/noise_raw_16tiles.png)  |  ![](docs/images/noise_smooth_16tiles.png)

Here we used _16 subdomains_, with _40x40 grid points_ each, covering a _doubly periodic global domain_. However, `MeshArrays` readily supports more elaborate grid configurations that are commonly used in global models.

![](docs/images/sphere_all.png)

### Jupyter Notebooks

The [JuliaCon-2018 presentation](https://youtu.be/RDxAy_zSUvg) corresponds to two of the [Jupyter notebooks](https://en.wikipedia.org/wiki/Project_Jupyter) found in the [GlobalOceanNotebooks repo](https://github.com/gaelforget/MeshArrayNotebooks.git) under `DataStructures/`. Other notebooks in this repo demonstrate how `MeshArrays` accurately computes ocean transport on native, cliamte model grids ([Arakawa C-grid](https://en.wikipedia.org/wiki/Arakawa_grids)). The [IndividualDisplacements.jl](https://github.com/gaelforget/IndividualDisplacements.jl) package in turn is able to compute trajectories (of plastic, plankton, etc) across the `C-grid configurations` supported in `MeshArrays`.



