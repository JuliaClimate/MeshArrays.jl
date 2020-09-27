# MeshArrays.jl


[![Travis Build Status](https://travis-ci.org/juliaclimate/MeshArrays.jl.svg?branch=master)](https://travis-ci.org/juliaclimate/MeshArrays.jl)
[![codecov](https://codecov.io/gh/juliaclimate/MeshArrays.jl/branch/master/graph/badge.svg)](https://codecov.io/gh/juliaclimate/MeshArrays.jl)
[![](https://img.shields.io/badge/docs-stable-blue.svg)](https://juliaclimate.github.io/MeshArrays.jl/stable)
[![](https://img.shields.io/badge/docs-dev-blue.svg)](https://juliaclimate.github.io/MeshArrays.jl/dev)
[![DOI](https://zenodo.org/badge/143987632.svg)](https://zenodo.org/badge/latestdoi/143987632)

`MeshArrays.jl` defines the `MeshArray` type that can contain, and distribute, collections of inter-connected arrays as generally done in climate models. This provides a simple yet efficient and general way to e.g. analyze climate system simulations.

```
struct gcmarray{T, N} <: AbstractMeshArray{T, N}
   grid::gcmgrid
   meta::varmeta
   f::Array{Array{T,2},N}
   fSize::Array{NTuple{2, Int}}
   fIndex::Array{Int,1}
   version::String
end
```

### Installation

```
using Pkg
Pkg.add("MeshArrays")
Pkg.test("MeshArrays")
```

### Use Examples

The example below (1) generates a grid decomposition, (2) seeds random noise everywhere, (3) smoothes out the noise, and (4) plots the (`outer`) array of subdomain (`inner`) arrays. The diffusion-based smoother illustrates how `MeshArrays.jl` computes partial derivatives over the whole domain by transfering data between neighboring subdomains. 

```
using MeshArrays; p=dirname(pathof(MeshArrays))
γ,Γ=GridOfOnes("PeriodicDomain",16,20)

include(joinpath(p,"../examples/Demos.jl"))
(xi,xo,_,_)=demo2(Γ);
show(xo)

using Plots; plotlyjs()
include(joinpath(p,"../examples/Plots.jl"))
heatmap(xo,clims=(-0.25,0.25))
```

<img src="docs/images/noise_smooth_16tiles.png" width="40%">

Above, we used _16 subdomains_, with _40x40 grid points_ each, covering a standard _doubly periodic domain_. However, `MeshArrays.jl` also readily supports elaborate grids commonly used in climate models, such as the ones shown below.

<img src="docs/images/sphere_all.png" width="40%">

### Jupyter Notebooks

The [Global Ocean Notebooks](https://github.com/JuliaClimate/GlobalOceanNotebooks.git) illustrate:

- Using `MeshArrays.jl` to accurately compute planetary transports on a ocean model [C-grid](https://en.wikipedia.org/wiki/Arakawa_grids).
- Using `MeshArrays.jl` with [IndividualDisplacements.jl](https://github.com/JuliaClimate/IndividualDisplacements.jl) to efficiently compute trajectories of ocean plastic, plankton, etc over the C-grid configurations supported by `MeshArrays.jl`.
- Support for [CF-compliant](http://cfconventions.org) [Netcdf](https://en.wikipedia.org/wiki/NetCDF) input / output of `MeshArray`s, with interpolation or domain decomposition, for `C-grid` variables as provided via [NCTiles.jl](https://gaelforget.github.io/NCTiles.jl/stable/).
- Support for [MITgcm](https://mitgcm.readthedocs.io/en/latest/) use cases and specificities is provided via [MITgcmTools.jl](https://github.com/gaelforget/MITgcmTools.jl).

[<img src="https://user-images.githubusercontent.com/20276764/84767001-b89a4400-af9f-11ea-956f-2e207f892c4f.png" width="40%">](https://youtu.be/M6vAUtIsIIY)

[(Jupyter notebook docs)](https://en.wikipedia.org/wiki/Project_Jupyter)

### JuliaCon 2018 Video

[<img src="https://user-images.githubusercontent.com/20276764/84893715-abe42180-b06d-11ea-92d3-173b678a701e.png" width="40%">](https://youtu.be/RDxAy_zSUvg)

where `MeshArrays.jl` was first introduced as as `gcmfaces.jl`. The [presentation](https://youtu.be/RDxAy_zSUvg) corresponds to [GlobalOceanNotebooks](https://github.com/JuliaClimate/GlobalOceanNotebooks.git) `/DataStructures/01_MeshArrays.ipynb`


