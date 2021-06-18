# MeshArrays.jl


[![Travis Build Status](https://travis-ci.com/juliaclimate/MeshArrays.jl.svg?branch=master)](https://travis-ci.com/juliaclimate/MeshArrays.jl)
[![](https://img.shields.io/badge/docs-stable-blue.svg)](https://juliaclimate.github.io/MeshArrays.jl/stable)
[![](https://img.shields.io/badge/docs-dev-blue.svg)](https://juliaclimate.github.io/MeshArrays.jl/dev)
[![codecov](https://codecov.io/gh/juliaclimate/MeshArrays.jl/branch/master/graph/badge.svg)](https://codecov.io/gh/juliaclimate/MeshArrays.jl)

[![DOI](https://zenodo.org/badge/143987632.svg)](https://zenodo.org/badge/latestdoi/143987632)
[![Binder](https://mybinder.org/badge_logo.svg)](https://mybinder.org/v2/gh/JuliaClimate/MeshArrays.jl/master)


**MeshArrays.jl** defines the `MeshArray` data structure (or type) that can contain, and distribute, collections of inter-connected arrays as generally done in climate models. This provides a simple yet efficient and general way to e.g. analyze climate system simulations.

### Installation

```
using Pkg
Pkg.add("MeshArrays")
Pkg.test("MeshArrays")
```

### Workflow Example

The diffusive smoother presented here as an example uses `MeshArrays.jl` to compute partial derivatives over a global domain / grid, which involves transfering data between neighboring subdomain arrays. In this workflow example, we 

1. generate a global grid decomposition
2. seed random noise across global domain
3. smooth out noise by applying diffusion globally
4. plots the `outer` array of subdomain / `inner` arrays

```
using MeshArrays; p=dirname(pathof(MeshArrays))
γ,Γ=GridOfOnes("PeriodicDomain",16,20)

include(joinpath(p,"../examples/Demos.jl"))
(xi,xo,_,_)=demo2(Γ);
show(xo)

using Plots
include(joinpath(p,"../examples/Plots.jl"))
heatmap(xo,clims=(-0.25,0.25),colorbar=false,tickfont = (4, :black))
```

<img src="https://user-images.githubusercontent.com/20276764/118325229-2d883d80-b4d1-11eb-953b-ddbb11bcfe1b.png" width="40%"> ===> <img src="https://user-images.githubusercontent.com/20276764/118325093-f31ea080-b4d0-11eb-8c6e-8cd0cc2cc255.png" width="40%">

### Global Grids

In the previous example we used a basic _doubly periodic_  domain with _16 subdomains_ of _40x40 grid points_ each. However, `MeshArrays.jl` also readily supports more elaborate global grid configurations, such as the ones shown below, which are commonly used in modeling climate.

<img src="docs/images/sphere_all.png" width="40%">

To be able to handle such climate model grids in practical and uniform fashion, `MeshArrays.jl` introduces custom array types that are geared towards climate science applications.

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

### Jupyter Notebooks

The [Global Ocean Notebooks](https://github.com/JuliaClimate/GlobalOceanNotebooks.git) illustrate standard workflows that use `MeshArrays.jl` to, e.g.:

- accurately compute planetary scale transports on a global ocean model [C-grid](https://en.wikipedia.org/wiki/Arakawa_grids)
- efficiently compute trajectories of ocean plastic, plankton, etc over any supported C-grid configuration using `MeshArrays.jl` along with [IndividualDisplacements.jl](https://github.com/JuliaClimate/IndividualDisplacements.jl)

Additional functionalities provided via related packages include:

- reading and writing files in [Netcdf](https://en.wikipedia.org/wiki/NetCDF) ([CF-compliant](http://cfconventions.org)), CSV, or binary formats often used in climate sciences. [NCTiles.jl](https://gaelforget.github.io/NCTiles.jl/stable/) readily supports domain decomposition and spatial interpolation when used along with `MeshArray.jl`
- support for the analysis of [MITgcm](https://mitgcm.readthedocs.io/en/latest/) model ouput and optimized, [state estimate](https://doi.org/10.5194/gmd-8-3071-2015) solutions are provided via [MITgcmTools.jl](https://github.com/gaelforget/MITgcmTools.jl) and [OceanStateEstimation.jl](https://github.com/gaelforget/OceanStateEstimation.jl), with interfaces in `MeshArray.jl`

_For information about Jupyter notebook, see [Jupyter docs](https://en.wikipedia.org/wiki/Project_Jupyter). Free apps like <https://nbviewer.jupyter.org> and <https://mybinder.org>  let you view them and try them out conveniently in the cloud._

[<img src="https://user-images.githubusercontent.com/20276764/84767001-b89a4400-af9f-11ea-956f-2e207f892c4f.png" width="40%">](https://youtu.be/M6vAUtIsIIY)

### JuliaCon 2018 Video

`MeshArrays.jl` was first introduced as as _gcmfaces.jl_ in [this presentation](https://youtu.be/RDxAy_zSUvg) at JuliaCon 2018 which corresponds to `DataStructures`/`01_MeshArrays.ipynb` in the [GlobalOceanNotebooks](https://github.com/JuliaClimate/GlobalOceanNotebooks.git)

[<img src="https://user-images.githubusercontent.com/20276764/84893715-abe42180-b06d-11ea-92d3-173b678a701e.png" width="40%">](https://youtu.be/RDxAy_zSUvg)



