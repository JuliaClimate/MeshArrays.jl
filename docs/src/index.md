# MeshArrays.jl documentation

`MeshArrays.jl` mainly defines an array type that can contain / organize / distribute collections of inter-connected arrays as done in climate models (see _Earth Model Grids_ below). The `MeshArray` type is a sub-type of `AbstractArray` with an `outer array` where each element is itself a 2D `inner array`. This setup allows different choices for the outer and inner arrays -- for example `DistributedArrays` and `AxisArrays`, respectively. `MeshArrays.jl` thus provides a simple but efficient and general solution to analyze and simulate climate system variables.

The internals of a `MeshArray` are simply regulated by index ranges, array sizes, and inter-connections that are encoded in the `gcmgrid` struct. Such an approach is often useful in climate modeling which can involve advanced domain decompositions (see _Earth Model Grids_). The `exchange` methods, however, readily transfer data between connected subdomains to allow for easy computaton of partial derivatives and related operators (e.g. gradients, curl, divergence) across subdomain edges. This allows precise transport, budget, etc computations using climate model output.

`MeshArrays.jl` was first introduced as as `gcmfaces.jl` in a [JuliaCon-2018 presentation](https://youtu.be/RDxAy_zSUvg). This [notebook folder](https://github.com/gaelforget/GlobalOceanNotebooks.git) demonstrates how its data structures can be used to accurately analyze the General Ocean Circulation. Examples include computations of [ocean heat transport](https://doi.org/10.1038/s41561-019-0333-7) and streamfunctions that are important and widely studied aspects of the climate system.

_Contents:_

```@contents
Pages = ["index.md","main.md","detail.md","API.md"]
Depth = 3
```

## Install & Test

```
using Pkg
Pkg.add("MeshArrays")
Pkg.test("MeshArrays")
```

`Julia`'s package manager, `Pkg.jl`, is documented the [main Julia doc](https://docs.julialang.org/en/v1/) and [here in details](https://julialang.github.io/Pkg.jl/v1/).

## Basic Examples

Examples below (1) generate a grid configuration, (2) seed a 2D field of random noise, (3) smooth out this field, and (4) plot subdomain arrays. Smoothing is done via a lateral diffusion equation through time to illustrate how `MeshArray` computes partial derivatives & transfers data between neighboring subdomains. Examples 2 & 3 illustrate grid configurations commonly used in global models.

**[A]** _16 subdomains_, with _40x40 grid points_ each, covering a _doubly periodic domain_

```
using MeshArrays; p=dirname(pathof(MeshArrays))
GridVariables=GridOfOnes("PeriodicDomain",16,20)

include(joinpath(p,"../examples/Demos.jl"))
(in,out,_,_)=demo2(GridVariables);
show(out)

using Plots; plotlyjs()
include(joinpath(p,"../examples/Plots.jl"))
heatmap(out,clims=(-0.25,0.25))
```

Grid scale noise           |  Smoothed noise
:------------------------------:|:---------------------------------:
![raw](https://raw.githubusercontent.com/gaelforget/MeshArrays.jl/master/docs/images/noise_raw_16tiles.png)  |  ![smooth](https://raw.githubusercontent.com/gaelforget/MeshArrays.jl/master/docs/images/noise_smooth_16tiles.png)

**[B]** _6 subdomains_, with _100x100 points_ each, covering the _six faces of a cube_

```
GridVariables=GridOfOnes("CubeSphere",6,100)
DemoVariables=demo2(GridVariables)
```

**[C]** Global Model Grid with _5 uneven subdomains_, _variable spacing_, & _continents_

_This requires downloading a pre-defined [global ocean grid](http://www.geosci-model-dev.net/8/3071/2015/) from the [MITgcm community](https://mitgcm.readthedocs.io/en/latest/)._

```
git clone https://github.com/gaelforget/GRID_LLC90
GridVariables=GridLoad(GridSpec("LatLonCap","GRID_LLC90/"))
show(GridVariables["Depth"])

DemoVariables=demo2(GridVariables)
heatmap(out,clims=(-0.25,0.25))
```

## Earth Model Grids

-![EarthGrids](https://raw.githubusercontent.com/gaelforget/MeshArrays.jl/master/docs/images/sphere_all.png)
