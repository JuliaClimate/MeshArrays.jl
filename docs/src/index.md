# MeshArrays.jl

`MeshArrays.jl` defines an array type that can contain / organize / distribute collections of inter-connected arrays as done in climate models (see below). 

As illustrated in the [Global Ocean Notebooks](https://github.com/JuliaClimate/GlobalOceanNotebooks.git), `MeshArrays`' data structures can be used to accurately analyze [ocean heat transport](https://doi.org/10.1038/s41561-019-0333-7), [material displacements](https://juliaclimate.github.io/IndividualDisplacements.jl/dev/), and many other important topics in climate science.

_Contents:_

```@contents
Pages = ["index.md","main.md","API.md","detail.md"]
Depth = 3
```

## Installation

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
γ,Γ=GridOfOnes("PeriodicDomain",16,20)

include(joinpath(p,"../examples/Demos.jl"))
(xi,xo,_,_)=demo2(Γ);
show(xo)

using Plots; plotlyjs()
include(joinpath(p,"../examples/Plots.jl"))
heatmap(xo,clims=(-0.25,0.25))
```

Grid scale noise           |  Smoothed noise
:------------------------------:|:---------------------------------:
![raw](https://raw.githubusercontent.com/gaelforget/MeshArrays.jl/master/docs/images/noise_raw_16tiles.png)  |  ![smooth](https://raw.githubusercontent.com/gaelforget/MeshArrays.jl/master/docs/images/noise_smooth_16tiles.png)

**[B]** _6 subdomains_, with _100x100 points_ each, covering the _six faces of a cube_

```
γ,Γ=GridOfOnes("CubeSphere",6,100)
D=demo2(Γ)
```

**[C]** Global Model Grid with _5 uneven subdomains_, _variable spacing_, & _continents_

_This requires downloading a pre-defined [global ocean grid](http://www.geosci-model-dev.net/8/3071/2015/) from the [MITgcm community](https://mitgcm.readthedocs.io/en/latest/)._

```
#run(`git clone https://github.com/gaelforget/GRID_LLC90`)
Γ=GridLoad(GridSpec("LatLonCap","GRID_LLC90/"))
D=demo2(Γ)
heatmap(D[2],clims=(-0.25,0.25))
```

## Earth Model Grids

-![EarthGrids](https://raw.githubusercontent.com/gaelforget/MeshArrays.jl/master/docs/images/sphere_all.png)

## JuliaCon 2018 Video

(where `MeshArrays.jl` was first introduced as as `gcmfaces.jl`)

[![JuliaCon-2018 presentation](https://user-images.githubusercontent.com/20276764/84893715-abe42180-b06d-11ea-92d3-173b678a701e.png)](https://youtu.be/W5DNqJG9jt0)


