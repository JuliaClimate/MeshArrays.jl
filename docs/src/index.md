# MeshArrays.jl documentation

`MeshArrays.jl` mainly defines an array type that can contain / organize / distribute collections of inter-connected arrays as done in climate models (see `EarthGrids` plot below). The `MeshArray` type is a sub-type `AbstractArray` with an `outer array` where each element is itself a 2D `inner array`. This setup allows different choices for the outer and inner arrays -- for example `DistributedArrays` and `AxisArrays`, respectively. `MeshArrays.jl` thus provides a simple but efficient and general solution to analyze and simulate climate system variables.

The internals of a `MeshArray` are simply regulated by index ranges, array sizes, and inter-connections that are encoded in the `gcmgrid` struct. Such an approach is often useful in climate modeling which can involve advanced domain decompositions (see `EarthGrids` plot). The `exchange` methods readily transfer data between connected subdomains to allow for easy computaton of partial derivatives and related operators (e.g. gradients, curl, divergence) across subdomain edges though. This allows precise transport, budget, etc computations in climate models.

`MeshArrays.jl` was first introduced as as `gcmfaces.jl` in a [JuliaCon-2018 presentation](https://youtu.be/RDxAy_zSUvg). [This notebook folder](https://github.com/gaelforget/GlobalOceanNotebooks.git) demonstrates how its data structures can be used to accurately analyze the General Ocean Circulation. Examples include computations of [ocean heat transport](https://doi.org/10.1038/s41561-019-0333-7) and streamfunctions that are impportant and widely studied aspects of the climate system.

_Contents:_

```@contents
Pages = ["index.md","main.md","detail.md","API.md"]
Depth = 3
```

!!! note

    `MeshArrays.jl` is registered, documented, archived, and routinely tested, but also still regarded as a **preliminary implementation**.

## Install & Test

```
using Pkg
Pkg.add("MeshArrays")
Pkg.test("MeshArrays")
```

`Julia`'s package manager is documented [here within docs.julialang.org](https://docs.julialang.org/en/stable/stdlib/Pkg/).

## Use Examples

```
using MeshArrays
GridVariables=GridOfOnes("cs",6,100)
DemoVariables=MeshArrays.demo2(GridVariables)
```

The above example integrates lateral diffusion over the surface of a cube. The overall grid has 6 subdomains with `100 x 100` grid points each and all grid scales set to `1.0`. 

Alternatively, we can download and use a pre-defined [global ocean grid](http://www.geosci-model-dev.net/8/3071/2015/) from the [MITgcm](https://mitgcm.readthedocs.io/en/latest/) community.

```
git clone https://github.com/gaelforget/GRID_LLC90
GridVariables=GridLoad(GridSpec("LLC90"))
DemoVariables= MeshArrays.demo2(GridVariables)
```

This grid has 5 subdomains of uneven sizes (`105300` grid points in total), variable grid scale factors, and a realistic representation of  continents. 

## Earth Model Grids

![EarthGrids](https://raw.githubusercontent.com/gaelforget/MeshArrays.jl/master/docs/images/sphere_all.png)
