# MeshArrays.jl documentation

`MeshArrays.jl` is a `Julia` package. It defines an Array type for collections of inter-connected arrays, and extends standard methods to readily operate on these `MeshArray`s. Its `exchange` methods transfer data between connected subdomains of the overall mesh. The internals of a `MeshArray` instance are regulated by index ranges, array sizes, and inter-connections that are encoded as a `gcmgrid` struct.

Such computational frameworks are commonly used in Earth System Modeling. `MeshArrays.jl` aims to provide a simple but versatile and powerful solution to this end. It was first introduced in this [JuliaCon-2018 presentation](https://youtu.be/RDxAy_zSUvg) as `gcmfaces.jl` (see [this other repo](https://github.com/gaelforget/JuliaCon2018Notebooks.git) for notebooks).

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
DemoVariables= MeshArrays.demo2(GridVariables)
```

The above will run an application which integrates lateral diffusion over the surface of a cube. This grid has 6 subdomains with `100 x 100` grid points each and all grid scales set to `1.0`. Alternatively,

```
git clone https://github.com/gaelforget/GRID_LLC90
GridVariables=GridLoad(GridSpec("LLC90"))
DemoVariables= MeshArrays.demo2(GridVariables)
```

would download and use a pre-defined [global ocean grid](http://www.geosci-model-dev.net/8/3071/2015/) from the [MITgcm](https://mitgcm.readthedocs.io/en/latest/) community. This grid has 5 subdomains of uneven sizes (`105300` grid points in total), variable grid scale factors, and a realistic representation of  continents. Results can be visualized as follows:

```
p=dirname(pathof(MeshArrays));
using Plots; include(joinpath(p,"Plots.jl"));
heatmap(Rend,title="smoothed noise",clims=(-0.5,0.5))
heatmap(Rini,title="raw noise",clims=(-0.5,0.5))
```