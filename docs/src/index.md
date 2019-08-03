# MeshArrays.jl documentation

`MeshArrays.jl` primarily defines composite types that embed inter-connected array collections within a `struct` and provides an `exchange` function that effectively transfers data between connected arrays. It was originally introduced, as `gcmfaces.jl`, in this [JuliaCon-2018 presentation](https://youtu.be/RDxAy_zSUvg) (see below for **notebooks**). _Note:_ `MeshArrays.jl` is registered, documented, archived, and routinely tested, but is still regarded as a **preliminary implementation**.

## Contents

```@contents
```

## Installation

```
using Pkg
Pkg.add("MeshArrays")
Pkg.test("MeshArrays")
```
_Note:_ `Julia`'s package manager is documented [here within docs.julialang.org](https://docs.julialang.org/en/stable/stdlib/Pkg/).

## Use examples

The [JuliaCon-2018 presentation](https://youtu.be/RDxAy_zSUvg) relied on two `Jupyter notebooks` that are available in the [JuliaCon2018Notebooks repo](https://github.com/gaelforget/JuliaCon2018Notebooks.git) and pre-defined grids that are available from [this website](http://mit.ecco-group.org/opendap/gforget/grid_examples/contents.html) and can also be downloaded from `github` as follows:

```
git clone https://github.com/gaelforget/GRID_CS32
git clone https://github.com/gaelforget/GRID_LLC90
```

After downloading `GRID_LLC90/`, you can also replicate the following example:

```
using MeshArrays

!isdir("GRID_LLC90") ? error("could not find GRID_LLC90/") : nothing
(D,Dexch,Darr,DD)=demo1("LLC90")
(Rini,Rend,DXCsm,DYCsm)=demo2()

include(joinpath(dirname(pathof(MeshArrays)),"gcmfaces_plot.jl"))
qwckplot(Rini)
qwckplot(Rend)
```

## Main Features

**MeshArrays.jl** composite types contain array collections where arrays are typically inter-connected at their edges. It provides `exchange functions` that transfer data between neighbor arrays so that the user can extend the domain of computation as needed e.g. for interpolation or to compute spatial derivatives.

The composite types specify how each array collection forms a mesh and provide information to allow `exchange` to dispatch to the appropriate method. Various configurations that are commonly used in `Earth System Models` are implemented using the concrete type called `gcmfaces`. This type contains a `gcmgrid` instance that provides grid specifications (incl. location on disk). Embedded arrays, or meshes, each represent a subdomain within, e.g., an Earth System Model grid.

The `gcmfaces` name derives from a [previous Matlab / Octave package](https://gcmfaces.readthedocs.io/en/latest/) that was introduced in [Forget et al., 2015](http://www.geosci-model-dev.net/8/3071/2015/), `doi:10.5194/gmd-8-3071-2015`, and inspired this `Julia` package. Here, `GCM` is an acronym for General Circulation Model, or Global Climate Model as described in [this wikipedia entry](https://en.wikipedia.org/wiki/General_circulation_model), and `faces` is just another name for meshes, arrays, facets, or subdomains.

## API index

```@index
```

## API details

```@autodocs
Modules = [MeshArrays]
Order   = [:type,:function]
```
