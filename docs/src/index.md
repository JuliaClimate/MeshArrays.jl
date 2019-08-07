# MeshArrays.jl documentation

`MeshArrays.jl` primarily defines composite types that embed inter-connected array collections within a `struct` and provides an `exchange` function that effectively transfers data between connected arrays. It was originally introduced, as `gcmfaces.jl`, in this [JuliaCon-2018 presentation](https://youtu.be/RDxAy_zSUvg) (see below for **notebooks**). _Note: even though `MeshArrays.jl` is registered, documented, archived, and routinely tested, it is still regarded as a **preliminary implementation**._

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

The [JuliaCon-2018 presentation](https://youtu.be/RDxAy_zSUvg) relied on two `Jupyter notebooks` privided in the [JuliaCon2018Notebooks repo](https://github.com/gaelforget/JuliaCon2018Notebooks.git) and pre-defined grids that that can be downloaded as follows:

```
git clone https://github.com/gaelforget/GRID_CS32
git clone https://github.com/gaelforget/GRID_LLC90
```

At the command line, one can reproduce the same computation as follows:

```
using MeshArrays
isdir("GRID_LLC90") ? GridVariables=GCMGridLoad(GCMGridSpec("LLC90")) : GridVariables=GCMGridOnes("cs",6,100);
                    
(Rini,Rend,DXCsm,DYCsm)= MeshArrays.demo2(GridVariables);
```

For plotting directions, see the help section (`?MeshArrays.demo2`). Additional demos are also provided (e.g., see `?MeshArrays.demo1`, `?MeshArrays.demo3`).

## Main Features

**MeshArrays.jl** composite types contain array collections where arrays are typically inter-connected at their edges. It provides `exchange functions` that transfer data between neighbor arrays so that the user can extend the domain of computation as needed e.g. for interpolation or to compute spatial derivatives.

The composite types specify how each array collection forms a mesh and provide information to allow `exchange` to dispatch to the appropriate method. Various configurations that are commonly used in `Earth System Models` are implemented using the concrete type called `gcmfaces`. This type contains a `gcmgrid` instance that provides grid specifications (incl. location on disk). Embedded arrays, or meshes, each represent a subdomain within, e.g., an Earth System Model grid.

_About names:_ the `gcmfaces` name derives from a [previous Matlab / Octave package](https://gcmfaces.readthedocs.io/en/latest/), introduced in [Forget et al., 2015](http://www.geosci-model-dev.net/8/3071/2015/) (`doi:10.5194/gmd-8-3071-2015`), which provided the inspiratoin for `MeshArrays.jl`. To break it down, `GCM` is an acronym for General Circulation Model, or Global Climate Model (see [this wikipedia entry](https://en.wikipedia.org/wiki/General_circulation_model)), and `faces` can mean meshes, arrays, facets, or subdomains.

## API index

```@index
```

## API details

```@autodocs
Modules = [MeshArrays]
Order   = [:type,:function]
```
