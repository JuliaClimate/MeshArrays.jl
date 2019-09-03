# MeshArrays.jl


[![Travis Build Status](https://travis-ci.org/gaelforget/MeshArrays.jl.svg?branch=master)](https://travis-ci.org/gaelforget/MeshArrays.jl)
[![codecov](https://codecov.io/gh/gaelforget/GCMFaces.jl/branch/master/graph/badge.svg)](https://codecov.io/gh/gaelforget/GCMFaces.jl)
[![Coverage Status](https://coveralls.io/repos/github/gaelforget/GCMFaces_jl/badge.svg?branch=master)](https://coveralls.io/github/gaelforget/GCMFaces_jl?branch=master)
[![DOI](https://zenodo.org/badge/143987632.svg)](https://zenodo.org/badge/latestdoi/143987632)
[![](https://img.shields.io/badge/docs-stable-blue.svg)](https://gaelforget.github.io/MeshArrays.jl/stable)
[![](https://img.shields.io/badge/docs-dev-blue.svg)](https://gaelforget.github.io/MeshArrays.jl/dev)

`MeshArrays.jl` primarily defines composite types that embed inter-connected array collections within a `struct` and provides an `exchange` function that effectively transfers data between connected arrays. It was originally introduced, as `gcmfaces.jl`, in this [JuliaCon-2018 presentation](https://youtu.be/RDxAy_zSUvg) (see below for **notebooks**). _Note:_ `MeshArrays.jl` is registered, documented, archived, and routinely tested, but is still regarded as a **preliminary implementation**.

### Installation

```
using Pkg
Pkg.add("MeshArrays")
Pkg.test("MeshArrays")
```

### Use example

```
using MeshArrays

isdir("GRID_LLC90") ? GridVariables=GCMGridLoad(GCMGridSpec("LLC90")) : GridVariables=GCMGridOnes("cs",6,100);

(Rini,Rend,DXCsm,DYCsm)= MeshArrays.demo2(GridVariables);

include(joinpath(dirname(pathof(MeshArrays)),"gcmfaces_plot.jl"))
qwckplot(Rini)
qwckplot(Rend)
```

### Notebooks

The [JuliaCon-2018 presentation](https://youtu.be/RDxAy_zSUvg) relied on two `Jupyter notebooks` that are available in the [JuliaCon2018Notebooks repo](https://github.com/gaelforget/JuliaCon2018Notebooks.git). Another notebook is included which illustrate how `MeshArrays.smooth` is used for unit testing purposes (`demo_smooth.ipynb`). These notebooks can readily be executed via the `launch binder` badge in the notebooks repo.
