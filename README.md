# MeshArrays.jl


[![Travis Build Status](https://travis-ci.org/gaelforget/MeshArrays.jl.svg?branch=master)](https://travis-ci.org/gaelforget/MeshArrays.jl)
[![codecov](https://codecov.io/gh/gaelforget/GCMFaces.jl/branch/master/graph/badge.svg)](https://codecov.io/gh/gaelforget/GCMFaces.jl)
[![Coverage Status](https://coveralls.io/repos/github/gaelforget/GCMFaces_jl/badge.svg?branch=master)](https://coveralls.io/github/gaelforget/GCMFaces_jl?branch=master)

[![DOI](https://zenodo.org/badge/143987632.svg)](https://zenodo.org/badge/latestdoi/143987632)
[![](https://img.shields.io/badge/docs-stable-blue.svg)](https://gaelforget.github.io/MeshArrays.jl/stable)
[![](https://img.shields.io/badge/docs-dev-blue.svg)](https://gaelforget.github.io/MeshArrays.jl/dev)

`MeshArrays.jl` primarily defines composite types that embed inter-connected array collections and provides `exchange` functions that transfer data between connected arrays. It was originally introduced in this [JuliaCon-2018 presentation](https://youtu.be/RDxAy_zSUvg) (see **notebooks** below). 

_Note:_ `MeshArrays.jl` is registered, documented, etc., but still regarded as a **preliminary implementation**.

### Installation

```
using Pkg
Pkg.add("MeshArrays")
Pkg.test("MeshArrays")
```

### Use example

Let's integrate a diffusion equation over the surface of a cube:

```
using MeshArrays; p=dirname(pathof(MeshArrays));

GridVariables=GridOfOnes("cs",6,100)
DemoVariables=MeshArrays.demo2(GridVariables)

using Plots; include(joinpath(p,"Plots.jl"));
heatmap(DemoVariables[2],clims=(-0.5,0.5))
```

Starting from a noisy `DemoVariables[1]`, this leads to a smoothed `DemoVariables[2]`.

### Notebooks

The [JuliaCon-2018 presentation](https://youtu.be/RDxAy_zSUvg) relied on two `Jupyter notebooks` that are available in the [JuliaCon2018Notebooks repo](https://github.com/gaelforget/JuliaCon2018Notebooks.git). Another notebook is included which illustrate how `MeshArrays.smooth` is used for unit testing purposes (`demo_smooth.ipynb`). These notebooks can readily be executed via the `launch binder` badge in the notebooks repo.
