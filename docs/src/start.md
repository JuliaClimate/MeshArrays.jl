# Get Started

## Install

To install `MeshArrays.jl` and verify it works as expected, open the `Julia` REPL and type:

```
using Pkg
Pkg.add("MeshArrays")
Pkg.test("MeshArrays")
```

## Use

To create your first `MeshArray`, open the `Julia` REPL and type:

```
using MeshArrays
tmp=MeshArray(randn(20,10))
```

## Grids

Below we visualize a subset of grid lines in a cube sphere (top right), LLC (bottom right), and two other grids. 

Three such grids are available directly via this package 
 the examples (`GRID_LL360`, `GRID_CS32`, and `GRID_LLC90`).

![EarthGrids](https://raw.githubusercontent.com/gaelforget/MeshArrays.jl/master/docs/images/sphere_all.png)


