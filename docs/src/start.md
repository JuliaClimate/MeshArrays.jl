# Get Started

To create your first `MeshArray`, open the `Julia` REPL and type:

```@example 1
using MeshArrays
C=MeshArray(randn(21,10))
```

```@example 1
using CairoMakie
heatmap(C,x=-10:10,y=1:10)
```

## Grids

Three grids are available directly via this package the examples (`GRID_LL360`, `GRID_CS32`, and `GRID_LLC90`).


| Lat-lon | Cube Sphere | Lat-Lon-Cap |
|:-------------------------------------:|:-------------------------------------:|:-------------------------------------:|
![Lat-Lon](https://user-images.githubusercontent.com/20276764/144249858-df986169-8f4a-4c42-bf64-45bc97c34ca8.png) | ![Cube Sphere](https://user-images.githubusercontent.com/20276764/144249876-a37ba2da-7258-4f01-b438-0b3efbf75c2d.png) | ![Lat-Lon-Cap](https://user-images.githubusercontent.com/20276764/144249899-4d94980a-87aa-4bfb-a6d6-6145f9f0324f.png)

In the example below we read and display `GRID_CS32`. This is a [cubed-sphere grid](https://en.wikipedia.org/wiki/Quadrilateralized_spherical_cube) where the EArth surface is projected on the six faces of a cube.

```@example
using MeshArrays,CairoMakie
γ=GridSpec("CubeSphere",MeshArrays.GRID_CS32)
Γ=GridLoad(γ)
heatmap(Γ.YC,title="grid point latitudes")
```

## Install

To install `MeshArrays.jl` and verify it works as expected, open the `Julia` REPL and type:

```
using Pkg
Pkg.add("MeshArrays")
Pkg.test("MeshArrays")
```
