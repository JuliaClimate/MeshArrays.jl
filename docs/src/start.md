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

Three grids are available directly via this package the examples (`GRID_LL360`, `GRID_CS32`, and `GRID_LLC90`).

```
using MeshArrays

#γ=GridSpec("PeriodicChannel",MeshArrays.GRID_LL360)
#γ=GridSpec("CubeSphere",MeshArrays.GRID_CS32)
γ=GridSpec("LatLonCap",MeshArrays.GRID_LLC90)

Γ=GridLoad(γ);
```

| Lat-lon | Cube Sphere | Lat-Lon-Cap |
|:-------------------------------------:|:-------------------------------------:|:-------------------------------------:|
![Lat-Lon](https://user-images.githubusercontent.com/20276764/144249858-df986169-8f4a-4c42-bf64-45bc97c34ca8.png) | ![Cube Sphere](https://user-images.githubusercontent.com/20276764/144249876-a37ba2da-7258-4f01-b438-0b3efbf75c2d.png) | ![Lat-Lon-Cap](https://user-images.githubusercontent.com/20276764/144249899-4d94980a-87aa-4bfb-a6d6-6145f9f0324f.png)

