# Get Started

## Install

```
using Pkg
Pkg.add("MeshArrays")
Pkg.test("MeshArrays")
```

!!! note
    The `Julia` package manager, `Pkg.jl`, is documented [here](https://docs.julialang.org/en/v1/) and [further here](https://julialang.github.io/Pkg.jl/v1/).

## Tutorial

The sequence of examples in this tutorial is as follows :

1. generate a grid configuration
1. create a map of random noise
1. apply diffusion (a transport process)
1. plot results for each subdomain array

These rely on grid configurations commonly used in global models (see [Earth Grids](@ref)). In step 3, the smoothing is based on integrating a lateral diffusion equation through time over the global domain. This illustrates how `MeshArrays.jl` computes partial derivatives between neighboring subdomains. 

### 1. Doubly Periodic Domain

Let's setup a doubly periodic domain with 16 subdomains. Each one contains an array of 20 by 20 grid points.

```
using MeshArrays; p=dirname(pathof(MeshArrays))
γ,Γ=GridOfOnes("PeriodicDomain",16,20)

include(joinpath(p,"../examples/Demos.jl"))
(xi,xo,_,_)=demo2(Γ);
show(xo)

using Plots
include(joinpath(p,"../examples/Plots.jl"))
heatmap(xo,clims=(-0.25,0.25),colorbar=false,tickfont = (4, :black))
```

Grid scale noise           |  Smoothed noise
:------------------------------:|:---------------------------------:
![raw](https://user-images.githubusercontent.com/20276764/118325229-2d883d80-b4d1-11eb-953b-ddbb11bcfe1b.png)  |  ![smooth](https://user-images.githubusercontent.com/20276764/118325093-f31ea080-b4d0-11eb-8c6e-8cd0cc2cc255.png)

### 2. Cube Sphere

This grid has 6 subdomains, 100x100 points each, covering the six faces of a cube.

```
γ,Γ=GridOfOnes("CubeSphere",6,100)
D=demo2(Γ)
```

### 3. LLC90 grid

The [Lat-Lon-Cap grid](http://www.geosci-model-dev.net/8/3071/2015/) (or LLC) is a global ocean model grid which is widely used in the [MITgcm user community](https://mitgcm.readthedocs.io/en/latest/). It has 5 uneven subdomains, variable grid spacing, and continents [(Forget et al 2015)](http://www.geosci-model-dev.net/8/3071/2015/). LLC90's resolution is one degree albeit with modications in the Arctic and along the Equator.

```
#run(`git clone https://github.com/gaelforget/GRID_LLC90`)
Γ=GridLoad(GridSpec("LatLonCap",MeshArrays.GRID_LLC90))
D=demo2(Γ)
heatmap(D[2],clims=(-0.25,0.25))
```

## Earth Grids

Here we visualize a subset of grid lines in a cube sphere (top right), LLC (bottom right), and two other grids.

![EarthGrids](https://raw.githubusercontent.com/gaelforget/MeshArrays.jl/master/docs/images/sphere_all.png)


