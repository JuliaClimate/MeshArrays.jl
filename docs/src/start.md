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
show(tmp)
```

## Tutorial

In this tutorial, the same example is repeated three times -- for three different grid configurations commonly used in numerical models (see [Earth Grids](@ref)). The example proceeds as follows in each grid case:

1. generate grid configuration
1. initialize map with random noise
1. apply smoothing across whole domain
1. plot results for the subdomain arrays

Step 3 illustrates how `MeshArrays.jl` computes partial derivatives between neighboring subdomains. The `MeshArrays.smooth()` method is indeed based on lateral diffusion, a transport process, which is integrated through time over the global domain. Other processes like e.g. advection by ocean currents would work similarly.

### 1. Doubly Periodic Domain

Let's start with a doubly periodic domain split into `nF=16` subdomains. Each subdomain corresponds to an array of `nP=20` by `nQ=20` grid points. 

```
(nP,nQ,nF)=(20,20,16)
facesSize=fill((nP,nQ),nF)
ioSize=[nP nQ*nF]
```

The `UnitGrid()` function is used to generate such a grid configuration.

```
using MeshArrays
γ=gcmgrid("","PeriodicDomain",
			nF,facesSize, ioSize,
			Float32, read, write)
Γ=UnitGrid(γ)
```

Then we do steps 2 (`zin`) and 3 (`zout`) as follows.

```
#initialize 2D field of random numbers
tmp1=randn(Float32,Tuple(γ.ioSize))
zin =γ.read(tmp1,MeshArray(γ,Float32))

#smoothing length scales in x, y directions
Lx=3*Γ.DXC; Ly=3*Γ.DYC

#apply smoother
zout=smooth(zin,Lx,Ly,Γ)
```

The `heatmap` function allows you to visualize that `zout` is indeed smoother (and therefore muted) than is `zin`.

```
using Plots
p=dirname(pathof(MeshArrays))
include(joinpath(p,"../examples/Plots.jl"))
heatmap(zout,clims=(-0.25,0.25),tickfont = (4, :black))
```

Grid scale noise           |  Smoothed noise
:------------------------------:|:---------------------------------:
![raw](https://user-images.githubusercontent.com/20276764/118325229-2d883d80-b4d1-11eb-953b-ddbb11bcfe1b.png)  |  ![smooth](https://user-images.githubusercontent.com/20276764/118325093-f31ea080-b4d0-11eb-8c6e-8cd0cc2cc255.png)

### 2. Cube Sphere

Now we instead use a grid that has 6 subdomains, 100x100 points each, covering the six faces of a cube. We note that this _cube sphere_ topology involves connections between subdomain that are slightly more complicated than in the first example.

```
(nP,nQ,nF)=(32,32,6)
facesSize=fill((nP,nQ),nF)
ioSize=[nP nQ*nF]

using MeshArrays
γ=gcmgrid("","CubeSphere",
			nF,facesSize, ioSize,
			Float32, read, write)
Γ=UnitGrid(γ)
```

Next we call `demo2()` which combines steps 2 and 3. 

```
p=dirname(pathof(MeshArrays))
include(joinpath(p,"../examples/Demos.jl"))
D=demo2(Γ)
```

### 3. LLC90 grid

The [Lat-Lon-Cap grid](http://www.geosci-model-dev.net/8/3071/2015/) (or LLC) is a global ocean model grid which is widely used in the [MITgcm user community](https://mitgcm.readthedocs.io/en/latest/). It has 5 uneven subdomains, variable grid spacing, and continents [(Forget et al 2015)](http://www.geosci-model-dev.net/8/3071/2015/). LLC90's resolution is one degree albeit with modications in the Arctic and along the Equator.

In this case, the grid variables are read from files found in the `MeshArrays.GRID_LLC90` folder by the `GridLoad` function. The `GridSpec` function provides the corresponding domain sizes for this commonly used, `LLC90`, grid. 

```
γ=GridSpec("LatLonCap",MeshArrays.GRID_LLC90)
Γ=GridLoad(γ)
D=demo2(Γ)
```

## Earth Grids

Here we visualize a subset of grid lines in a cube sphere (top right), LLC (bottom right), and two other grids.

![EarthGrids](https://raw.githubusercontent.com/gaelforget/MeshArrays.jl/master/docs/images/sphere_all.png)


