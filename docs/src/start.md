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

[The extended tutorial](basics.html) further illustrate how `MeshArrays.jl` let's the user write generic code that is readily applicable to whole grid families. 

It focuses on a global, smoothing, workflow that requires communication across the entire gridded domain -- a key feature provided by `MeshArrays.jl`.

The same workflow is repeated three times, for different grid configurations commonly used in numerical models (see Earth Grids below). 

Grid scale noise           |  Smoothed noise
:------------------------------:|:---------------------------------:
![raw](https://user-images.githubusercontent.com/20276764/118325229-2d883d80-b4d1-11eb-953b-ddbb11bcfe1b.png)  |  ![smooth](https://user-images.githubusercontent.com/20276764/118325093-f31ea080-b4d0-11eb-8c6e-8cd0cc2cc255.png)

## Earth Grids

Here we visualize a subset of grid lines in a cube sphere (top right), LLC (bottom right), and two other grids.

![EarthGrids](https://raw.githubusercontent.com/gaelforget/MeshArrays.jl/master/docs/images/sphere_all.png)


