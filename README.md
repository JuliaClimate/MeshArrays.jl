# MeshArrays.jl

[![](https://img.shields.io/badge/documentation-blue.svg)](https://juliaclimate.github.io/MeshArrays.jl/dev)

**MeshArrays.jl** defines the `MeshArray` data structure (or type) that can contain collections of inter-connected arrays as generally done in climate models. This provides a simple yet efficient and general way to e.g. analyze climate system simulations and manipulate their output.

Please refer to the [Docs](https://juliaclimate.github.io/MeshArrays.jl/dev/), [Tutorials](https://juliaclimate.github.io/MeshArrays.jl/dev/tutorials/), and [videos](https://juliaclimate.github.io/MeshArrays.jl/dev/videos/) for more information.

[![CI](https://github.com/juliaclimate/MeshArrays.jl/actions/workflows/ci.yml/badge.svg)](https://github.com/juliaclimate/MeshArrays.jl/actions/workflows/ci.yml)
[![codecov](https://codecov.io/gh/juliaclimate/MeshArrays.jl/branch/master/graph/badge.svg)](https://codecov.io/gh/juliaclimate/MeshArrays.jl)
[![Binder](https://mybinder.org/badge_logo.svg)](https://mybinder.org/v2/gh/JuliaClimate/MeshArrays.jl/master)
[![DOI](https://zenodo.org/badge/143987632.svg)](https://zenodo.org/badge/latestdoi/143987632)

Some features and related packages:

- interpolate to geographic coordinates and visualize as maps
- visualize on the sphere, in geographic coordinates, in projection, or grid coordinate
- examples include interaction and animation using `Makie.jl` and `Pluto.jl`
- accurately derive planetary scale transports on a global [climate model grid](https://en.wikipedia.org/wiki/Arakawa_grids)
- compute trajectories of materials over any supported grid via [IndividualDisplacements.jl](https://github.com/JuliaClimate/IndividualDisplacements.jl)
- read/write [NetCDF](https://en.wikipedia.org/wiki/NetCDF), files and other formats often used in climate sciences
- support domain decomposition via [NCTiles.jl](https://gaelforget.github.io/NCTiles.jl/stable/) that 
- support for [MITgcm](https://mitgcm.readthedocs.io/en/latest/) model ouput via [MITgcm.jl](https://github.com/gaelforget/MITgcm.jl)
- support for [ECCO](https://doi.org/10.5194/gmd-8-3071-2015) ocean estimates via [OceanStateEstimation.jl](https://github.com/gaelforget/OceanStateEstimation.jl)

`MeshArrays.jl` was first introduced in [this presentation](https://youtu.be/RDxAy_zSUvg) at JuliaCon 2018.

| | | |
|:-------------------------------------:|:-------------------------------------:|:-------------------------------------:|
| <img src="https://user-images.githubusercontent.com/20276764/144332405-ed8d163f-04b9-408a-8fd0-08d91e9be91b.png" width="180" height="120"> | <img src="https://user-images.githubusercontent.com/20276764/144878637-1412679c-f1e6-4491-a8f1-43d729aa224d.png" width="150" height="150"> | <img src="https://user-images.githubusercontent.com/20276764/215533819-d0fe6709-6040-4a71-ad50-cfd5c43e6030.png" width="140" height="120">


