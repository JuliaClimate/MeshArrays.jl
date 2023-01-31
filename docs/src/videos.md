# Video Examples

The videos below introduce `MeshArrays.jl` and present two use cases :

1. introduction of `MeshArrays.jl` at [JuliaCon 2018](https://www.youtube.com/live/RDxAy_zSUvg?feature=share), for context.
1. interactive analysis of gridded ocean model output for [ECCO](https://ecco-group.org/storymaps.htm?id=69).
1. simulation of plastic particles following [ocean currents](https://youtu.be/6pvKW1hcghg).

## JuliaCon 2018

`MeshArrays.jl` was first introduced at [JuliaCon in 2018](https://juliacon.org/2018/). This video provides some motivation and context for the project.

[![JuliaCon-2018 presentation](https://user-images.githubusercontent.com/20276764/215771788-c52feaae-1257-4525-aa7e-1ccdc175df30.png)](https://youtu.be/RDxAy_zSUvg)

## Interactive Model Analysis

In this demo of interactive visualization of ocean variables, the notebook provides various options for choosing variables and viewing them. Plots etc react to these user selections as illustrated in the video. Code is based on  MeshArrays.jl (grids, arrays), [OceanStateEstimation.jl](https://github.com/gaelforget/OceanStateEstimation.jl) (velocity fields, etc), [Pluto.jl](https://github.com/fonsp/Pluto.jl#readme) (notebook), [Makie.jl](https://docs.makie.org/stable/) (plotting), and other open source packages.

[![ocean state estimate analysis](https://user-images.githubusercontent.com/20276764/144332405-ed8d163f-04b9-408a-8fd0-08d91e9be91b.png)](https://youtu.be/UEmBnzspSRg)

Another method for interacting with data sets is to use `GLMakie.jl` or `WGLMakie.jl`. For example, [Tyler.jl](https://github.com/MakieOrg/Tyler.jl) can be used to explore high-resolution model output.

[![high resolution model analysis](https://user-images.githubusercontent.com/20276764/215533819-d0fe6709-6040-4a71-ad50-cfd5c43e6030.png)](https://youtu.be/TftqT7oZ0Bs)

## Particle Tracking and Modeling

Here we visualize a simulation of particles moving at a fixed depth in the Ocean (300m depth). This uses [IndividualDisplacements.jl](https://github.com/JuliaClimate/IndividualDisplacements.jl), and `MeshArrays.jl` underneath, to simulate particle trajectories. 

[![simulated particle movie (300m)](https://user-images.githubusercontent.com/20276764/84767001-b89a4400-af9f-11ea-956f-2e207f892c4f.png)](https://youtu.be/M6vAUtIsIIY)

More examples like this, using Julia to run models, and related work in [JuliaOcean](https://github.com/JuliaOcean) are available in the longer video below. 

[![Modeling Marine Ecosystems](https://user-images.githubusercontent.com/20276764/132381907-1ab7d682-ea3d-4db7-b245-3cdb9dd2dcd3.png)](https://www.youtube.com/watch?v=UCIRrXz2ZS0)