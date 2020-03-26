## Main Features

An instance of the `MeshArray` type has arrays as elements. Elementary arrays within a `MeshArray` may represent subdomains inter-connected at their edges. The organization and connections between subdomains is determined by a user-specified `gcmgrid` which is embeded in `MeshArray` instances. `Exchange` methods transfer data between neighboring arrays to extend computational subdomains -- this is often needed in analyses of climate or ocean model output.

The current default for `MeshArray` is the `gcmarray` type and an instance `H` is shown below. This example is based on a grid known as `LatLonCap` where each global map is associated with 5 subdomains. Hence, `H.f` is a `(5, 50)` array when `H` represents a gridded variable on `50` depth levels, and elements of  `H.f` are arrays of size `(90, 270)`, `(90, 90)`, or `(270, 90)`. 

```
julia> show(D)
  name        = Depth
  unit        = m
  data type   = Float64
  cell pos.   = [0.5, 0.5]

  tile array  = (5,)
  tile sizes  = (90, 270)
                (90, 270)
                (90, 90)
                (270, 90)
                (270, 90)

  grid class  = LatLonCap
  MeshArray   = gcmarray 
  version     = 0.2.7 
```

The underlying, `MeshArray`, data structure is:

```
struct gcmarray{T, N} <: AbstractMeshArray{T, N}
   grid::gcmgrid
   meta::varmeta
   f::Array{Array{T,2},N}
   fSize::Array{NTuple{2, Int}}
   fIndex::Array{Int,1}
   version::String
end
```

A `MeshArray` generally behaves just like an `Array` including for operations listed below. The _broadcasting_ function has been customized so that it reaches elements of each elementary array.

```
size(D)
eltype(D)
view(D,:)

D .* 1.0
D .* D
1000*D
D*1000

D[findall(D .> 300.)] .= NaN
D[findall(D .< 1.)] .= NaN

D[1]=0.0 .+ D[1]
tmp=cos.(D)
```

In addition, `Mesharray` specific functions like `exchange` can alter the internal structure of a `MeshArray`. Elementary array sizes are thus larger in `show(exchange(D))` than `show(D)`.

```
julia> show(exchange(D))
  tile sizes  = (92, 272)
                (92, 272)
                (92, 92)
                (272, 92)
                (272, 92)
```

### Embedded Meta Data

A `MeshArray` includes a `gcmgrid` specification which can be constructed as outlined below.

```
gcmgrid(path::String, class::String, 
        nFaces::Int, fSize::Array{NTuple{2, Int},1}, 
        ioSize::Array{Int64,2}, ioPrec::Type, 
        read::Function, write::Function)
```

Importantly, a `gcmgrid` does **not** contain any actual grid data -- hence its memory footprint is minimal. Grid variables are instead read to memory only when needed e.g. as shown below. To make this easy, each `gcmgrid` includes a pair of `read` / `write` methods to allow for basic `I/O` at any time. These methods are typically specified by the user although defaults are provided. 

```
using MeshArrays, Unitful
γ=GridSpec("LatLonCap","GRID_LLC90/")
m=MeshArrays.varmeta(u"m",fill(0.5,2),"Depth","Depth")
D=γ.read(γ.path*"Depth.data",MeshArray(γ,Float64;meta=m))
```

The above commands define a `MeshArray` called `D` which is the one displayed at the top of this section. A definition of the `varmeta` structure is reported below. The `position` of a `D` point within its grid cell is given as `x ∈ [0. 1.]` in each direction.

```
varmeta(unit::Union{Unitful.AbstractQuantity,Number},
        position::Array{Float64,1},
        name::String,long_name::String)
```

### Examples

The [JuliaCon-2018 presentation](https://youtu.be/RDxAy_zSUvg) relied on two `Jupyter` notebooks that are available in the [MeshArrayNotebooks repo](https://github.com/gaelforget/JuliaCon2018Notebooks.git). `demo1` and `demo2` are very similar. 

Standard oceanography examples are also provided in [MeshArrayNotebooks](https://github.com/gaelforget/JuliaCon2018Notebooks.git) (e.g., `04_transports.ipynb`) and / or in `MeshArrays.jl` (e.g., `demo3`).

A simple way to plot a `MeshArray` consists in plotting each elementary array separately (see below). Other methods that e.g. produce global maps and projections are illustrated in the notebooks. 

```
p=dirname(pathof(MeshArrays));
using Plots; include(joinpath(p,"Plots.jl"));
heatmap(D,title="Ocean Depth",clims=(0.,6000.))
```

![OceanDepthMap](https://raw.githubusercontent.com/gaelforget/MeshArrays.jl/master/docs/images/ocean_depth.png)
