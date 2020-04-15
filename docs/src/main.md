## Main Features

The elements of a `MeshArray` are arrays. These elementary arrays typically represent subdomains inter-connected at their edges. The organization and connections between subdomains is determined by a user-specified `gcmgrid` which is embeded inside each `MeshArray` instance. 

`InterpolationFactors` ... `Exchange` methods transfer data between neighboring arrays to extend computational subdomains -- this is often needed in analyses of climate or ocean model output. 

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

### Examples, Interpolation, & Plots

The [JuliaCon-2018 presentation](https://youtu.be/RDxAy_zSUvg) relied on two `Jupyter` notebooks available in [GlobalOceanNotebooks](https://github.com/juliaclimate/GlobalOceanNotebooks.git)/DataStructures which are also included here in `examples/Demos.jl`. [GlobalOceanNotebooks](https://github.com/juliaclimate/GlobalOceanNotebooks.git)/OceanTransports provides use case examples related to Earth System transports.


A simple way to plot a `MeshArray` consists in plotting each elementary array separately. Other methods that e.g. produce global maps and projections are illustrated in the notebooks. A simple one is shown below that demonstrates the included interpolation scheme.

```
p=dirname(pathof(MeshArrays));
using Plots; include(joinpath(p,"../examples/Plots.jl"));
heatmap(D,title="Ocean Depth",clims=(0.,6000.))
```

![OceanDepthMap](https://raw.githubusercontent.com/juliaclimate/MeshArrays.jl/master/docs/images/ocean_depth.png)

```
lon=[i for i=-179.5:1.0:179.5, j=-89.5:1.0:89.5]
lat=[j for i=-179.5:1.0:179.5, j=-89.5:1.0:89.5]

Γ=GridLoad(GridSpec("LatLonCap","GRID_LLC90/"))
(f,i,j,w)=InterpolationFactors(Γ,vec(lon),vec(lat))
DD=Interpolate(Γ["Depth"],f,i,j,w)

contourf(vec(lon[:,1]),vec(lat[1,:]),DD,clims=(0.,6000.))
```

![OceanDepthMap](https://raw.githubusercontent.com/juliaclimate/MeshArrays.jl/master/docs/images/interp_depth.png)
