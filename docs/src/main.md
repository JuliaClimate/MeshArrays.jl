## Main Features

The [JuliaCon-2018 presentation](https://youtu.be/RDxAy_zSUvg) relied on two `Jupyter notebooks` available in [this other repo](https://github.com/gaelforget/JuliaCon2018Notebooks.git) which correspond to `MeshArrays.demo1` and `MeshArrays.demo2`. Additional demo functions are provided within `MeshArrays` and as notebooks (e.g., `demo3` & `demo_trsp.ipynb`).

`MeshArrays.jl` composite array types have array elements. The elementary arrays within a `MeshArray` are typically inter-connected at their edges according to a user-specified `gcmgrid` specification. `exchange` methods transfer data between neighboring arrays to extend their computational domains, as often needed to derive e.g. planetary transports in the climate system.

The current default for `MeshArray` is the `gcmarray` type and an instance `H` is shown below. This example is based on a grid known as `LLC90` where each global map is associated with 5 subdomains. Hence, `H.f` is a `5x50` array when `H` represents a gridded variable on `50` depth levels, and the elements of  `H.f` are of size `(90, 270)`, `(90, 90)`, or `(270, 90)`. 

```
julia> show(H)
 gcmarray 
  grid type   = llc
  array size  = (5, 50)
  data type   = Float64
  face sizes  = (90, 270)
                (90, 270)
                (90, 90)
                (270, 90)
                (270, 90)
```

The underlying, `MeshArray`, data structure is:

```
struct gcmarray{T, N} <: AbstractMeshArray{T, N}
   grid::gcmgrid
   f::Array{Array{T,2},N}
   fSize::Array{NTuple{2, Int}}
   fIndex::Array{Int,1}
end
```

And the embedded `gcmgrid` specification is constructed as: 

```
gcmgrid(path::String, class::String, nFaces::Int,
        fSize::Array{NTuple{2, Int},1}, ioSize::Array{Int64,2},
        ioPrec::Type, read::Function, write::Function)
```

