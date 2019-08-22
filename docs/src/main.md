## Main Features

**MeshArrays'** composite types contain array collections. Elementary arrays in a `MeshArray` are typically inter-connected at their edges according to a user-specified global structure. `exchange` methods transfer data between neighboring arrays to extend their computational domains, as often needed to compute e.g. planetary transports in the climate system.

The **current default** for `MeshArray` is `gcmfaces`. Instead, a `gcmarray` instance `H` is depicted below. `H.f` is a `5x50` array of elementary arrays which can be of size `(90, 270)`, `(90, 90)`, or `(270, 90)`. This grid known as `LLC90` is such that each global map is associated with 5 subdomain arrays like this. 

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

The corresponding data structure is:

```
struct gcmarray{T, N} <: AbstractMeshArray{T, N}
   grid::gcmgrid
   f::Array{Array{T,2},N}
   fSize::Array{NTuple{2, Int}}
   fIndex::Array{Int,1}
end
```

And the embedded grid specification is: 

```
gcmgrid(path::String, class::String, nFaces::Int,
        fSize::Array{NTuple{2, Int},1}, ioSize::Array{Int64,2},
        ioPrec::Type, read::Function, write::Function)
```

