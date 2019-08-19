
"""
    AbstractMeshArray{T, N}

Subtype of AbstractArray{T, N}
"""
abstract type AbstractMeshArray{T, N} <: AbstractArray{T, N} end

"""
    gcmgrid

gcmgrid data structure. Available constructors:

```
gcmgrid(path::String, class::String, nFaces::Int,
        fSize::Array{NTuple{2, Int},1}, ioSize::Array{Int64,2},
        ioPrec::Type, read::Function, write::Function)
```
"""
struct gcmgrid
  path::String
  class::String
  nFaces::Int
  fSize::Array{NTuple{2, Int},1}
#  ioSize::NTuple{2, Int}
  ioSize::Array{Int64,2}
  ioPrec::Type
  read::Function
  write::Function
end

## concrete type and MeshArray alias:

#include("gcmfaces_type.jl"); MeshArray=gcmfaces
#include("gcmarray_type.jl"); MeshArray=gcmarray
#include("gcmarray_type01.jl"); MeshArray=gcmarray
include("gcmarray_type02.jl"); MeshArray=gcmarray
