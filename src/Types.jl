
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

## concrete types and MeshArray alias:

include("gcmfaces_type.jl");
include("gcmarray_type.jl");
include("gcmvector_type.jl");

#MeshArray=gcmfaces
MeshArray=gcmarray

## Methods that apply to all AbstractMeshArray types

import Base: maximum, minimum, sum, fill, fill!

function maximum(a::AbstractMeshArray)
  c=-Inf;
  for I in eachindex(a)
    c = max(c,maximum(a[I]))
  end
  return c
end

function minimum(a::AbstractMeshArray)
  c=-Inf;
  for I in eachindex(a)
    c = min(c,minimum(a[I]))
  end
  return c
end

function sum(a::AbstractMeshArray)
  c=0.0
  for I in eachindex(a)
    c = c + sum(a[I])
  end
  return c
end

function fill(val::Any,a::AbstractMeshArray)
  c=similar(a)
  for I in eachindex(a.f)
    c.f[I] = fill(val,size(a.f[I]))
  end
  return c
end

function fill!(a::AbstractMeshArray,val::Any)
  for I in eachindex(a.f)
    fill!(a.f[I],val)
  end
  return a
end

import Base: +, -, *, /

function +(a::AbstractMeshArray,b::AbstractMeshArray)
  c=similar(a)
  for I in eachindex(a.f)
    c.f[I] = a.f[I] + b.f[I]
  end
  return c
end

function -(a::AbstractMeshArray,b::AbstractMeshArray)
  c=similar(a)
  for I in eachindex(a.f)
    c.f[I] = a.f[I] - b.f[I]
  end
  return c
end

function /(a::AbstractMeshArray,b::AbstractMeshArray)
  c=similar(a)
  for I in eachindex(a.f)
    c.f[I] = a.f[I] ./ b.f[I]
  end
  return c
end

function *(a::AbstractMeshArray,b::AbstractMeshArray)
  c=similar(a)
  for I in eachindex(a.f)
    c.f[I] = a.f[I] .* b.f[I]
  end
  return c
end
