
using Unitful, Dates

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
        fSize::Array{NTuple{2, Int},1},
        ioSize::Union{NTuple{2, Int},Array{Int64,2}},
        ioPrec::Type, read::Function, write::Function)
```

The `class` can be set to "LatLonCap", "CubeSphere", "PeriodicChannel", "PeriodicDomain".

For example, A periodic channel (periodic in the x direction) of size 360 by 160, can be defined as follows.

```
pth=MeshArrays.GRID_LL360
class="PeriodicChannel"
ioSize=(360, 160)
ioPrec=Float32

γ=gcmgrid(pth, class, 1, [ioSize], ioSize, ioPrec, read, write)

Γ=GridLoad(γ)    
```

Please refer to `GridSpec` and `UnitGrid` for more info related to `gcmgrid` options.
"""
struct gcmgrid
  path::String
  class::String
  nFaces::Int
  fSize::Array{NTuple{2, Int},1}
  ioSize::Union{NTuple{2, Int},Array{Int64,2}}
  ioPrec::Type
  read::Function
  write::Function
end

copy_if_isarray(x) = isa(x,Array) ? copy(x) : x

Base.similar(g::gcmgrid)=gcmgrid(g.path, g.class, g.nFaces, copy(g.fSize), copy_if_isarray(g.ioSize), g.ioPrec, g.read, g.write)

"""
    varmeta

varmeta data structure. By default, `unit` is `missing` (non-dimensional), `position`
is `fill(0.5,3)` (cell center), time is `missing`, and `name` / `long_name` is `unknown`.

Available constructors:

```
varmeta(unit::Union{Unitful.Units,Number,Missing},position::Array{Float64,1},
        time::Union{DateTime,Missing,Array{DateTime,1}},name::String,long_name::String)
```

And:

```defaultmeta = varmeta(missing,fill(0.5,3),missing,"unknown","unknown")```

"""
struct varmeta
  unit::Union{Unitful.Units,Number,Missing}
  position::Array{Float64,1}
  time::Union{DateTime,Missing,Array{DateTime,1}}
  name::String
  long_name::String
end

defaultmeta = varmeta(missing,fill(0.5,3),missing,"unknown","unknown")

## concrete types and MeshArray alias:

include("Type_gcmfaces.jl");
OuterArray{T,N}=Array{T,N} where {T,N}
InnerArray{T,N}=Array{T,N} where {T,N}
include("Type_gcmarray.jl");
include("Type_gcmvector.jl");

"""
    MeshArray

Alias to `gcmarray` or `gcmfaces` concrete type
"""
MeshArray=gcmarray

## Methods that apply to all AbstractMeshArray types

import Base: maximum, minimum, extrema, sum, fill, fill!

extrema(a::AbstractMeshArray) = [minimum(a) maximum(a)]

function maximum(a::AbstractMeshArray)
  c=-Inf;
  for I in eachindex(a)
    c = max(c,maximum(a[I]))
  end
  return c
end

function minimum(a::AbstractMeshArray)
  c=Inf;
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

fill(val::Any,a::gcmgrid,args...) = val .* ones(a,args...)

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

function *(a::AbstractMeshArray,b::Array)
  c=MeshArray(a.grid,eltype(a[eachindex(a.f)[1]]),size(b)...)
  for f in 1:length(a.f)
    for bb in eachindex(IndexCartesian(),b)
      c.f[f,bb.I...].=b[bb.I...]*a[f]
    end
  end
  return c
end

import Base: ones, zeros

function zeros(a::gcmgrid,args...)
  b=MeshArray(a)
  [b.f[c].=0.0 for c in eachindex(b.f)]
  (length(args)>0 ? b*ones(args...) : b) 
end

function ones(a::gcmgrid,args...)
  b=MeshArray(a)
  [b.f[c].=1.0 for c in eachindex(b.f)]
  (length(args)>0 ? b*ones(args...) : b) 
end

function zeros(a::AbstractMeshArray)
  b=similar(a)
  [b.f[c].=0.0 for c in eachindex(b.f)]
  b
end

function ones(a::AbstractMeshArray)
  b=similar(a)
  [b.f[c].=1.0 for c in eachindex(b.f)]
  b
end

## derivative types


"""
    gridpath

gridpath data structure.

```
using MeshArrays
γ=GridSpec("LatLonCap",MeshArrays.GRID_LLC90)
Γ=GridLoad(γ;option="light")
lons=[-68 -63]; lats=[-54 -66]; name="Drake Passage"
Trsct=Transect(name,lons,lats,Γ)
```
"""
Base.@kwdef struct gridpath
#  grid::gcmgrid
  grid::NamedTuple
  name::String
  C::Matrix
  W::Matrix
  S::Matrix
end

"""
    gridmask

gridmask data structure.
"""
Base.@kwdef struct gridmask
  map::MeshArray
  names::Array
  depths::Array
  h_sum::Array
  v_int::Array
  tmp2d::MeshArray
  tmp3d::MeshArray
end



