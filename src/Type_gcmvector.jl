
using CatViews

"""
    gcmvector{T, N}

gcmvector data structure that can be used for
  subsetting and indexing into a gcmarray.

```
gcmvector{T,N}(grid::gcmgrid,f::Array{Array{T,1},N},
         fSize::Array{NTuple{N, Int}},fIndex::Array{Int,1})

```
"""
struct gcmvector{T,N} <: AbstractMeshArray{T,N}
   grid::gcmgrid
   f::Array{Array{T,1},N}
   fSize::Array{NTuple{2, Int}}
   fIndex::Array{Int,1}
end

## Customize standard functions

import Base: findall, size, show, getindex, setindex!, copyto!

Base.size(A::gcmvector) = size(A.f)
Base.size(A::gcmvector, dim::Integer) = size(A)[dim]

function Base.getindex(A::gcmvector{T,N}, I::Vararg{Union{Int,Array{Int},AbstractUnitRange,Colon}, N}) where {T,N}
  return A.f[I...]
end

function Base.setindex!(A::gcmvector{T,N}, v, I::Vararg{Int, N}) where {T,N}
  return (A.f[I...] = v)
end

function Base.show(io::IO, z::gcmvector{T,N}) where {T,N}
    printstyled(io, " gcmvector \n",color=:normal)
    printstyled(io, "  grid type   = ",color=:normal)
    printstyled(io, "$(z.grid.class)\n",color=:blue)
    printstyled(io, "  data type   = ",color=:normal)
    printstyled(io, "$(eltype(z))\n",color=:blue)
    printstyled(io, "  tile array  = ",color=:normal)
    printstyled(io, "$(size(z))\n",color=:blue)
    printstyled(io, "  tile sizes  = ",color=:normal)
    printstyled(io, "$(size(z.f[1]))\n",color=:blue)
    for iFace=2:length(z.fIndex)
      printstyled(io, "                ",color=:normal)
      printstyled(io, "$(size(z.f[iFace]))\n",color=:blue)
    end
  return
end

function Base.view(A::gcmvector{T, N}, I::Vararg{Union{Int,AbstractUnitRange,Colon}, N}) where {T,N}
  tmpf=view(A.f,I...)
  return gcmvector{eltype(A),ndims(tmpf)}(A.grid,tmpf,A.fSize,A.fIndex)
end

## Customize broadcasting

Base.BroadcastStyle(::Type{<:gcmvector}) = Broadcast.ArrayStyle{gcmvector}()

function Base.similar(bc::Broadcast.Broadcasted{Broadcast.ArrayStyle{gcmvector}}, ::Type{ElType}) where ElType
  # Scan the inputs for the gcmarray:
  A = find_gcmvector(bc)
  # Create the gcmvector output:
  return gcmvector{ElType,ndims(A)}(A.grid,similar(A.f),A.fSize,A.fIndex)
end

find_gcmvector(bc::Base.Broadcast.Broadcasted) = find_gcmvector(bc.args)
find_gcmvector(args::Tuple) = find_gcmvector(find_gcmvector(args[1]), Base.tail(args))
find_gcmvector(x) = x
find_gcmvector(a::gcmvector, rest) = a
find_gcmvector(::Any, rest) = find_gcmvector(rest)

# Specialize this method if all you want to do is specialize on typeof(dest)
@inline function copyto!(dest::gcmvector, bc::Broadcast.Broadcasted{Nothing})
    axes(dest) == axes(bc) || throwdm(axes(dest), axes(bc))
    # Performance optimization: broadcast!(identity, dest, A) is equivalent to copyto!(dest, A) if indices match
    if bc.f === identity && bc.args isa Tuple{AbstractArray} # only a single input argument to broadcast!
        A = bc.args[1]
        if axes(dest) == axes(A)
            return copyto!(dest, A)
        end
    end
    bc′ = Broadcast.preprocess(dest, bc)
    @simd for I in eachindex(bc′)
        #@inbounds dest[I] = bc′[I]
        @inbounds dest[I] = gcmarray_getindex_evalf(bc′,I)
    end
    return dest
end

##

"""
    findall(A::gcmarray{Bool})

Return a gcmvector of the true indices in A. This allows:
```
findall(A.<0) #gcmvector of CartesianIndex{2}
A[findall(A.<0)] #gcmvector of eltype(A)
view(A,findall(A.<0)) #CatView of eltype(A)

A[findall(A.<0)]=B[findall(A.<0)]
A[findall(A.<0)].=view(B,findall(A.<0))
A[findall(A.<0)].=NaN
```
"""
function findall(A::gcmarray{Bool,N}) where {N}
  tmpOut=Array{Array{CartesianIndex{2},1},N}(undef,size(A))
  for a in eachindex(A)
    tmpOut[a]=findall(A[a])
  end
  return gcmvector{CartesianIndex{2},N}(A.grid,tmpOut,A.fSize,A.fIndex)
end

function Base.getindex(A::gcmarray{T,N}, B::gcmvector{CartesianIndex{2},N}) where {T,N}
  tmpOut=Array{Array{T,1},N}(undef,size(A))
  for a in eachindex(A)
    tmpOut[a]=A[a][B[a]]
  end
  return gcmvector{T,N}(A.grid,tmpOut,A.fSize,A.fIndex)
end

function Base.view(A::gcmarray{T,N}, B::gcmvector{CartesianIndex{2},N}) where {T,N}
   tmpOut=missing
   for a in eachindex(A)
     ismissing(tmpOut) ? tmpOut=view(A[a],B[a]) : tmpOut=CatView(tmpOut,view(A[a],B[a]))
   end
   return tmpOut
end

function Base.setindex!(A::gcmarray{T,N}, B::gcmvector{T,N}, C::gcmvector{CartesianIndex{2},N}) where {T,N}
  D=A
  for a in eachindex(A)
    D[a][C[a]]=B[a]
  end
  return (A = D)
end

function Base.setindex!(A::gcmarray{T,N}, v::Number, C::gcmvector{CartesianIndex{2},N}) where {T,N}
  D=A
  for a in eachindex(A)
    D[a][C[a]].=v
  end
  return (A = D)
end
