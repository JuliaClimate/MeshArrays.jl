
"""
    gcmarray{T, N}

gcmarray data structure. Available constructors:

```
gcmarray{T,N}(grid::gcmgrid,f::Array{Array{T,N},1},
         fSize::Array{NTuple{N, Int}},fIndex::Array{Int,1})

gcmarray(grid::gcmgrid,f::Array{Array{T,2},N}) where {T,N}
gcmarray(grid::gcmgrid,f::Array{Array{T,N},1}) where {T,N}

gcmarray(grid::gcmgrid,fSize::Array{NTuple{N, Int}},fIndex::Array{Int,1})
gcmarray(<same as above>,n3::Int)
gcmarray(<same as above>,n3::Int,n4::Int)

gcmarray(grid::gcmgrid)
gcmarray(grid::gcmgrid,::Type{T})
gcmarray(grid::gcmgrid,::Type{T},n3::Int)
gcmarray(grid::gcmgrid,::Type{T},n3::Int,n4::Int)
```
"""
struct gcmarray{T, N} <: AbstractMeshArray{T, N}
   grid::gcmgrid
   f::Array{Array{T,2},N}
   fSize::Array{NTuple{2, Int}}
   fIndex::Array{Int,1}
end

function gcmarray(grid::gcmgrid,f::Array{Array{T,2},N}) where {T, N}
  gcmarray{T,N}(grid,f,grid.fSize,collect(1:grid.nFaces))
end

function gcmarray(grid::gcmgrid,f::Array{Array{T,N},1}) where {T, N}
  nFaces=grid.nFaces
  if N>2
    n3=size(f[1],3); n4=size(f[1],4);
    g=Array{Array{T,2},3}(undef,nFaces,n3,n4)
    for I in eachindex(view(g,1:nFaces,1:n3,1:n4))
      g[I]=view(f[I[1]],:,:,I[2],I[3])
    end
    n4==1 ? g=dropdims(g,dims=3) : nothing
    gcmarray{T,ndims(g)}(grid,g,grid.fSize,collect(1:nFaces))
  else
    gcmarray{T,1}(grid,f,grid.fSize,collect(1:nFaces))
  end
end

function gcmarray(grid::gcmgrid,::Type{T},
        fSize::Union{Array{NTuple{2, Int},1},NTuple{2, Int}},
        fIndex::Union{Array{Int,1},Int}) where {T}
  nFaces=length(fIndex)
  f=Array{Array{T,2},1}(undef,nFaces)
  isa(fSize,NTuple) ? fSize=[fSize] : nothing
  isa(fIndex,Int) ? fIndex=[fIndex] : nothing
  for a=1:nFaces
    f[a]=Array{T}(undef,fSize[a])
  end
  gcmarray{T,1}(grid,f,fSize,fIndex)
end

function gcmarray(grid::gcmgrid,::Type{T},
        fSize::Union{Array{NTuple{2, Int},1},NTuple{2, Int}},
        fIndex::Union{Array{Int,1},Int},n3::Int) where {T}
  nFaces=length(fIndex)
  f=Array{Array{T,2},2}(undef,nFaces,n3)
  isa(fSize,NTuple) ? fSize=[fSize] : nothing
  isa(fIndex,Int) ? fIndex=[fIndex] : nothing
  for a=1:nFaces; for i3=1:n3;
    f[a,i3]=Array{T}(undef,fSize[a]...)
  end; end;
  gcmarray{T,2}(grid,f,fSize,fIndex)
end

function gcmarray(grid::gcmgrid,::Type{T},
        fSize::Union{Array{NTuple{2, Int},1},NTuple{2, Int}},
        fIndex::Union{Array{Int,1},Int},n3::Int,n4::Int) where {T}
  nFaces=length(fIndex)
  f=Array{Array{T,2},3}(undef,nFaces,n3,n4)
  isa(fSize,NTuple) ? fSize=[fSize] : nothing
  isa(fIndex,Int) ? fIndex=[fIndex] : nothing
  for a=1:nFaces; for i4=1:n4; for i3=1:n3;
    f[a,i3,i4]=Array{T}(undef,fSize[a]...)
  end; end; end;
  gcmarray{T,3}(grid,f,fSize,fIndex)
end

# +
function gcmarray(grid::gcmgrid)
  nFaces=grid.nFaces
  fSize=grid.fSize
  fIndex=collect(1:grid.nFaces)
  T=grid.ioPrec
  gcmarray(grid,T,fSize,fIndex)
end

function gcmarray(grid::gcmgrid,::Type{T}) where {T}
  nFaces=grid.nFaces
  fSize=grid.fSize
  fIndex=collect(1:grid.nFaces)
  gcmarray(grid,T,fSize,fIndex)
end

function gcmarray(grid::gcmgrid,::Type{T},n3::Int) where {T}
  nFaces=grid.nFaces
  fSize=grid.fSize
  fIndex=collect(1:grid.nFaces)
  gcmarray(grid,T,fSize,fIndex,n3)
end

function gcmarray(grid::gcmgrid,::Type{T},n3::Int,n4::Int) where {T}
  nFaces=grid.nFaces
  fSize=grid.fSize
  fIndex=collect(1:grid.nFaces)
  gcmarray(grid,T,fSize,fIndex,n3,n4)
end

# -

# # Interface Methods

# +
Base.size(A::gcmarray) = size(A.f)
Base.size(A::gcmarray, dim::Integer) = size(A)[dim]

# +
function Base.getindex(A::gcmarray{T, N}, I::Vararg{Union{Int,AbstractUnitRange,Colon}, N}) where {T,N}
    return A.f[I...]
end

"""
    getindexetc(A::gcmarray, I::Vararg{_}) where {T,N}

Same as getindex but also returns the face size and index

(needed in other types?)
"""
function getindexetc(A::gcmarray{T, N}, I::Vararg{Union{Int,AbstractUnitRange,Colon}, N}) where {T,N}
    f=A.f[I...]
    fSize=A.fSize[I[1]]
    fIndex=A.fIndex[I[1]]
    return f,fSize,fIndex
end
# -

function Base.setindex!(A::gcmarray{T, N}, v, I::Vararg{Int, N}) where {T,N}
  return (A.f[I...] = v)
end

function Base.view(A::gcmarray{T, N}, I::Vararg{Union{Int,AbstractUnitRange,Colon}, N}) where {T,N}
  J=1:length(A.fIndex)
  !isa(I[1],Colon) ? J=J[I[1]] : nothing
  nFaces=length(J)

  tmpf=view(A.f,I...)
  n3=Int(length(tmpf)/nFaces) #length(tmpf)>nFaces ? n3=Int(length(tmpf)/nFaces) : n3=1

  K=(A.grid,eltype(A),A.fSize[J],A.fIndex[J])
  n3>1 ? tmp=gcmarray(K...,n3) : tmp=gcmarray(K...)
  for I in eachindex(tmpf); tmp.f[I] = view(tmpf[I],:,:); end

  return tmp
end

# ### Custom pretty-printing, similar, and broadcast

function Base.show(io::IO, z::gcmarray{T, N}) where {T,N}
    printstyled(io, " gcmarray \n",color=:normal)
    fs=z.fSize
    printstyled(io, "  grid type   = ",color=:normal)
    printstyled(io, "$(z.grid.class)\n",color=:blue)
    printstyled(io, "  array size  = ",color=:normal)
    printstyled(io, "$(size(z))\n",color=:blue)

    if ~isassigned(z.f);
      printstyled(io, "  data type   = ",color=:normal)
      printstyled(io, "unassigned\n",color=:green)
      printstyled(io, "  tile sizes  = ",color=:normal)
      printstyled(io, "unassigned\n",color=:green)
    else
      printstyled(io, "  data type   = ",color=:normal)
      printstyled(io, "$(typeof(z.f[1][1]))\n",color=:blue)
      printstyled(io, "  face sizes  = ",color=:normal)
      printstyled(io, "$(fs[1])\n",color=:blue)
      for iFace=2:length(fs)
        printstyled(io, "                ",color=:normal)
        printstyled(io, "$(fs[iFace])\n",color=:blue)
      end
    end

  return
end

function Base.similar(A::gcmarray)
    if ndims(A)==1
        B=gcmarray(A.grid,eltype(A),A.fSize,A.fIndex)
    else
        B=gcmarray(A.grid,eltype(A),A.fSize,A.fIndex,size(A,2))
    end
    return B
end

# ### Customize broadcasting

Base.BroadcastStyle(::Type{<:gcmarray}) = Broadcast.ArrayStyle{gcmarray}()

function Base.similar(bc::Broadcast.Broadcasted{Broadcast.ArrayStyle{gcmarray}}, ::Type{ElType}) where ElType
  # Scan the inputs for the gcmarray:
  A = find_gcmarray(bc)
  # Create the gcmarray output:
  #similar(A)
  if ndims(A)==1
        B=gcmarray(A.grid,eltype(A),A.fSize,A.fIndex)
  else
        B=gcmarray(A.grid,eltype(A),A.fSize,A.fIndex,size(A,2))
  end
  return B
end

find_gcmarray(bc::Base.Broadcast.Broadcasted) = find_gcmarray(bc.args)
find_gcmarray(args::Tuple) = find_gcmarray(find_gcmarray(args[1]), Base.tail(args))
find_gcmarray(x) = x
find_gcmarray(a::gcmarray, rest) = a
find_gcmarray(::Any, rest) = find_gcmarray(rest)

####

import Base: copyto!

# Specialize this method if all you want to do is specialize on typeof(dest)
@inline function copyto!(dest::MeshArrays.gcmarray, bc::Broadcast.Broadcasted{Nothing})
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
        @inbounds dest[I] = my_getindex_evalf(bc′,I)
    end
    return dest
end

function my_getindex_evalf(bc,I)
  @boundscheck checkbounds(bc, I)
  args = Broadcast._getindex(bc.args, I)
  return bc.f.(args...)
end

###

import Base: +, -, *, /
import Base: isnan, isinf, isfinite
import Base: maximum, minimum, sum, fill

#+(a::gcmarray,b::gcmarray) = a.+b
function +(a::gcmarray,b::gcmarray)
  c=similar(a)
  for I in eachindex(a)
    c[I] = a[I] + b[I]
  end
  return c
end

function -(a::gcmarray,b::gcmarray)
  c=similar(a)
  for I in eachindex(a)
    c[I] = a[I] - b[I]
  end
  return c
end

function /(a::gcmarray,b::gcmarray)
  c=similar(a)
  for I in eachindex(a)
    c[I] = a[I] ./ b[I]
  end
  return c
end

function *(a::gcmarray,b::gcmarray)
  c=similar(a)
  for I in eachindex(a)
    c[I] = a[I] .* b[I]
  end
  return c
end

isnan(a::gcmarray) = isnan.(a)
isinf(a::gcmarray) = isinf.(a)
isfinite(a::gcmarray) = isfinite.(a)

function maximum(a::gcmarray)
  c=-Inf;
  for I in eachindex(a)
    c = max(c,maximum(a[I]))
  end
  return c
end

function minimum(a::gcmarray)
  c=-Inf;
  for I in eachindex(a)
    c = min(c,minimum(a[I]))
  end
  return c
end

function sum(a::gcmarray)
  c=0.0
  for I in eachindex(a)
    c = c + sum(a[I])
  end
  return c
end

function fill(val::Any,a::gcmarray)
  c=similar(a)
  for I in eachindex(a)
    c[I] = fill(val,size(a[I]))
  end
  return c
end

###

"""
    nFacesEtc(a::gcmarray)

Return nFaces, n3 (1 in 2D case; >1 otherwise)
"""
function nFacesEtc(a::gcmarray)
  nFaces=length(a.fIndex)
  ndims(a.f)>1 ? n3=size(a.f,2) : n3=1
  ndims(a.f)>2 ? n4=size(a.f,3) : n4=1
  return nFaces, n3, n4
end

"""
    fijind(A::gcmarray,ij::Int)

Compute face and local indices (f,j,k?) from global index (ij?).
"""
#missing functionality
