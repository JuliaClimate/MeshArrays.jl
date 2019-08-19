
# ### gcmarray{T, N} and constructors
#
# gcmarray data structure. Available constructors:
#
# ```
# gcmarray{T,N}(grid::gcmgrid,f::Array{Array{T,N},1},
#          fSize::Array{NTuple{N, Int}},fIndex::Array{Int,1})
# gcmarray(grid::gcmgrid,
#          fSize::Array{NTuple{N, Int}},fIndex::Array{Int,1})
#
# gcmarray(grid::gcmgrid,::Type{T})
# gcmarray(grid::gcmgrid,::Type{T},n3::Int)
# ```

# +
struct gcmarray{T, N} <: AbstractMeshArray{T, N}
   grid::gcmgrid
   f::Array{Array{T,2},N}
   fSize::Array{NTuple{2, Int}}
   fIndex::Array{Int,1}
end

function gcmarray(grid::gcmgrid,::Type{T},
        fSize::Array{NTuple{2, Int}},
        fIndex::Array{Int,1}) where {T}
  nFaces=length(fIndex)
  f=Array{Array{T,2},1}(undef,nFaces)
  for a=1:nFaces
    f[a]=Array{T}(undef,fSize[a])
  end
  gcmarray{T,1}(grid,f,fSize,fIndex)
end

function gcmarray(grid::gcmgrid,::Type{T},
        fSize::Array{NTuple{2, Int}},
        fIndex::Array{Int,1},n3::Int) where {T}
  nFaces=length(fIndex)
  f=Array{Array{T,2},2}(undef,nFaces,n3)
  for a=1:nFaces; for i3=1:n3;
    f[a,i3]=Array{T}(undef,fSize[a]...)
  end; end;
  gcmarray{T,2}(grid,f,fSize,fIndex)
end

# +
function gcmarray(grid::gcmgrid,::Type{T}) where {T}
  nFaces=grid.nFaces
  fSize=grid.fSize
  fIndex=collect(1:grid.nFaces)
  gcmarray(grid,T,fSize,fIndex)
end

## this one needs to be redone: fSize 2D, f 3D
function gcmarray(grid::gcmgrid,::Type{T},n3::Int) where {T}
  nFaces=grid.nFaces
  fSize=grid.fSize
  fIndex=collect(1:grid.nFaces)
  gcmarray(grid,T,fSize,fIndex,n3)
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
    getindexetc(A::gcmarray{T, N}, I::Vararg{Union{Int,AbstractUnitRange,Colon}, N}) where {T,N}

Same as getindex but also returns the face size and index
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
  return view(A.f,I...)
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
