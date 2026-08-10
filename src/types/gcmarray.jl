
"""
    gcmarray{T, N, AT}

gcmarray data structure. Available constructors:

```
gcmarray{T,N,AT}(grid::gcmgrid,meta::varmeta,f::Array{AT,N},
         fSize::Array{NTuple{N, Int}},fIndex::Array{Int,1},v::String)

gcmarray(grid::gcmgrid,f::Array{Array{T,2},N}) where {T,N}
gcmarray(grid::gcmgrid,f::Array{Array{T,N},1}) where {T,N}

gcmarray(grid::gcmgrid,fSize::Array{NTuple{N, Int}},fIndex::Array{Int,1})
gcmarray(<same as above>,n3::Int)
gcmarray(<same as above>,n3::Int,n4::Int)

gcmarray(grid::gcmgrid)
gcmarray(grid::gcmgrid,::Type{T})
gcmarray(grid::gcmgrid,::Type{T},n3::Int)
gcmarray(grid::gcmgrid,::Type{T},n3::Int,n4::Int)

gcmarray(A::Array{T,N};meta::varmeta=defaultmeta)
```
"""
struct gcmarray{T, N, AT} <: AbstractMeshArray{T, N}
   grid::gcmgrid
   meta::varmeta
   f::OuterArray{AT,N}
   fSize::OuterArray{NTuple{2, Int}}
   fIndex::OuterArray{Int,1}
   version::String
end

function gcmarray(A::Array{T,N};
  meta::varmeta=defaultmeta) where {T, N}
  s=size(A)[1:2]
  fs=Array{NTuple{2, Int},1}(undef,1)
  fs[1]=s
  ios=[s[1] s[2]]
  γ=gcmgrid("","PeriodicDomain",1, fs, ios, T, read, write)
  if N==2
    B=gcmarray(γ,T;meta=meta)
  else
    B=gcmarray(γ,T,size(A)[3];meta=meta)
  end
  read!(A,B)
  B
end

function gcmarray(grid::gcmgrid,f::OuterArray{InnerArray{T,2},N};
                  meta::varmeta=defaultmeta) where {T, N}
  gcmarray{T,N,InnerArray{T,2}}(grid,meta,f,grid.fSize,collect(1:grid.nFaces),thisversion)
end

function gcmarray(grid::gcmgrid,f::OuterArray{InnerArray{T,N},1};
                  meta::varmeta=defaultmeta) where {T, N}
  nFaces=grid.nFaces
  if N>2
    n3=size(f[1],3); n4=size(f[1],4)
    g=OuterArray{InnerArray{T,2},3}(undef,nFaces,n3,n4)
    for I in eachindex(view(g,1:nFaces,1:n3,1:n4))
      g[I]=view(f[I[1]],:,:,I[2],I[3])
    end
    n4==1 ? g=dropdims(g,dims=3) : nothing
    gcmarray{T,ndims(g),InnerArray{T,2}}(grid,meta,g,grid.fSize,collect(1:nFaces),thisversion)
  else
    gcmarray{T,1,InnerArray{T,2}}(grid,meta,f,grid.fSize,collect(1:nFaces),thisversion)
  end
end

function gcmarray(grid::gcmgrid,::Type{T},
        fSize::Union{OuterArray{NTuple{2, Int},1},NTuple{2, Int}},
        fIndex::Union{OuterArray{Int,1},Int};
        meta::varmeta=defaultmeta) where {T}
  nFaces=length(fIndex)
  f=OuterArray{InnerArray{T,2},1}(undef,nFaces)
  isa(fSize,NTuple) ? fSize=[fSize] : nothing
  isa(fIndex,Int) ? fIndex=[fIndex] : nothing
  for a=1:nFaces
    f[a]=InnerArray{T}(undef,fSize[a])
  end
  gcmarray{T,1,InnerArray{T,2}}(grid,meta,f,fSize,fIndex,thisversion)
end

function gcmarray(grid::gcmgrid,::Type{T},
        fSize::Union{OuterArray{NTuple{2, Int},1},NTuple{2, Int}},
        fIndex::Union{OuterArray{Int,1},Int},n3::Int;
        meta::varmeta=defaultmeta) where {T}
  nFaces=length(fIndex)
  f=OuterArray{InnerArray{T,2},2}(undef,nFaces,n3)
  isa(fSize,NTuple) ? fSize=[fSize] : nothing
  isa(fIndex,Int) ? fIndex=[fIndex] : nothing
  for a=1:nFaces; for i3=1:n3
    f[a,i3]=InnerArray{T,2}(undef,fSize[a]...)
  end; end
  gcmarray{T,2,InnerArray{T,2}}(grid,meta,f,fSize,fIndex,thisversion)
end

function gcmarray(grid::gcmgrid,::Type{T},
        fSize::Union{OuterArray{NTuple{2, Int},1},NTuple{2, Int}},
        fIndex::Union{OuterArray{Int,1},Int},n3::Int,n4::Int;
        meta::varmeta=defaultmeta) where {T}
  nFaces=length(fIndex)
  f=OuterArray{InnerArray{T,2},3}(undef,nFaces,n3,n4)
  isa(fSize,NTuple) ? fSize=[fSize] : nothing
  isa(fIndex,Int) ? fIndex=[fIndex] : nothing
  for a=1:nFaces; for i4=1:n4; for i3=1:n3
    f[a,i3,i4]=InnerArray{T,2}(undef,fSize[a]...)
  end; end; end
  gcmarray{T,3,InnerArray{T,2}}(grid,meta,f,fSize,fIndex,thisversion)
end

# +
function gcmarray(grid::gcmgrid; meta::varmeta=defaultmeta)
  nFaces=grid.nFaces
  fSize=grid.fSize
  fIndex=collect(1:grid.nFaces)
  T=grid.ioPrec
  gcmarray(grid,T,fSize,fIndex,meta=meta)
end

function gcmarray(grid::gcmgrid,::Type{T};
                  meta::varmeta=defaultmeta) where {T}
  nFaces=grid.nFaces
  fSize=grid.fSize
  fIndex=collect(1:grid.nFaces)
  gcmarray(grid,T,fSize,fIndex,meta=meta)
end

function gcmarray(grid::gcmgrid,::Type{T},n3::Int;
                  meta::varmeta=defaultmeta) where {T}
  nFaces=grid.nFaces
  fSize=grid.fSize
  fIndex=collect(1:grid.nFaces)
  gcmarray(grid,T,fSize,fIndex,n3,meta=meta)
end

function gcmarray(grid::gcmgrid,::Type{T},n3::Int,n4::Int;
                  meta::varmeta=defaultmeta) where {T}
  nFaces=grid.nFaces
  fSize=grid.fSize
  fIndex=collect(1:grid.nFaces)
  gcmarray(grid,T,fSize,fIndex,n3,n4,meta=meta)
end

# -

# # Interface Methods

# +
#Base.dataids(A::gcmarray) = Base.dataids(A.f)
Base.dataids(A::gcmarray) = (Base.dataids(A.f)..., Base.dataids(A.fSize)..., Base.dataids(A.fIndex)...)

# +

function Base.getindex(A::AbstractMeshArray{T, N}, I::Vararg{Union{Int,Array{Int},AbstractUnitRange,Colon}, N}) where {T,N}
    J = 1:length(A.fIndex)
    !isa(I[1], Colon) ? J = J[I[1]] : nothing
    nFaces = length(J)
    
    tmpf = A.f[I...]
    
    if isa(tmpf, Array{eltype(A), 2})
        return tmpf  # Scalar index → single face
    else
        fSize_sub = A.fSize[J]
        fIndex_sub = A.fIndex[J]
        n3 = Int(length(tmpf) / nFaces)
        
        if n3 == 1
            # 1D result
            f_new = OuterArray{InnerArray{T,2},1}(undef, nFaces)
            for I_iter in eachindex(tmpf)
                f_new[I_iter] = view(tmpf[I_iter], :, :)
            end
            B = gcmarray{T, 1, InnerArray{T,2}}(
                A.grid, A.meta, f_new, fSize_sub, fIndex_sub, thisversion)
        else
            # 2D result
            f_new = OuterArray{InnerArray{T,2}, 2}(undef, nFaces, n3)
            for I_iter in eachindex(tmpf)
                f_new[I_iter] = view(tmpf[I_iter], :, :)
            end
            B = gcmarray{T, 2, InnerArray{T,2}}(
                A.grid, A.meta, f_new, fSize_sub, fIndex_sub, thisversion)
        end
        return B
    end
end

"""
    getindexetc(A::AbstractMeshArray{T, N}, I::Vararg{_}) where {T,N}

Same as getindex but also returns the face size and index
"""
function getindexetc(A::AbstractMeshArray{T, N}, I::Vararg{Union{Int,Array{Int},AbstractUnitRange,Colon}, N}) where {T,N}
    f=A[I...]
    fSize=A.fSize[I[1]]
    fIndex=A.fIndex[I[1]]
    return f,fSize,fIndex
end
# -

function Base.setindex!(A::AbstractMeshArray{T, N}, v, I::Vararg{Int, N}) where {T,N}
  return (A.f[I...] = v)
end

function Base.view(A::AbstractMeshArray{T, N}, I::Vararg{Union{Int,AbstractUnitRange,Colon}, N}) where {T,N}
    J = 1:length(A.fIndex)
    !isa(I[1], Colon) ? J = J[I[1]] : nothing
    nFaces = length(J)
    
    tmpf = view(A.f, I...)
    n3 = Int(length(tmpf) / nFaces)
    
    fSize_sub = A.fSize[J]
    fIndex_sub = A.fIndex[J]
    
    if n3 == 1
        # 1D result
        f_new = OuterArray{InnerArray{T,2},1}(undef, nFaces)
        for I_iter in eachindex(tmpf)
            f_new[I_iter] = view(tmpf[I_iter], :, :)
        end
        B = gcmarray{T, 1, InnerArray{T,2}}(
            A.grid, A.meta, f_new, fSize_sub, fIndex_sub, thisversion)
    else
        # 2D result
        f_new = OuterArray{InnerArray{T,2}, 2}(undef, nFaces, n3)
        for I_iter in eachindex(tmpf)
            f_new[I_iter] = view(tmpf[I_iter], :, :)
        end
        B = gcmarray{T, 2, InnerArray{T,2}}(
            A.grid, A.meta, f_new, fSize_sub, fIndex_sub, thisversion)
    end
    return B
end

# ### Custom pretty-printing, similar, and broadcast

function Base.show(io::IO, z::AbstractMeshArray{T, N}) where {T,N}
    if ~isa(z,MeshArrays.gcmfaces)&&~isa(z.meta.unit,Missing)
      printstyled(io, "  name        = ",color=:normal)
      printstyled(io, "$(z.meta.name)\n",color=:magenta)
      printstyled(io, "  unit        = ",color=:normal)
      printstyled(io, "$(z.meta.unit)\n",color=:magenta)
    end
    printstyled(io, "  data type   = ",color=:normal)
    printstyled(io, "$(eltype(z))\n",color=:magenta)
    if ~isa(z,MeshArrays.gcmfaces)
      printstyled(io, "  cell pos.   = ",color=:normal)
      printstyled(io, "$(z.meta.position)\n",color=:magenta)
    end
    printstyled(io, "  tile array  = ",color=:normal)
    printstyled(io, "$(size(z))\n",color=:cyan)
    printstyled(io, "  tile sizes  = ",color=:normal)
    #fs=(~isa(z,MeshArrays.gcmfaces) ? z.fSize : fsize(z))
    fs=z.fSize
    printstyled(io, "$(fs[1])\n",color=:cyan)
    for iFace=2:length(fs)
      printstyled(io, "                ",color=:normal)
      printstyled(io, "$(fs[iFace])\n",color=:cyan)
    end
    printstyled(io, "  grid class  = ",color=:normal)
    printstyled(io, "$(z.grid.class)\n",color=:green)
    printstyled(io, "  MeshArray   = ",color=:normal)
    printstyled(io, "$(typeof(z)) \n",color=:green)
    if ~isa(z,MeshArrays.gcmfaces)
     printstyled(io, "  version     = ",color=:normal)
     printstyled(io, "$(z.version) \n",color=:green)
    end
  return
end

import Base: display; display(X::AbstractMeshArray)=show(X)

"""
    similar(A::gcmarray; m::varmeta=defaultmeta, allocate=false)

Create a gcmarray with the same structure and type as `A`.

If `allocate=true`, eagerly allocates all face arrays sized according to `A.fSize`.
If `allocate=false` (default), creates lazy empty placeholders; use when the array
will be filled in-place (e.g., `read!`, `readtiles`).
"""
function Base.similar(A::gcmarray; m::varmeta=defaultmeta,allocate=false)
    ElType = eltype(A)
    nFaces = length(A.fIndex)

    if ndims(A) == 1
      if allocate
        # 1D case: (nFaces,)
        f = OuterArray{InnerArray{ElType,2},1}(undef, nFaces)
        for a in 1:nFaces
            f[a] = InnerArray{ElType}(undef, A.fSize[a]...)
        end
      else
        f = fill(InnerArray{ElType}(undef, 0,0),nFaces)
      end
      B = gcmarray{ElType, 1, InnerArray{ElType,2}}(
          A.grid, m, f, A.fSize, A.fIndex, thisversion)
    elseif ndims(A) == 2
        # 2D case: (nFaces, n3)
        n3 = size(A, 2)
        if allocate
        f = OuterArray{InnerArray{ElType,2}, 2}(undef, nFaces, n3)
        for a in 1:nFaces
            for i3 in 1:n3
                f[a, i3] = InnerArray{ElType}(undef, A.fSize[a]...)
            end
        end
        else
        f = fill(InnerArray{ElType}(undef, 0,0),nFaces,n3)
        end
        B = gcmarray{ElType, 2, InnerArray{ElType,2}}(
            A.grid, m, f, A.fSize, A.fIndex, thisversion)
    else  # ndims(A) == 3
        # 3D case: (nFaces, n3, n4)
        n3 = size(A, 2)
        n4 = size(A, 3)
        if allocate
        f = OuterArray{InnerArray{ElType,2}, 3}(undef, nFaces, n3, n4)
        for a in 1:nFaces
            for i4 in 1:n4
                for i3 in 1:n3
                    f[a, i3, i4] = InnerArray{ElType}(undef, A.fSize[a]...)
                end
            end
        end
        else
        f = fill(InnerArray{ElType}(undef, 0,0),nFaces,n3,n4)
        end
        B = gcmarray{ElType, 3, InnerArray{ElType,2}}(
            A.grid, m, f, A.fSize, A.fIndex, thisversion)

    end
    return B
end

# ### Customize broadcasting

Base.BroadcastStyle(::Type{<:AbstractMeshArray}) = Broadcast.ArrayStyle{AbstractMeshArray}()

function Base.similar(bc::Broadcast.Broadcasted{Broadcast.ArrayStyle{AbstractMeshArray}},
                      ::Type{ElType}) where ElType
    A = find_gcmarray(bc)
    nFaces = length(A.fIndex)
    if ndims(A) == 1
        f = fill(InnerArray{ElType}(undef, 0, 0), nFaces)
        return gcmarray{ElType, 1, InnerArray{ElType,2}}(
            A.grid, defaultmeta, f, A.fSize, A.fIndex, thisversion)
    elseif ndims(A) == 2
        n3 = size(A, 2)
        f = fill(InnerArray{ElType}(undef, 0, 0), nFaces, n3)
        return gcmarray{ElType, 2, InnerArray{ElType,2}}(
            A.grid, defaultmeta, f, A.fSize, A.fIndex, thisversion)
    else
        n3, n4 = size(A, 2), size(A, 3)
        f = fill(InnerArray{ElType}(undef, 0, 0), nFaces, n3, n4)
        return gcmarray{ElType, 3, InnerArray{ElType,2}}(
            A.grid, defaultmeta, f, A.fSize, A.fIndex, thisversion)
    end
end

find_gcmarray(bc::Base.Broadcast.Broadcasted) = find_gcmarray(bc.args)
find_gcmarray(args::Tuple) = find_gcmarray(find_gcmarray(args[1]), Base.tail(args))
find_gcmarray(x) = x
find_gcmarray(a::AbstractMeshArray, rest) = a
find_gcmarray(::Any, rest) = find_gcmarray(rest)

# Recursively rewrite a Broadcasted tree for face I:
# replace each AbstractMeshArray leaf with its per-face InnerArray;
# scalars and plain Arrays pass through unchanged.
function _bc_face(bc::Broadcast.Broadcasted, I)
    new_args = map(a -> _bc_face(a, I), bc.args)
    Broadcast.Broadcasted(bc.f, new_args)
end
_bc_face(a::AbstractMeshArray, I) = ndims(a.f) == 1 ? a.f[I[1]] : a.f[I]
_bc_face(a, I) = a

import Base: copyto!

@inline function Base.copyto!(dest::AbstractMeshArray,
                               bc::Broadcast.Broadcasted{Broadcast.ArrayStyle{AbstractMeshArray}})
    for I in CartesianIndices(dest.f)
        face_bc = _bc_face(bc, I)
        if size(dest.f[I]) == (0, 0)
            dest.f[I] = Base.Broadcast.materialize(face_bc)
        else
            Base.Broadcast.materialize!(dest.f[I], face_bc)
        end
    end
    return dest
end
