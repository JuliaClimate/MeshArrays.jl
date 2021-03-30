## Module associated with the MeshArrays package

module MeshArrays

using Pkg, Pkg.Artifacts
thistoml=joinpath(dirname(pathof(MeshArrays)), "..", "Project.toml")
thisversion=Pkg.TOML.parsefile(thistoml)["version"]

p=dirname(pathof(MeshArrays))
artifact_toml = joinpath(p, "../Artifacts.toml")

GRID_LLC90_hash = artifact_hash("GRID_LLC90", artifact_toml)
GRID_LLC90 = joinpath(artifact_path(GRID_LLC90_hash)*"/","GRID_LLC90-1.1/")
GRID_LL360_hash = artifact_hash("GRID_LL360", artifact_toml)
GRID_LL360 = joinpath(artifact_path(GRID_LL360_hash)*"/","GRID_LL360-1.0/")
GRID_CS32_hash = artifact_hash("GRID_CS32", artifact_toml)
GRID_CS32 = joinpath(artifact_path(GRID_CS32_hash)*"/","GRID_CS32-1.1/")

#include("Types.jl");
using Unitful

"""
    AbstractMeshArray{T, N}

Subtype of AbstractArray{T, N}
"""
abstract type AbstractMeshArray{T, N} <: AbstractArray{T, N} end

"""
    gcmgrid

gcmgrid data structure. Available constructors:

```
gcmgrid(path::String, class::String,
        nFaces::Int, fSize::Array{NTuple{2, Int},1},
        ioSize::Array{Int64,2}, ioPrec::Type,
        read::Function, write::Function)
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

"""
    varmeta

varmeta data structure. By default, `unit` is `1.0` (non-dimensional), `position`
is `fill(0.5,3)` (cell center), and `name` / `long_name` is unknown.

Available constructors:

```
varmeta(unit::Union{Unitful.AbstractQuantity,Number},position::Array{Float64,1},
        name::String,long_name::String)
```

And:

```defaultmeta = varmeta(1.0,fill(0.5,3),"unknown","unknown")```

"""
struct varmeta
  unit::Union{Unitful.Units,Number,Missing}
  position::Array{Float64,1}
  name::String
  long_name::String
end

defaultmeta = varmeta(missing,fill(0.5,3),"unknown","unknown")

## concrete types and MeshArray alias:

#include("Type_gcmfaces.jl");
## gcmfaces type definition + methods

## type definition

"""
    gcmfaces{T, N}

gcmfaces data structure. Available constructors:

```
gcmfaces{T,N}(grid::gcmgrid,f::Array{Array{T,N},1},
         fSize::Array{NTuple{N, Int}}, aSize::NTuple{N,Int})

gcmfaces(grid::gcmgrid,v1::Array{Array{T,N},1}) where {T,N}
gcmfaces(grid::gcmgrid,::Type{T},
         fSize::Array{NTuple{N, Int}}, aSize::NTuple{N,Int}) where {T,N}

gcmfaces(grid::gcmgrid)
gcmfaces(grid::gcmgrid,::Type{T})
gcmfaces(grid::gcmgrid,::Type{T},n3::Int)
```
"""
struct gcmfaces{T, N} <: AbstractMeshArray{T, N}
   grid::gcmgrid
   f::Array{Array{T,N},1}
   fSize::Array{NTuple{N, Int}}
   aSize::NTuple{N, Int}
end

"""
    gcmsubset{T, N}

gcmsubset data structure for subsets of gcmfaces. Available constructors:

```
gcmsubset{T,N}(grid::gcmgrid,f::Array{Array{T,N},1},
               fSize::Array{NTuple{N, Int}},aSize::NTuple{N, Int},
               i::Array{Array{T,N},1},iSize::Array{NTuple{N, Int}})
gcmsubset(grid::gcmgrid,::Type{T},fSize::Array{NTuple{N, Int}},
          aSize::NTuple{N,Int},dims::NTuple{N,Int}) where {T,N}
```
"""
struct gcmsubset{T, N} <: AbstractMeshArray{T, N}
   grid::gcmgrid
   f::Array{Array{T,N},1}
   fSize::Array{NTuple{N, Int}}
   aSize::NTuple{N, Int}
   i::Array{Array{T,N},1}
   iSize::Array{NTuple{N, Int}}
end

## additional constructors for gcmfaces

function gcmfaces(grid::gcmgrid,::Type{T},
  fSize::Array{NTuple{N, Int}},
  aSize::NTuple{N,Int}) where {T,N}
  nFaces=grid.nFaces
  f=Array{Array{T,N},1}(undef,nFaces)
  for a=1:nFaces
    f[a]=Array{T}(undef,fSize[a])
  end
  gcmfaces{T,N}(grid,f,fSize,aSize)
end

function gcmfaces(grid::gcmgrid,::Type{T}) where {T,N}
  nFaces=grid.nFaces
  fSize=grid.fSize
  aSize=(prod(grid.ioSize),1)
  gcmfaces(grid,T,fSize,aSize)
end

function gcmfaces(grid::gcmgrid,::Type{T},n3::Int) where {T,N}
  nFaces=grid.nFaces
  fSize=Array{NTuple{3, Int},1}(undef,nFaces)
  for a=1:nFaces
    fSize[a]=(grid.fSize[a][1],grid.fSize[a][2],n3)
  end
  aSize=(prod(grid.ioSize),1,n3)
  gcmfaces(grid,T,fSize,aSize)
end

#other possibilities:
#gcmfaces{T,N}(grid::gcmgrid)
#gcmfaces(grid::gcmgrid,::Type{T}) where {T}
#gcmfaces(grid::gcmgrid,::Type{T},n3::Int) where {T}

function gcmfaces(grid::gcmgrid,
  v1::Array{Array{T,N},1}) where {T,N}
  fSize=fsize(v1)
  aSize=fsize(v1,0)
  gcmfaces{T,N}(grid,deepcopy(v1),fSize,aSize)
end

#should this be called similar? deepcopy?
#function gcmfaces(A::AbstractMeshArray{T, N}) where {T,N}
#  fSize=fsize(A)
#  aSize=size(A)
#  grid=A.grid
#  gcmfaces{T,N}(grid,deepcopy(A.f),fSize,aSize)
#end

function gcmfaces(grid::gcmgrid)
  T=grid.ioPrec
  fSize=grid.fSize
  aSize=(prod(grid.ioSize),1)
  gcmfaces(grid,T,fSize,aSize)
end

#function gcmfaces()
#  T=Float64
#  fSize=[(90, 270), (90, 270), (90, 90), (270, 90), (270, 90)]
#  aSize=(105300, 1)
#  grid=gcmgrid("GRID_LLC90/", "LatLonCap", 5, fSize, [90 1170], T, read, write)
#
#  gcmfaces(grid,T,fSize,aSize)
#end

## additional constructors for gcmsubset

#maybe: replace this constructor with one that gets A and sets f to view(A.f)
function gcmsubset(grid::gcmgrid,::Type{T},
  fSize::Array{NTuple{N, Int}},aSize::NTuple{N,Int},
  dims::NTuple{N,Int}) where {T,N}
  nFaces=grid.nFaces
  f=Array{Array{T,N},1}(undef,nFaces)
  i=Array{Array{T,N},1}(undef,nFaces)
  iSize=Array{NTuple{N, Int},1}(undef,nFaces)
  for a=1:nFaces
    f[a]=Array{T}(undef,fSize[a])
    #eventually I will distribute across faces; for now I just use face 1:
    a==1 ? nloc=dims[1] : nloc=0
    tmp1=Base.tail(dims)
    iSize[a]=(nloc,tmp1...)
    i[a]=Array{T}(undef,iSize[a])
  end
  gcmsubset{T,N}(grid,f,fSize,aSize,i,iSize)
end

## Convenience functions

"""
    fijind(A::gcmfaces,ij::Int)

Compute face and local indices (f,j,k) from global index (ij).

(needed in other types?)
"""
function fijind(A::gcmfaces,ij::Int)
  f=0
  j=0
  k=0
  tmp1=0
  for iFace=1:A.grid.nFaces
    tmpsize=fsize(A,iFace)
    tmp11=tmpsize[1]*tmpsize[2]
    tmp2=tmp1+tmp11
    if tmp1<ij<=tmp2
      f=iFace;
      tmp3=(ij-tmp1);
      k=Int(ceil(tmp3/tmpsize[1]))
      j=Int(tmp3-tmpsize[1]*(k-1))
    end;
    tmp1=tmp1+tmp11
  end
  return (f,j,k)
end

"""
    fsize(A::Union{gcmfaces{T, N},gcmsubset{T, N}}) where {T,N}

Return vector of face array sizes. Other methods:
```
fsize(A::Union{gcmfaces{T, N},gcmsubset{T, N}},i::Int) where {T,N}
fsize(A::Array{Array{T,N},1}) where {T,N}
fsize(A::Array{Array{T,N},1},i::Int) where {T,N}
```
"""
#deprecate documentation
function fsize(A::Union{gcmfaces{T, N},gcmsubset{T, N}}) where {T,N}
  fs=Array{NTuple{N, Int}}(undef,A.grid.nFaces)
  for i=1:A.grid.nFaces
    fs[i]=size(A.f[i]);
  end
  return fs
end

function fsize(A::Union{gcmfaces{T, N},gcmsubset{T, N}},i::Int) where {T,N}
  if i>0
    fs=size(A.f[i])
  else
    tmp1=0
    for i=1:A.grid.nFaces
      tmp1=tmp1+size(A.f[i],1)*size(A.f[i],2)
    end
    tmp2=size(A.f[1])
    fs=(tmp1,1,tmp2[3:end]...)
  end
end

function fsize(A::Array{Array{T,N},1}) where {T,N}
  fs=Array{NTuple{N, Int}}(undef,length(A))
  for i=1:length(A)
    fs[i]=size(A[i]);
  end
  return fs
end

function fsize(A::Array{Array{T,N},1},i::Int) where {T,N}
  if i>0
    fs=size(A[i])
  else
    tmp1=0
    for i=1:length(A)
      tmp1=tmp1+size(A[i],1)*size(A[i],2)
    end
    tmp2=size(A[1])
    fs=(tmp1,1,tmp2[3:end]...)
  end
end

## Interface Methods

Base.size(A::gcmfaces) = fsize(A, 0)
Base.size(A::gcmfaces, dim::Integer) = fsize(A, 0)[dim]
Base.size(A::gcmsubset) = fsize(A.i, 0)
Base.size(A::gcmsubset, dim::Integer) = fsize(A.i, 0)[dim]

#

function Base.getindex(A::Union{gcmfaces{T, N},gcmsubset{T, N}}, I::Vararg{Union{Int,AbstractUnitRange,Colon}, N}) where {T,N}
  if typeof(I[1])<:Int
    (f,i,j)=fijind(A,I[1])
    J=Base.tail(Base.tail(I))
    J=(i,j,J...)
    val=A.f[f][J...]
  elseif typeof(I[1])<:AbstractUnitRange
    val=similar(A,eltype(A),length.(I))
    for iFace=1:A.grid.nFaces
      @views val.f[iFace]=A.f[iFace]
    end
    #eventually I will distribute across faces; for now I just use face 1:
    k=0
    J=Base.tail(Base.tail(I))
    for kk=I[1]
      k+=1
      (f,i,j)=fijind(A,kk)
      tmp1=(k,1,J...)
      tmp2=(i,j,J...)
      val.i[1][tmp1...]=A.f[f][tmp2...]
    end
  elseif typeof(I[1])<:Colon
    #should this rather be a copy as the above?
    val=view(A,I...)
  else
    er1=typeof(A)
    er2=typeof(I[1])
    error("getindex not yet implemented for "*"$er1"*" and "*"$er2"*" indices")
  end
  return val
end

function Base.getindex(A::gcmsubset{T, N}, I::Vararg{Union{Int,AbstractUnitRange,Colon}, N}) where {T,N}
  #eventually I will distribute across faces; for now I just use face 1:
  return getindex(A.i[1],I...)
end

#

function Base.setindex!(A::Union{gcmfaces{T, N},gcmsubset{T, N}}, v, I::Vararg{Int, N}) where {T,N}
  (f,i,j)=fijind(A,I[1])
  J=Base.tail(Base.tail(I))
  J=(i,j,J...)
  return (A.f[f][J...] = v)
end

function Base.setindex!(A::gcmsubset{T, N}, v, I::Vararg{Int, N}) where {T,N}
  #eventually I will distribute across faces; for now I just use face 1:
  return (A.i[1][I...] = v)
end

## view

function Base.view(a::Union{gcmfaces{T, N},gcmsubset{T, N}}, I::Vararg{Union{Int,AbstractUnitRange,Colon}, N}) where {T,N}
  nFaces=a.grid.nFaces
  grTopo=a.grid.class
  if !isa(I[1],Colon)|!isa(I[2],Colon)
    J=Base.tail(Base.tail(I))
    J=(:,:,J...)
  else
    J=I
  end
  Nout=length(size(view(a.f[1],J...)));
  v1=Array{Array{T,Nout}}(undef,nFaces);
  for iFace=1:nFaces
    v1[iFace]=view(a.f[iFace],J...);
  end
  c=gcmfaces(a.grid,v1);
  return c;
end

# Custom pretty-printing

function Base.show(io::IO, z::Union{gcmfaces{T, N},gcmsubset{T, N}}) where {T,N}

#    @printf io " MeshArrays instance with \n"
    if isa(z,gcmfaces)
      printstyled(io, " gcmfaces array \n",color=:normal)
      nm="face"
      fs=fsize(z.f)
    elseif isa(z,gcmsubset)
      printstyled(io, " gcmsubset array \n",color=:normal)
      fs=fsize(z.i)
      nm="subset"
    else
      error("unknown type")
    end
    printstyled(io, "  grid type   = ",color=:normal)
    printstyled(io, "$(z.grid.class)\n",color=:blue)
    printstyled(io, "  # of faces  = ",color=:normal)
    printstyled(io, "$(z.grid.nFaces)\n",color=:blue)
    if ~isassigned(z.f);
      printstyled(io, "  data type   = ",color=:normal)
      printstyled(io, "unassigned\n",color=:green)
      printstyled(io, "  face sizes  = ",color=:normal)
      printstyled(io, "unassigned\n",color=:green)
    else
      printstyled(io, "  data type   = ",color=:normal)
      printstyled(io, "$(typeof(z.f[1][1]))\n",color=:blue)
      printstyled(io, "  $(nm) sizes  = ",color=:normal)
      printstyled(io, "$(fs[1])\n",color=:blue)
      for iFace=2:z.grid.nFaces
        printstyled(io, "                ",color=:normal)
        printstyled(io, "$(fs[iFace])\n",color=:blue)
      end
    end

    return
end

#

function Base.similar(A::gcmfaces, ::Type{T}, dims::Dims) where {T}
  if prod(dims)==length(A)
    B=gcmfaces(A.grid,T,A.fSize,A.aSize)
  else
    B=gcmsubset(A.grid,T,A.fSize,A.aSize,dims)
  end
end

Base.BroadcastStyle(::Type{<:gcmfaces}) = Broadcast.ArrayStyle{gcmfaces}()

function Base.similar(bc::Broadcast.Broadcasted{Broadcast.ArrayStyle{gcmfaces}}, ::Type{ElType}) where ElType
    # Scan the inputs for the gcmfaces:
    A = find_gcmfaces(bc)
    # Create the gcmfaces output:
    similar(A,ElType,A.aSize)
end

find_gcmfaces(bc::Base.Broadcast.Broadcasted) = find_gcmfaces(bc.args)
find_gcmfaces(args::Tuple) = find_gcmfaces(find_gcmfaces(args[1]), Base.tail(args))
find_gcmfaces(x) = x
find_gcmfaces(a::gcmfaces, rest) = a
find_gcmfaces(::Any, rest) = find_gcmfaces(rest)

#

function Base.similar(A::gcmsubset, ::Type{T}, dims::Dims) where {T}
    B=gcmsubset(A.grid,T,A.fSize,A.aSize,dims[1])
end

Base.BroadcastStyle(::Type{<:gcmsubset}) = Broadcast.ArrayStyle{gcmsubset}()

function Base.similar(bc::Broadcast.Broadcasted{Broadcast.ArrayStyle{gcmsubset}}, ::Type{ElType}) where ElType
    # Scan the inputs for the gcmsubset:
    A = find_gcmsubset(bc)
    # Create the gcmsubset output:
    similar(A)
end

find_gcmsubset(bc::Base.Broadcast.Broadcasted) = find_gcmsubset(bc.args)
find_gcmsubset(args::Tuple) = find_gcmsubset(find_gcmsubset(args[1]), Base.tail(args))
find_gcmsubset(x) = x
find_gcmsubset(a::gcmsubset, rest) = a
find_gcmsubset(::Any, rest) = find_gcmsubset(rest)

###

function nFacesEtc(a::gcmfaces)
  nFaces=length(a.f)
  ndims(a.f[1])>2 ? n3=size(a.f[1],3) : n3=1
  return nFaces, n3
end
OuterArray{T,N}=Array{T,N} where {T,N}
InnerArray{T,N}=Array{T,N} where {T,N}
#include("Type_gcmarray.jl");

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

function gcmarray(grid::gcmgrid,f::OuterArray{InnerArray{T,2},N};
                  meta::varmeta=defaultmeta) where {T, N}
  gcmarray{T,N,InnerArray{T,2}}(grid,meta,f,grid.fSize,collect(1:grid.nFaces),thisversion)
end

function gcmarray(grid::gcmgrid,f::OuterArray{InnerArray{T,N},1};
                  meta::varmeta=defaultmeta) where {T, N}
  nFaces=grid.nFaces
  if N>2
    n3=size(f[1],3); n4=size(f[1],4);
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
  for a=1:nFaces; for i3=1:n3;
    f[a,i3]=InnerArray{T,2}(undef,fSize[a]...)
  end; end;
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
  for a=1:nFaces; for i4=1:n4; for i3=1:n3;
    f[a,i3,i4]=InnerArray{T,2}(undef,fSize[a]...)
  end; end; end;
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
Base.size(A::gcmarray) = size(A.f)
Base.size(A::gcmarray, dim::Integer) = size(A)[dim]

# +
function Base.getindex(A::gcmarray{T, N, Array{T,2}}, I::Vararg{Union{Int,Array{Int},AbstractUnitRange,Colon}, N}) where {T,N}
  J=1:length(A.fIndex)
  !isa(I[1],Colon) ? J=J[I[1]] : nothing
  nFaces=length(J)

  tmpf=A.f[I...]
  if isa(tmpf,Array{eltype(A),2})
    tmp=tmpf
  else
    n3=Int(length(tmpf)/nFaces)
    K=(A.grid,eltype(A),A.fSize[J],A.fIndex[J])
    n3>1 ? tmp=gcmarray(K...,n3) : tmp=gcmarray(K...)
    for I in eachindex(tmpf); tmp.f[I] = view(tmpf[I],:,:); end
  end

  return tmp
end

"""
    getindexetc(A::gcmarray, I::Vararg{_}) where {T,N}

Same as getindex but also returns the face size and index
"""
function getindexetc(A::gcmarray{T, N,InnerArray{T,2}}, I::Vararg{Union{Int,Array{Int},AbstractUnitRange,Colon}, N}) where {T,N}
    f=A[I...]
    fSize=A.fSize[I[1]]
    fIndex=A.fIndex[I[1]]
    return f,fSize,fIndex
end
# -

function Base.setindex!(A::gcmarray{T, N,InnerArray{T,2}}, v, I::Vararg{Int, N}) where {T,N}
  return (A.f[I...] = v)
end

function Base.view(A::gcmarray{T, N,InnerArray{T,2}}, I::Vararg{Union{Int,AbstractUnitRange,Colon}, N}) where {T,N}
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

function Base.show(io::IO, z::gcmarray{T, N, Array{T,2}}) where {T,N}
    if ~isa(z.meta.unit,Missing)
      printstyled(io, "  name        = ",color=:normal)
      printstyled(io, "$(z.meta.name)\n",color=:blue)
      printstyled(io, "  unit        = ",color=:normal)
      printstyled(io, "$(z.meta.unit)\n",color=:blue)
    end
    printstyled(io, "  data type   = ",color=:normal)
    printstyled(io, "$(eltype(z))\n",color=:blue)
    printstyled(io, "  cell pos.   = ",color=:normal)
    printstyled(io, "$(z.meta.position)\n",color=:blue)
    printstyled(io, "  tile array  = ",color=:normal)
    printstyled(io, "$(size(z))\n",color=:cyan)
    printstyled(io, "  tile sizes  = ",color=:normal)
    printstyled(io, "$(size(z[1]))\n",color=:cyan)
    for iFace=2:length(z.fIndex)
      printstyled(io, "                ",color=:normal)
      printstyled(io, "$(size(z[iFace]))\n",color=:cyan)
    end
    printstyled(io, "  grid class  = ",color=:normal)
    printstyled(io, "$(z.grid.class)\n",color=:green)
    printstyled(io, "  MeshArray   = ",color=:normal)
    printstyled(io, "gcmarray \n",color=:green)
    printstyled(io, "  version     = ",color=:normal)
    printstyled(io, "$(z.version) \n",color=:green)
  return
end

function Base.similar(A::gcmarray;m::varmeta=defaultmeta)
    if ndims(A)==1
        B=gcmarray(A.grid,eltype(A),A.fSize,A.fIndex; meta=m)
    else
        B=gcmarray(A.grid,eltype(A),A.fSize,A.fIndex,size(A,2); meta=m)
    end
    return B
end

# ### Customize broadcasting

Base.BroadcastStyle(::Type{<:gcmarray}) = Broadcast.ArrayStyle{gcmarray}()

function Base.similar(bc::Broadcast.Broadcasted{Broadcast.ArrayStyle{gcmarray}}, ::Type{ElType}) where ElType
  # Scan the inputs for the gcmarray:
  A = find_gcmarray(bc)
  # Create the gcmarray output:
  if ndims(A)==1
        B=gcmarray(A.grid,ElType,A.fSize,A.fIndex)
  else
        B=gcmarray(A.grid,ElType,A.fSize,A.fIndex,size(A,2))
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
        @inbounds dest[I] = gcmarray_getindex_evalf(bc′,I)
    end
    return dest
end

function gcmarray_getindex_evalf(bc,I)
  @boundscheck checkbounds(bc, I)
  args = Broadcast._getindex(bc.args, I)
  return bc.f.(args...)
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

#include("Type_gcmvector.jl");
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

#include("Grids.jl");
"""
    simple_periodic_domain(np::Integer,nq=missing)

Set up a simple periodic domain of size np x nq

```jldoctest
using MeshArrays
np=16 #domain size is np x np
Γ=simple_periodic_domain(np)
isa(Γ["XC"],MeshArray)

# output

true
```

"""
function simple_periodic_domain(np::Integer,nq=missing)
    ismissing(nq) ? nq=np : nothing

    nFaces=1
    ioSize=[np nq]
    facesSize=Array{NTuple{2, Int},1}(undef,nFaces)
    facesSize[:].=[(np,nq)]
    ioPrec=Float32
    γ=gcmgrid("","PeriodicDomain",1,facesSize, ioSize, ioPrec, read, write)

    Γ=Dict()
    list_n=("XC","XG","YC","YG","RAC","RAW","RAS","RAZ","DXC","DXG","DYC","DYG","Depth","hFacC","hFacS","hFacW");
    list_u=(u"m",u"m",u"m",u"m",u"m^2",u"m^2",u"m^2",u"m^2",u"m",u"m",u"m",u"m",u"m",1.0,1.0,1.0)
    pc=fill(0.5,2); pg=fill(0.0,2); pu=[0.,0.5]; pv=[0.5,0.];
    list_p=(pc,pg,pc,pg,pc,pu,pv,pg,pu,pv,pv,pu,pc,fill(0.5,3),[0.,0.5,0.5],[0.5,0.,0.5])
    for ii=1:length(list_n);
        tmp1=fill(1.,np,nq*nFaces);
        m=varmeta(list_u[ii],list_p[ii],list_n[ii],list_n[ii]);
        tmp1=γ.read(tmp1,MeshArray(γ,Float64;meta=m));
        tmp2=Symbol(list_n[ii]);
        @eval (($tmp2) = ($tmp1))
        Γ[list_n[ii]]=tmp1
    end

    Γ["XC"][1]=vec(0.5:1.0:np-0.5)*ones(1,nq)
    Γ["XG"][1]=vec(0.0:1.0:np-1.0)*ones(1,nq)
    Γ["YC"][1]=ones(np,1)*transpose(vec(0.5:1.0:nq-0.5))
    Γ["YG"][1]=ones(np,1)*transpose(vec(0.0:1.0:nq-1.0))
    return Γ
end

## GridSpec function with default GridName argument:

GridSpec() = GridSpec("PeriodicDomain","./")

## GridSpec function with GridName argument:

"""
    GridSpec(GridName,GridParentDir="./")

Select one of the pre-defined grids (by `GridName`) and return 
the corresponding `gmcgrid` -- a global grid specification 
which contains the grid files location (`GridParentDir`).
    

Possible choices for `GridName`:

- `"PeriodicDomain"`
- `"PeriodicChannel"`
- `"CubeSphere"`
- `"LatLonCap"``

```jldoctest
using MeshArrays
g = GridSpec()
g = GridSpec("PeriodicChannel",MeshArrays.GRID_LL360)
g = GridSpec("CubeSphere",MeshArrays.GRID_CS32)
g = GridSpec("LatLonCap",MeshArrays.GRID_LLC90)
isa(g,gcmgrid)

# output

true
```
"""
function GridSpec(GridName,GridParentDir="./")

grDir=GridParentDir
if GridName=="LatLonCap"
    nFaces=5
    grTopo="LatLonCap"
    ioSize=[90 1170]
    facesSize=[(90, 270), (90, 270), (90, 90), (270, 90), (270, 90)]
    ioPrec=Float64
elseif GridName=="CubeSphere"
    nFaces=6
    grTopo="CubeSphere"
    ioSize=[32 192]
    facesSize=[(32, 32), (32, 32), (32, 32), (32, 32), (32, 32), (32, 32)]
    ioPrec=Float32
elseif GridName=="PeriodicChannel"
    nFaces=1
    grTopo="PeriodicChannel"
    ioSize=[360 160]
    facesSize=[(360, 160)]
    ioPrec=Float32
elseif GridName=="PeriodicDomain"
    nFaces=4
    grTopo="PeriodicDomain"
    ioSize=[80 42]
    facesSize=[(40, 21), (40, 21), (40, 21), (40, 21)]
    ioPrec=Float32
else
    error("unknown GridName case")
end

return gcmgrid(grDir,grTopo,nFaces,facesSize, ioSize, ioPrec, read, write)

end

## GridLoad function

"""
    GridLoad(γ::gcmgrid)

Return a `Dict` of grid variables read from files located in `γ.path` (see `?GridSpec`).

Based on the MITgcm naming convention, grid variables are:

- XC, XG, YC, YG, AngleCS, AngleSN, hFacC, hFacS, hFacW, Depth.
- RAC, RAW, RAS, RAZ, DXC, DXG, DYC, DYG.
- DRC, DRF, RC, RF (one-dimensional)

```jldoctest
using MeshArrays
γ = GridSpec("CubeSphere",MeshArrays.GRID_CS32)
Γ = GridLoad(γ)

isa(Γ["XC"],MeshArray)

# output

true
```
"""
function GridLoad(γ::gcmgrid)

    Γ=Dict()

    list_n=("XC","XG","YC","YG","RAC","RAW","RAS","RAZ","DXC","DXG","DYC","DYG","Depth");
    list_u=(u"°",u"°",u"°",u"°",u"m^2",u"m^2",u"m^2",u"m^2",u"m",u"m",u"m",u"m",u"m")
    pc=fill(0.5,2); pg=fill(0.0,2); pu=[0.,0.5]; pv=[0.5,0.];
    list_p=(pc,pg,pc,pg,pc,pu,pv,pg,pu,pv,pv,pu,pc)

    if !isempty(filter(x -> occursin("AngleCS",x), readdir(γ.path)))
        list_n=(list_n...,"AngleCS","AngleSN");
        list_u=(list_u...,1.0,1.0)
        list_p=(list_p...,pc,pc)
    end

    for ii=1:length(list_n)
        m=varmeta(list_u[ii],list_p[ii],list_n[ii],list_n[ii])
        tmp1=γ.read(γ.path*list_n[ii]*".data",MeshArray(γ,γ.ioPrec;meta=m))
        tmp2=Symbol(list_n[ii])
        @eval (($tmp2) = ($tmp1))
        Γ[list_n[ii]]=tmp1
    end

    γ.ioPrec==Float64 ? reclen=8 : reclen=4

    list_n=("DRC","DRF","RC","RF")
    for ii=1:length(list_n)
        fil=γ.path*list_n[ii]*".data"
        tmp1=stat(fil)
        n3=Int64(tmp1.size/reclen)

        fid = open(fil)
        tmp1 = Array{γ.ioPrec,1}(undef,n3)
        read!(fid,tmp1)
        tmp1 = hton.(tmp1)

        tmp2=Symbol(list_n[ii])
        @eval (($tmp2) = ($tmp1))
        Γ[list_n[ii]]=tmp1
    end

    list_n=("hFacC","hFacS","hFacW");
    list_u=(1.0,1.0,1.0)
    list_p=(fill(0.5,3),[0.,0.5,0.5],[0.5,0.,0.5])
    n3=length(Γ["RC"])
    for ii=1:length(list_n)
        m=varmeta(list_u[ii],list_p[ii],list_n[ii],list_n[ii]);
        tmp1=γ.read(γ.path*list_n[ii]*".data",MeshArray(γ,γ.ioPrec,n3;meta=m))
        tmp2=Symbol(list_n[ii])
        @eval (($tmp2) = ($tmp1))
        Γ[list_n[ii]]=tmp1
    end

    return Γ

end

"""
    GridOfOnes(grTp,nF,nP)

Define all-ones grid variables instead of using `GridSpec` & `GridLoad`. E.g.

```
γ,Γ=GridOfOnes("CubeSphere",6,20);
```
"""
function GridOfOnes(grTp,nF,nP)

    grDir=""
    grTopo=grTp
    nFaces=nF
    if grTopo=="LatLonCap"
        ioSize=[nP nP*nF]
    elseif grTopo=="CubeSphere"
        ioSize=[nP nP*nF]
    elseif grTopo=="PeriodicChannel"
        ioSize=[nP nP]
    elseif grTopo=="PeriodicDomain"
        nFsqrt=Int(sqrt(nF))
        ioSize=[nP*nFsqrt nP*nFsqrt]
    end
    facesSize=Array{NTuple{2, Int},1}(undef,nFaces)
    facesSize[:].=[(nP,nP)]
    ioPrec=Float32

    γ=gcmgrid(grDir,grTopo,nFaces,facesSize, ioSize, ioPrec, read, write)

    Γ=Dict()
    list_n=("XC","XG","YC","YG","RAC","RAW","RAS","RAZ","DXC","DXG","DYC","DYG","Depth","hFacC","hFacS","hFacW");
    list_u=(u"m",u"m",u"m",u"m",u"m^2",u"m^2",u"m^2",u"m^2",u"m",u"m",u"m",u"m",u"m",1.0,1.0,1.0)
    pc=fill(0.5,2); pg=fill(0.0,2); pu=[0.,0.5]; pv=[0.5,0.];
    list_p=(pc,pg,pc,pg,pc,pu,pv,pg,pu,pv,pv,pu,pc,fill(0.5,3),[0.,0.5,0.5],[0.5,0.,0.5])

    for ii=1:length(list_n);
        tmp1=fill(1.,nP,nP*nF); m=varmeta(list_u[ii],list_p[ii],list_n[ii],list_n[ii]);
        tmp1=γ.read(tmp1,MeshArray(γ,Float64;meta=m));
        tmp2=Symbol(list_n[ii]);
        @eval (($tmp2) = ($tmp1))
        Γ[list_n[ii]]=tmp1
    end

    return γ, Γ

end

"""
    Tiles(γ::gcmgrid,ni::Int,nj::Int)

Define sudomain `tiles` of size `ni,nj`. Each tile is defined by a `Dict` where
`tile,face,i,j` correspond to tile ID, face ID, index ranges.

```jldoctest
using MeshArrays
γ=GridSpec("LatLonCap",MeshArrays.GRID_LLC90)
τ=Tiles(γ,30,30)

isa(τ[1],Dict)

# output

true
```
"""
function Tiles(γ::gcmgrid,ni::Int,nj::Int)
    nt=Int(prod(γ.ioSize)/ni/nj)
    τ=Array{Dict,1}(undef,nt)
    #
    cnt=0
    for iF=1:γ.nFaces
        for jj=Int.(1:γ.fSize[iF][2]/nj)
            for ii=Int.(1:γ.fSize[iF][1]/ni)
                cnt=cnt+1
                i=(1:ni).+ni*(ii-1)
                j=(1:nj).+nj*(jj-1)
                τ[cnt]=Dict("tile" => cnt, "face" => iF, "i" => i, "j" => j)
            end
        end
    end
    #
    return τ
end

"""
    Tiles(τ::Array{Dict},x::MeshArray)

Return an `Array` of tiles which cover `x` according to tile partition `τ`.

```jldoctest
using MeshArrays
γ=GridSpec("LatLonCap",MeshArrays.GRID_LLC90)
d=γ.read(γ.path*"Depth.data",MeshArray(γ,γ.ioPrec))
τ=Tiles(γ,30,30)
td=Tiles(τ,d)

isa(td[1],Array)

# output

true
```
"""
function Tiles(τ::Array{Dict},x::MeshArray)
    nt=length(τ)
    tx=Array{typeof(x[1]),1}(undef,nt)
    dn=size(x[1],1)-x.fSize[1][1]
    for ii=1:nt
        f=τ[ii]["face"]
        i0=minimum(τ[ii]["i"])
        i1=maximum(τ[ii]["i"])
        j0=minimum(τ[ii]["j"])
        j1=maximum(τ[ii]["j"])
        tx[ii]=view(x[f],i0:i1+dn,j0:j1+dn)
    end
    return tx
end

"""
    GridAddWS!(Γ::Dict)

Compute XW, YW, XS, and YS (vector field locations) from XC, YC (tracer
field locations) and add them to Γ.

```jldoctest
using MeshArrays
γ=GridSpec("LatLonCap",MeshArrays.GRID_LLC90)
Γ=GridLoad(γ)
GridAddWS!(Γ)

isa(Γ["XC"],MeshArray)

# output

true
```
"""
function GridAddWS!(Γ::Dict)

    XC=exchange(Γ["XC"])
    YC=exchange(Γ["YC"])
    nFaces=XC.grid.nFaces
    uX=XC.meta.unit
    uY=YC.meta.unit

    XW=MeshArrays.gcmarray(XC.grid,eltype(XC);meta=varmeta(uX,[0.,0.5],"XW","XW"))
    YW=MeshArrays.gcmarray(XC.grid,eltype(XC);meta=varmeta(uY,[0.,0.5],"YW","YW"))
    XS=MeshArrays.gcmarray(XC.grid,eltype(XC);meta=varmeta(uX,[0.5,0.],"XS","XS"))
    YS=MeshArrays.gcmarray(XC.grid,eltype(XC);meta=varmeta(uY,[0.5,0.],"YS","YS"))

    for ff=1:nFaces
        tmp1=XC[ff][1:end-2,2:end-1]
        tmp2=XC[ff][2:end-1,2:end-1]
        tmp2[tmp2.-tmp1.>180]=tmp2[tmp2.-tmp1.>180].-360;
        tmp2[tmp1.-tmp2.>180]=tmp2[tmp1.-tmp2.>180].+360;
        XW[ff]=(tmp1.+tmp2)./2;
       #
        tmp1=XC[ff][2:end-1,1:end-2]
        tmp2=XC[ff][2:end-1,2:end-1]
        tmp2[tmp2.-tmp1.>180]=tmp2[tmp2.-tmp1.>180].-360;
        tmp2[tmp1.-tmp2.>180]=tmp2[tmp1.-tmp2.>180].+360;
        XS[ff]=(tmp1.+tmp2)./2;
       #
        tmp1=YC[ff][1:end-2,2:end-1]
        tmp2=YC[ff][2:end-1,2:end-1]
        YW[ff]=(tmp1.+tmp2)./2;
       #
        tmp1=YC[ff][2:end-1,1:end-2]
        tmp2=YC[ff][2:end-1,2:end-1]
        YS[ff]=(tmp1.+tmp2)./2;
    end;

    Xmax=180; Xmin=-180;
    XS[findall(XS.<Xmin)]=XS[findall(XS.<Xmin)].+360;
    XS[findall(XS.>Xmax)]=XS[findall(XS.>Xmax)].-360;
    XW[findall(XW.<Xmin)]=XW[findall(XW.<Xmin)].+360;
    XW[findall(XW.>Xmax)]=XW[findall(XW.>Xmax)].-360;

    Γ["XW"]=XW
    Γ["XS"]=XS
    Γ["YW"]=YW
    Γ["YS"]=YS
    return Γ
end

#include("Operations.jl");

## gradient methods

"""
    gradient(inFLD::MeshArray,Γ::Dict)

Compute spatial derivatives. Other methods:
```
gradient(inFLD::MeshArray,Γ::Dict,doDIV::Bool)
gradient(inFLD::MeshArray,iDXC::MeshArray,iDYC::MeshArray)
```
"""
function gradient(inFLD::MeshArray,Γ::Dict)
(dFLDdx, dFLDdy)=gradient(inFLD,Γ,true)
return dFLDdx, dFLDdy
end

function gradient(inFLD::MeshArray,Γ::Dict,doDIV::Bool)

exFLD=exchange(inFLD,1)
dFLDdx=similar(inFLD)
dFLDdy=similar(inFLD)

for a=1:inFLD.grid.nFaces
  (s1,s2)=size(exFLD.f[a])
  tmpA=view(exFLD.f[a],2:s1-1,2:s2-1)
  tmpB=tmpA-view(exFLD.f[a],1:s1-2,2:s2-1)
  tmpC=tmpA-view(exFLD.f[a],2:s1-1,1:s2-2)
  if doDIV
    dFLDdx.f[a]=tmpB./Γ["DXC"].f[a]
    dFLDdy.f[a]=tmpC./Γ["DYC"].f[a]
  else
    dFLDdx.f[a]=tmpB
    dFLDdy.f[a]=tmpC
  end
end

return dFLDdx, dFLDdy
end

function gradient(inFLD::MeshArray,iDXC::MeshArray,iDYC::MeshArray)

exFLD=exchange(inFLD,1)
dFLDdx=similar(inFLD)
dFLDdy=similar(inFLD)

for a=1:inFLD.grid.nFaces
  (s1,s2)=size(exFLD.f[a])
  tmpA=view(exFLD.f[a],2:s1-1,2:s2-1)
  tmpB=tmpA-view(exFLD.f[a],1:s1-2,2:s2-1)
  tmpC=tmpA-view(exFLD.f[a],2:s1-1,1:s2-2)
  dFLDdx.f[a]=tmpB.*iDXC.f[a]
  dFLDdy.f[a]=tmpC.*iDYC.f[a]
end

return dFLDdx, dFLDdy
end

## mask methods

function mask(fld::MeshArray)
fldmsk=mask(fld,NaN)
return fldmsk
end

"""
    mask(fld::MeshArray, val::Number)

Replace non finite values with val. Other methods:
```
mask(fld::MeshArray)
mask(fld::MeshArray, val::Number, noval::Number)
```
"""
function mask(fld::MeshArray, val::Number)
  fldmsk=similar(fld)
  for a=1:fld.grid.nFaces
    tmp1=copy(fld.f[a])
    replace!(x -> !isfinite(x) ? val : x, tmp1 )
    fldmsk.f[a]=tmp1
  end
  return fldmsk
end

function mask(fld::MeshArray, val::Number, noval::Number)
  fldmsk=similar(fld)
  for a=1:fld.grid.nFaces
    tmp1=copy(fld.f[a])
    replace!(x -> x==noval ? val : x, tmp1  )
    fldmsk.f[a]=tmp1
  end
  return fldmsk
end

## convergence methods

"""
    convergence(uFLD::MeshArray,vFLD::MeshArray)

Compute convergence of a vector field
"""
function convergence(uFLD::MeshArray,vFLD::MeshArray)

#important note:
#  Normally uFLD, vFLD should not contain any NaN;
#  if otherwise then something this may be needed:
#  uFLD=mask(uFLD,0.0); vFLD=mask(vFLD,0.0);

CONV=similar(uFLD)

(tmpU,tmpV)=exch_UV(uFLD,vFLD)
for a=1:tmpU.grid.nFaces
  (s1,s2)=size(uFLD.f[a])
  tmpU1=view(tmpU.f[a],1:s1,1:s2)
  tmpU2=view(tmpU.f[a],2:s1+1,1:s2)
  tmpV1=view(tmpV.f[a],1:s1,1:s2)
  tmpV2=view(tmpV.f[a],1:s1,2:s2+1)
  CONV.f[a]=tmpU1-tmpU2+tmpV1-tmpV2
end

return CONV
end

## smooth function

"""
    smooth(FLD::MeshArray,DXCsm::MeshArray,DYCsm::MeshArray,Γ::Dict)

Smooth out scales below DXCsm / DYCsm via diffusion
"""
function smooth(FLD::MeshArray,DXCsm::MeshArray,DYCsm::MeshArray,Γ::Dict)

#important note:
#input FLD should be land masked (NaN/1) by caller if needed

#get land masks (NaN/1):
mskC=fill(1.0,FLD) + 0.0 * mask(FLD)
(mskW,mskS)=gradient(FLD,Γ,false)
mskW=fill(1.0,FLD) + 0.0 * mask(mskW)
mskS=fill(1.0,FLD) + 0.0 * mask(mskS)

#replace NaN with 0. in FLD and land masks:
FLD=mask(FLD,0.0)
mskC=mask(mskC,0.0)
mskW=mask(mskW,0.0)
mskS=mask(mskS,0.0)

#get inverse grid spacing:
iDXC=similar(FLD)
iDYC=similar(FLD)
for a=1:FLD.grid.nFaces
  iDXC.f[a]=1.0./Γ["DXC"].f[a]
  iDYC.f[a]=1.0./Γ["DYC"].f[a]
end

#Before scaling the diffusive operator ...
tmp0=DXCsm*iDXC*mskW;
tmp00=maximum(tmp0);
tmp0=DYCsm*iDYC*mskS;
tmp00=max(tmp00,maximum(tmp0));

#... determine a suitable time period:
nbt=ceil(1.1*2*tmp00^2);
dt=1.;
T=nbt*dt;
#println("nbt="*"$nbt")

#diffusion operator times DYG / DXG:
KuxFac=mskW*DXCsm*DXCsm/T/2.0*Γ["DYG"];
KvyFac=mskS*DYCsm*DYCsm/T/2.0*Γ["DXG"];

#time steping factor:
dtFac=dt*mskC/Γ["RAC"];

#loop:
for it=1:nbt
  (dTdxAtU,dTdyAtV)=gradient(FLD,iDXC,iDYC);
  tmpU=similar(FLD)
  tmpV=similar(FLD)
  for a=1:FLD.grid.nFaces
      tmpU.f[a]=dTdxAtU.f[a].*KuxFac.f[a];
      tmpV.f[a]=dTdyAtV.f[a].*KvyFac.f[a];
  end
  tmpC=convergence(tmpU,tmpV);
  for a=1:FLD.grid.nFaces
      FLD.f[a]=FLD.f[a]-dtFac.f[a].*tmpC.f[a];
  end
end

#Apply land mask (NaN/1) to end result:
mskC=mask(mskC,NaN,0.0)
FLD=mskC*FLD

return FLD

end

## ThroughFlow function

"""
    ThroughFlow(VectorField,IntegralPath,Γ::Dict)

Compute transport through an integration path
"""
function ThroughFlow(VectorField,IntegralPath,Γ::Dict)

    #Note: vertical intergration is not always wanted; left for user to do outside

    U=VectorField["U"]
    V=VectorField["V"]

    nd=ndims(U)
    #println("nd=$nd and d=$d")

    n=fill(1,4)
    tmp=size(U)
    n[1:nd].=tmp[1:nd]

    haskey(VectorField,"factors") ? f=VectorField["factors"] : f=Array{String,1}(undef,0)
    haskey(VectorField,"dimensions") ? d=VectorField["dimensions"] : d=Array{String,1}(undef,nd)

    #a bit of a hack to distinguish gcmfaces v gcmarray indexing:
    isdefined(U,:fIndex) ? ndoffset=1 : ndoffset=0
    length(d)!=nd+ndoffset ? error("inconsistent specification of dims") : nothing

    trsp=Array{Float64}(undef,1,n[3],n[4])
    do_dz=sum(f.=="dz")
    do_dxory=sum(f.=="dxory")

    for i3=1:n[3]
        #method 1: quite slow
        #mskW=IntegralPath["mskW"]
        #do_dxory==1 ? mskW=mskW*Γ["DYG"] : nothing
        #do_dz==1 ? mskW=Γ["DRF"][i3]*mskW : nothing
        #mskS=IntegralPath["mskS"]
        #do_dxory==1 ? mskS=mskS*Γ["DXG"] : nothing
        #do_dz==1 ? mskS=Γ["DRF"][i3]*mskS : nothing
        #
        #method 2: less slow
        tabW=IntegralPath["tabW"]
        tabS=IntegralPath["tabS"]
        for i4=1:n[4]
            #method 1: quite slow
            #trsp[1,i3,i4]=sum(mskW*U[:,:,i3,i4])+sum(mskS*V[:,:,i3,i4])
            #
            #method 2: less slow
            trsp[1,i3,i4]=0.0
            for k=1:size(tabW,1)
                (a,i1,i2,w)=tabW[k,:]
                do_dxory==1 ? w=w*Γ["DYG"].f[a][i1,i2] : nothing
                do_dz==1 ? w=w*Γ["DRF"][i3] : nothing
                isdefined(U,:fIndex) ? u=U.f[a,i3,i4][i1,i2] : u=U.f[a][i1,i2,i3,i4]
                trsp[1,i3,i4]=trsp[1,i3,i4]+w*u
            end
            for k=1:size(tabS,1)
                (a,i1,i2,w)=tabS[k,:]
                do_dxory==1 ? w=w*Γ["DXG"].f[a][i1,i2] : nothing
                do_dz==1 ? w=w*Γ["DRF"][i3] : nothing
                isdefined(V,:fIndex) ? v=V.f[a,i3,i4][i1,i2] : v=V.f[a][i1,i2,i3,i4]
                trsp[1,i3,i4]=trsp[1,i3,i4]+w*v
            end
        end
    end

    nd+ndoffset<4 ? trsp=dropdims(trsp,dims=3) : nothing
    nd+ndoffset<3 ? trsp=dropdims(trsp,dims=2) : nothing
    nd+ndoffset==2 ? trsp=trsp[1] : nothing

    return trsp
end

## LatitudeCircles function

"""
    LatitudeCircles(LatValues,Γ::Dict)

Compute integration paths that follow latitude circles
"""
function LatitudeCircles(LatValues,Γ::Dict)

    LatitudeCircles=Array{Dict}(undef,length(LatValues))

    for j=1:length(LatValues)
        mskCint=1*(Γ["YC"] .>= LatValues[j])
        mskC=similar(mskCint)
        mskW=similar(mskCint)
        mskS=similar(mskCint)

        mskCint=exchange(mskCint,1)

        for i=1:mskCint.grid.nFaces
            tmp1=mskCint.f[i]
            # tracer mask:
            tmp2=tmp1[2:end-1,1:end-2]+tmp1[2:end-1,3:end]+
            tmp1[1:end-2,2:end-1]+tmp1[3:end,2:end-1]
            mskC.f[i]=1((tmp2.>0).&(tmp1[2:end-1,2:end-1].==0))
            # velocity masks:
            mskW.f[i]=tmp1[2:end-1,2:end-1] - tmp1[1:end-2,2:end-1]
            mskS.f[i]=tmp1[2:end-1,2:end-1] - tmp1[2:end-1,1:end-2]
        end

        function MskToTab(msk::MeshArray)
          n=Int(sum(msk .!= 0)); k=0
          tab=Array{Int,2}(undef,n,4)
          for i=1:msk.grid.nFaces
            a=msk.f[i]
            b=findall( a .!= 0)
            for ii in eachindex(b)
              k += 1
              tab[k,:]=[i,b[ii][1],b[ii][2],a[b[ii]]]
            end
          end
          return tab
        end

        LatitudeCircles[j]=Dict("lat"=>LatValues[j],
        "tabC"=>MskToTab(mskC),"tabW"=>MskToTab(mskW),"tabS"=>MskToTab(mskS))
    end

    return LatitudeCircles

end
##

#include("Exchanges.jl");

## This file contains the exchange and exch_UV functions
# along with grid-specific methods (exch_T_N.jl, etc.)

## User Front Ends

"""
    exchange(fld::MeshArray)

Exchange / transfer data between neighboring arrays. Other methods are

    exchange(fld::MeshArray,N::Integer)
    exchange(u::MeshArray,v::MeshArray)
    exchange(u::MeshArray,v::MeshArray,N::Integer)
"""
function exchange(fld::MeshArray)
  FLD=exch_T_N(fld,1);
end

function exchange(fld::MeshArray,N::Integer)
  FLD=exch_T_N(fld,N);
end

function exchange(u::MeshArray,v::MeshArray)
  (uex,vex)=exch_UV_N(u,v,1);
end

function exchange(u::MeshArray,v::MeshArray,N::Integer)
  (uex,vex)=exch_UV_N(u,v,N);
end

## dispatch over grid types

#note: the "CubeSphere" implementation covers both cs and llc

function exch_T_N(fld,N)

if fld.grid.class=="LatLonCap"
  FLD=exch_T_N_cs(fld,N)
elseif fld.grid.class=="CubeSphere"
  FLD=exch_T_N_cs(fld,N)
elseif fld.grid.class=="PeriodicChannel"
  FLD=exch_T_N_ll(fld,N)
elseif fld.grid.class=="PeriodicDomain"
  FLD=exch_T_N_dpdo(fld,N)
else
  error("unknown grid.class case")
end

return FLD

end

function exch_UV_N(u,v,N)

if u.grid.class=="LatLonCap"
  (uex,vex)=exch_UV_N_cs(u,v,N)
elseif u.grid.class=="CubeSphere"
  (uex,vex)=exch_UV_N_cs(u,v,N)
elseif u.grid.class=="PeriodicChannel"
  (uex,vex)=exch_UV_N_ll(u,v,N)
elseif u.grid.class=="PeriodicDomain"
  (uex,vex)=exch_UV_N_dpdo(u,v,N)
else
  error("unknown grid.class case")
end

return uex,vex

end

function exch_UV(u,v)

if u.grid.class=="LatLonCap"
  (uex,vex)=exch_UV_cs(u,v)
elseif u.grid.class=="CubeSphere"
  (uex,vex)=exch_UV_cs(u,v)
elseif u.grid.class=="PeriodicChannel"
  (uex,vex)=exch_UV_ll(u,v)
elseif u.grid.class=="PeriodicDomain"
  (uex,vex)=exch_UV_dpdo(u,v)
else
  error("unknown grid.class case")
end

return uex,vex

end

## Grid-specific implementations: PeriodicDomain case

function exch_T_N_dpdo(fld::MeshArray,N::Integer)

fillval=0.0

ni,nj=Int.(transpose(fld.grid.ioSize)./fld.grid.fSize[1])
s=fld.fSize
FLD=similar(fld;m=fld.meta);

for i=1:ni
  for j=1:nj
    k=i+ni*(j-1)
    kS=j-1; kS==0 ? kS=nj : nothing; kS=i+ni*(kS-1)
    kN=j+1; kN==nj+1 ? kN=1 : nothing; kN=i+ni*(kN-1)
    kW=i-1; kW==0 ? kW=ni : nothing; kW=kW+ni*(j-1)
    kE=i+1; kE==ni+1 ? kE=1 : nothing; kE=kE+ni*(j-1)

    #step 1

    FLD.f[k]=fill(fillval,s[k].+2N);
    @views FLD.f[k][N+1:N+s[k][1],N+1:N+s[k][2]]=fld.f[k];

    #step 2

    iW=(s[k][1]-N+1:s[k][1],1:s[k][2]);
    iE=(1:N,1:s[k][2]);
    jW=(1:N,N+1:N+s[k][2]);
    jE=(N+1+s[k][1]:2N+s[k][1],N+1:N+s[k][2]);
    FLD.f[k][jW[1],jW[2]]=view(fld.f[kW],iW[1],iW[2])
    FLD.f[k][jE[1],jE[2]]=view(fld.f[kE],iE[1],iE[2])

    #step 3

    iS=(1:s[k][1],s[k][2]-N+1:s[k][2]);
    iN=(1:s[k][1],1:N);
    jS=(N+1:N+s[k][1],1:N);
    jN=(N+1:N+s[k][1],N+1+s[k][2]:2N+s[k][2]);
    FLD.f[k][jS[1],jS[2]]=view(fld.f[kS],iS[1],iS[2])
    FLD.f[k][jN[1],jN[2]]=view(fld.f[kN],iN[1],iN[2])

  end
end

return FLD

end

##

function exch_UV_N_dpdo(fldU,fldV,N);

FLDU=exch_T_N_dpdo(fldU,N);
FLDV=exch_T_N_dpdo(fldV,N);

return FLDU,FLDV

end

##

function exch_UV_dpdo(fldU,fldV);

fillval=0.0

ni,nj=Int.(transpose(fldU.grid.ioSize)./fldU.grid.fSize[1])
s=fldU.fSize
FLDU=similar(fldU;m=fldU.meta)
FLDV=similar(fldV;m=fldV.meta)

for i=1:ni
  for j=1:nj
    k=i+ni*(j-1)
    kS=j-1; kS==0 ? kS=nj : nothing; kS=i+ni*(kS-1)
    kN=j+1; kN==nj+1 ? kN=1 : nothing; kN=i+ni*(kN-1)
    kW=i-1; kW==0 ? kW=ni : nothing; kW=kW+ni*(j-1)
    kE=i+1; kE==ni+1 ? kE=1 : nothing; kE=kE+ni*(j-1)

    #step 1

    FLDU.f[k]=fill(fillval,s[k][1]+1,s[k][2]);
    FLDV.f[k]=fill(fillval,s[k][1],s[k][2]+1);
    @views FLDU.f[k][1:s[k][1],1:s[k][2]]=fldU.f[k];
    @views FLDV.f[k][1:s[k][1],1:s[k][2]]=fldV.f[k];

    #step 2

    FLDU.f[k][s[k][1]+1,1:s[k][2]]=view(fldU.f[kE],1,1:s[k][2])

    #step 3

    FLDV.f[k][1:s[k][1],s[k][2]+1]=view(fldV.f[kN],1:s[k][1],1)

  end
end

return FLDU,FLDV

end

## Grid-specific implementations: PeriodicChannel case

function exch_T_N_ll(fld::MeshArray,N::Integer)

fillval=0.0

#step 1

s=size.(fld.f);
FLD=similar(fld;m=fld.meta);
FLD.f[1]=fill(fillval,s[1].+2N);
@views FLD.f[1][N+1:N+s[1][1],N+1:N+s[1][2]]=fld.f[1];

#step 2

iW=(s[1][1]-N+1:s[1][1],1:s[1][2]);
iE=(1:N,1:s[1][2]);
jW=(1:N,N+1:N+s[1][2]);
jE=(N+1+s[1][1]:2N+s[1][1],N+1:N+s[1][2]);
FLD.f[1][jW[1],jW[2]]=view(fld.f[1],iW[1],iW[2])
FLD.f[1][jE[1],jE[2]]=view(fld.f[1],iE[1],iE[2])

return FLD

end

##

function exch_UV_N_ll(fldU,fldV,N);

FLDU=exch_T_N_ll(fldU,N);
FLDV=exch_T_N_ll(fldV,N);

return FLDU,FLDV

end

##

function exch_UV_ll(fldU,fldV);

fillval=0.0

#step 1

s=size.(fldU.f);
FLDU=similar(fldU;m=fldU.meta);
FLDV=similar(fldV;m=fldV.meta);

FLDU.f[1]=fill(fillval,s[1][1]+1,s[1][2]);
FLDV.f[1]=fill(fillval,s[1][1],s[1][2]+1);
@views FLDU.f[1][1:s[1][1],1:s[1][2]]=fldU.f[1];
@views FLDV.f[1][1:s[1][1],1:s[1][2]]=fldV.f[1];

#step 2

FLDU.f[1][s[1][1]+1,1:s[1][2]]=view(fldU.f[1],1,1:s[1][2])

return FLDU,FLDV

end

## Grid-specific implementations: CubeSphere & LatLonCap case

#note: the "CubeSphere" implementation covers both cs and llc

function exch_T_N_cs(fld::MeshArray,N::Integer)

fillval=0.0

#step 1

s=size.(fld.f)
nf=fld.grid.nFaces
nf==5 ? s=vcat(s,s[3]) : nothing
tp=fld.grid.class
FLD=similar(fld;m=fld.meta)

for i=1:nf; FLD.f[i]=fill(fillval,s[i].+2N); end;
#code below yields strange, seemingly incorrect results:
#for i=1:nf; FLD.f[i]=Array{eltype(fld.f[i])}(undef,s[i].+2N); end;

#all versions below yield same @time and memory (despite diff in allocs)
for i=1:nf;
# FLD.f[i][N+1:end-N,N+1:end-N]=fld.f[i];
 @views FLD.f[i][N+1:N+s[i][1],N+1:N+s[i][2]]=fld.f[i];
end;

#step 2

(ovfW,ovfE,ovfS,ovfN,evfW,evfE,evfS,evfN)=exch_cs_viewfunctions();

for a=1:nf
(jW, jE, jS, jN)=exch_cs_target(s[a],N)
(aW,aE,aS,aN,iW,iE,iS,iN)=exch_cs_sources(a,s,N)
if !iseven(a)
 aW <= nf ? FLD.f[a][jW[1],jW[2]]=ovfW(fld.f[aW],iW[1],iW[2]) : nothing
 aE <= nf ? FLD.f[a][jE[1],jE[2]]=ovfE(fld.f[aE],iE[1],iE[2]) : nothing
 aS <= nf ? FLD.f[a][jS[1],jS[2]]=ovfS(fld.f[aS],iS[1],iS[2]) : nothing
 aN <= nf ? FLD.f[a][jN[1],jN[2]]=ovfN(fld.f[aN],iN[1],iN[2]) : nothing
else
 aW <= nf ? FLD.f[a][jW[1],jW[2]]=evfW(fld.f[aW],iW[1],iW[2]) : nothing
 aE <= nf ? FLD.f[a][jE[1],jE[2]]=evfE(fld.f[aE],iE[1],iE[2]) : nothing
 aS <= nf ? FLD.f[a][jS[1],jS[2]]=evfS(fld.f[aS],iS[1],iS[2]) : nothing
 aN <= nf ? FLD.f[a][jN[1],jN[2]]=evfN(fld.f[aN],iN[1],iN[2]) : nothing
end
end

return FLD

end

##

function exch_UV_N_cs(fldU::MeshArray,fldV::MeshArray,N::Integer)

fillval=0.0

#step 1

s=size.(fldU.f)
nf=fldU.grid.nFaces
nf==5 ? s=vcat(s,s[3]) : nothing
tp=fldU.grid.class
FLDU=similar(fldU;m=fldU.meta)
FLDV=similar(fldV;m=fldV.meta)

for i=1:nf;
 FLDU.f[i]=fill(fillval,s[i].+2N);
 FLDV.f[i]=fill(fillval,s[i].+2N);
 @views FLDU.f[i][N+1:N+s[i][1],N+1:N+s[i][2]]=fldU.f[i];
 @views FLDV.f[i][N+1:N+s[i][1],N+1:N+s[i][2]]=fldV.f[i];
end;

#step 2

(ovfW,ovfE,ovfS,ovfN,evfW,evfE,evfS,evfN)=exch_cs_viewfunctions();

for a=1:nf
(jW, jE, jS, jN)=exch_cs_target(s[a],N)
(aW,aE,aS,aN,iW,iE,iS,iN)=exch_cs_sources(a,s,N)
if !iseven(a)
 aW <= nf ? FLDU.f[a][jW[1],jW[2]]=ovfW(fldV.f[aW],iW[1],iW[2]) : nothing
 aE <= nf ? FLDU.f[a][jE[1],jE[2]]=ovfE(fldU.f[aE],iE[1],iE[2]) : nothing
 aS <= nf ? FLDU.f[a][jS[1],jS[2]]=ovfS(fldU.f[aS],iS[1],iS[2]) : nothing
 aN <= nf ? FLDU.f[a][jN[1].+1,jN[2]]=-ovfN(fldV.f[aN],iN[1],iN[2]) : nothing
 aW <= nf ? FLDV.f[a][jW[1],jW[2].+1]=-ovfW(fldU.f[aW],iW[1],iW[2]) : nothing
 aE <= nf ? FLDV.f[a][jE[1],jE[2]]=ovfE(fldV.f[aE],iE[1],iE[2]) : nothing
 aS <= nf ? FLDV.f[a][jS[1],jS[2]]=ovfS(fldV.f[aS],iS[1],iS[2]) : nothing
 aN <= nf ? FLDV.f[a][jN[1],jN[2]]=ovfN(fldU.f[aN],iN[1],iN[2]) : nothing
else
 aW <= nf ? FLDU.f[a][jW[1],jW[2]]=evfW(fldU.f[aW],iW[1],iW[2]) : nothing
 aE <= nf ? FLDU.f[a][jE[1],jE[2]]=evfE(fldV.f[aE],iE[1],iE[2]) : nothing
 aS <= nf ? FLDU.f[a][jS[1].+1,jS[2]]=-evfS(fldV.f[aS],iS[1],iS[2]) : nothing
 aN <= nf ? FLDU.f[a][jN[1],jN[2]]=evfN(fldU.f[aN],iN[1],iN[2]) : nothing
 aW <= nf ? FLDV.f[a][jW[1],jW[2]]=evfW(fldV.f[aW],iW[1],iW[2]) : nothing
 aE <= nf ? FLDV.f[a][jE[1],jE[2].+1]=-evfE(fldU.f[aE],iE[1],iE[2]) : nothing
 aS <= nf ? FLDV.f[a][jS[1],jS[2]]=evfS(fldU.f[aS],iS[1],iS[2]) : nothing
 aN <= nf ? FLDV.f[a][jN[1],jN[2]]=evfN(fldV.f[aN],iN[1],iN[2]) : nothing
end
end

return FLDU,FLDV

end

##

function exch_UV_cs(fldU::MeshArray,fldV::MeshArray)

fillval=0.0

#step 1

s=size.(fldU.f)
nf=fldU.grid.nFaces
nf==5 ? s=vcat(s,s[3]) : nothing
tp=fldU.grid.class
FLDU=similar(fldU;m=fldU.meta)
FLDV=similar(fldV;m=fldV.meta)

for i=1:nf
  FLDU.f[i]=fill(fillval,s[i][1]+1,s[i][2]);
  FLDV.f[i]=fill(fillval,s[i][1],s[i][2]+1);
  @views FLDU.f[i][1:s[i][1],1:s[i][2]]=fldU.f[i];
  @views FLDV.f[i][1:s[i][1],1:s[i][2]]=fldV.f[i];
end

 #step 2

(ovfW,ovfE,ovfS,ovfN,evfW,evfE,evfS,evfN)=exch_cs_viewfunctions();

for a=1:nf
(jW, jE, jS, jN)=exch_cs_target(s[a],1)
(aW,aE,aS,aN,iW,iE,iS,iN)=exch_cs_sources(a,s,1)
if !iseven(a)
 aE <= nf ? FLDU.f[a][jE[1].-1,jE[2].-1]=ovfE(fldU.f[aE],iE[1],iE[2]) : nothing
 aN <= nf ? FLDV.f[a][jN[1].-1,jN[2].-1]=ovfN(fldU.f[aN],iN[1],iN[2]) : nothing
else
 aE <= nf ? FLDU.f[a][jE[1].-1,jE[2].-1]=evfE(fldV.f[aE],iE[1],iE[2]) : nothing
 aN <= nf ? FLDV.f[a][jN[1].-1,jN[2].-1]=evfN(fldV.f[aN],iN[1],iN[2]) : nothing
end
end

return FLDU,FLDV

end

## Convenience functions used in the cs & llc case

function exch_cs_target(sa::Tuple{Int64,Int64},N::Integer)

    #target array indices
    jW=(1:N,N+1:N+sa[2]);
    jE=(N+1+sa[1]:2N+sa[1],N+1:N+sa[2]);
    jS=(N+1:N+sa[1],1:N);
    jN=(N+1:N+sa[1],N+1+sa[2]:2N+sa[2]);

    return jW, jE, jS, jN

end

function exch_cs_sources(a::Integer,s::Array{Tuple{Int64,Int64},1},N::Integer)

#source array IDs
aW=0; aE=0; aS=0; aN=0;
if a==1;     aW=5; aE=2; aS=6; aN=3;
elseif a==2; aW=1; aE=4; aS=6; aN=3;
elseif a==3; aW=1; aE=4; aS=2; aN=5;
elseif a==4; aW=3; aE=6; aS=2; aN=5;
elseif a==5; aW=3; aE=6; aS=4; aN=1;
elseif a==6; aW=5; aE=2; aS=4; aN=1;
else; error("Array index is out of bounds.");
end;

if !iseven(a)
    #source array indices
    iW=(1:s[aW][1],s[aW][2]-N+1:s[aW][2]);
    iE=(1:N,1:s[aE][2]);
    iS=(1:s[aS][1],s[aS][2]-N+1:s[aS][2]);
    iN=(1:N,1:s[aN][2]);
else
    #source array indices
    iW=(s[aW][1]-N+1:s[aW][1],1:s[aW][2]);
    iE=(1:s[aE][1],1:N);
    iS=(s[aS][1]-N+1:s[aS][1],1:s[aS][2]);
    iN=(1:s[aN][1],1:N);
end

return aW,aE,aS,aN,iW,iE,iS,iN
end

function exch_cs_viewfunctions()
#view functions for odd numbered arrays
ovfW(x,i,j)=PermutedDimsArray(view(x,reverse(i),j),(2,1))
ovfE(x,i,j)=view(x,i,j)
ovfS(x,i,j)=view(x,i,j)
ovfN(x,i,j)=PermutedDimsArray(view(x,i,reverse(j)),(2,1))
#view functions for even numbered arrays
evfW(x,i,j)=view(x,i,j)
evfE(x,i,j)=PermutedDimsArray(view(x,reverse(i),j),(2,1))
evfS(x,i,j)=PermutedDimsArray(view(x,i,reverse(j)),(2,1))
evfN(x,i,j)=view(x,i,j)
#
return ovfW,ovfE,ovfS,ovfN,evfW,evfE,evfS,evfN
end

#include("ReadWrite.jl");

import Base: read, write

"""
    read(fil::String,x::MeshArray)

Read binary file to MeshArray. Other methods:

```
read(xx::Array,x::MeshArray) #from Array
```
"""
function read(fil::String,x::MeshArray)

  (n1,n2)=x.grid.ioSize
  (nFaces,n3)=nFacesEtc(x)

  fid = open(fil)
  xx = Array{eltype(x),2}(undef,(n1*n2,n3))
  read!(fid,xx)
  xx = hton.(xx)
  close(fid)

  return read(xx,x)

end

function read(xx::Array,x::MeshArray)

  facesSize=x.grid.fSize
  (n1,n2)=x.grid.ioSize
  (nFaces,n3)=nFacesEtc(x)

  size(xx)!=(n1*n2,n3) ? xx=reshape(xx,(n1*n2,n3)) : nothing

  y=similar(x; m=x.meta)
  i0=0; i1=0;
  for iFace=1:nFaces
    i0=i1+1;
    nn=facesSize[iFace][1]; mm=facesSize[iFace][2];
    i1=i1+nn*mm;

    if n3>1 && ndims(x.f)==1
      y.f[iFace]=reshape(xx[i0:i1,:],(nn,mm,n3))
    elseif n3>1 && ndims(x.f)==2
      for i3=1:n3
        y.f[iFace,i3]=reshape(xx[i0:i1,i3],(nn,mm))
      end
    else
      y.f[iFace]=reshape(xx[i0:i1,:],(nn,mm))
    end
  end

  return y

end


"""
    write(fil::String,x::MeshArray)

Write MeshArray to binary file. Other methods:

```
write(xx::Array,x::MeshArray) #to Array
```
"""
function write(fil::String,x::MeshArray)
  y=write(x)
  fid = open(fil,"w")
  write(fid,ntoh.(y))
  close(fid)
end

function write(x::MeshArray)

  facesSize=x.grid.fSize
  (n1,n2)=x.grid.ioSize
  (nFaces,n3)=nFacesEtc(x)

  y = Array{eltype(x),2}(undef,(n1*n2,n3))
  i0=0; i1=0;
  for iFace=1:nFaces;
    i0=i1+1;
    nn=facesSize[iFace][1];
    mm=facesSize[iFace][2];
    i1=i1+nn*mm;
    if n3>1 && ndims(x.f)==2
      for i3=1:n3
        y[i0:i1,i3]=reshape(x.f[iFace,i3],(nn*mm,1))
      end
    else
      y[i0:i1,:]=reshape(x.f[iFace],(nn*mm,n3))
    end
  end

  y=reshape(y,(n1,n2,n3));
  n3==1 ? y=dropdims(y,dims=3) : nothing

  return y

end

#include("Solvers.jl");
# # Poisson And Laplace Equation Solvers On The Sphere
#
# Compute scalar potential, and divergent component, from vertically integrated transport.

using SparseArrays, Statistics

"""
    MaskWetPoints(TrspCon)

Mask land points with NaN.
```
(TrspCon, mskWet, mskDry)=MaskWetPoints(TrspCon)
```
"""
function MaskWetPoints(TrspCon)
    mskWet=1.0 .+ 0.0 * TrspCon
    mskDry=1.0 * isnan.(mskWet)
    mskDry=mask(mskDry,NaN,0.0)
    #
    tmp1=fill(1.0,mskWet); tmp2=exchange(tmp1);
    for I=1:size(tmp1.f,1)
        tmp3=mskWet[I]; tmp4=tmp2[I];
        tmp4=tmp4[2:end-1,1:end-2]+tmp4[2:end-1,3:end]+tmp4[1:end-2,2:end-1]+tmp4[3:end,2:end-1];
        !isempty(findall(isnan.(tmp4) .& (!isnan).(tmp3))) ? println("warning: modified mask") : nothing
        tmp3[findall(isnan.(tmp4))] .= NaN
        mskWet[I]=tmp3
    end
    #
    TrspCon=mask(TrspCon,0.0)*mskWet;
    return TrspCon, mskWet, mskDry
end

"""
    MapWetPoints(mskWet)

Mapping from global array to global ocean vector.
```
(Kvec,Lvec,Kmap,Lmap)=MapWetPoints(mskWet)
```
"""
function MapWetPoints(mskWet)
    tmp1=write(mskWet)[:]
    kk=findall((!isnan).(tmp1))
    nn=length(kk); s0=size(tmp1); s1=mskWet.grid.ioSize;
    Kvec=fill(0.0,s0...); Kvec[kk]=kk; Kmap=read(reshape(Kvec,s1...),mskWet) #global array indices
    Lvec=fill(0.0,s0...); Lvec[kk]=1:nn; Lmap=read(reshape(Lvec,s1...),mskWet) #global vector indices
    return Kvec,Lvec,Kmap,Lmap
end

"""
    SeedWetPoints(tmp::MeshArray,Kmap::MeshArray,Lmap::MeshArray,I...)

Seed a subset of grid points.
```
(FLDones,FLDkkFROM)=SeedWetPoints(tmp::MeshArray,Kmap::MeshArray,Lmap::MeshArray,I...)
```
"""
function SeedWetPoints(tmp::MeshArray,Kmap::MeshArray,Lmap::MeshArray,I...)
    aa=I[1]
    ii=I[2]
    jj=I[3]

    #1) seed points (FLDones) and neighborhood of influence (FLDkkFROM)
    FLDones=fill(0.0,tmp)
    FLDones[aa][ii:3:end,jj:3:end].=1.0
    FLDones[aa][findall(Kmap[aa].==0.0)].=0.0

    FLDkkFROMtmp=fill(0.0,tmp)
    FLDkkFROMtmp[aa][ii:3:end,jj:3:end]=Kmap[aa][ii:3:end,jj:3:end]
    FLDkkFROMtmp[aa][findall(isnan.(tmp[aa]))].=0.0

    FLDkkFROM=exchange(FLDkkFROMtmp)
    FLDkkFROM=mask(FLDkkFROM,0.0)

    for bb in 1:tmp.grid.nFaces
        tmp1=FLDkkFROM[bb]
        tmp2=zeros(size(tmp1) .- 2)
        for ii2 in 1:3
            for jj2 in 1:3
                tmp2=tmp2+tmp1[ii2:end-3+ii2,jj2:end-3+jj2]
            end
        end
        FLDkkFROM[bb]=tmp2
    end

    return FLDones,FLDkkFROM
end

"""
    MatrixForPoisson(TrspCon,mskWet,mskDry,Kvec,Lvec,Kmap,Lmap)

Assemble sparse matrix using mskWet, Kvec, Lvec directly and Kmap, Lmap via SeedWetPoints
```
A=MatrixForPoisson(TrspCon,mskWet,mskDry,Kvec,Lvec,Kmap,Lmap)
```
"""
function MatrixForPoisson(TrspCon,mskWet,mskDry,Kvec,Lvec,Kmap,Lmap)
    I=Array{Int}(undef,0)
    J=Array{Int}(undef,0)
    V=Array{Float64}(undef,0)

    for aa=1:TrspCon.grid.nFaces;
        for ii=1:3; for jj=1:3;
            #1) compute effect of each point on neighboring target point:
            (FLDones,FLDkkFROM)=SeedWetPoints(TrspCon,Kmap,Lmap,aa,ii,jj);
            (tmpU,tmpV)=gradient(FLDones,Dict(),false)
            dFLDdt=convergence(tmpU,tmpV);
            #2) mask `dFLDdt` since we use a **Neumann** boundary condition.
            #Extrapolation uses a **Dirichlet** boundary condition, so mskFreeze should not be applied then.
            isa(FLDkkFROM,MeshArray) ? FLDkkFROM=write(FLDkkFROM)[:] : nothing;
            #3.1) For wet points -- add contributions in matrix:
            dFLDdtWet=write(dFLDdt.*mskWet)[:];
            tmp1=findall( (dFLDdtWet .!= 0.0) .* (!isnan).(dFLDdtWet));
            tmpV=dFLDdtWet[tmp1]; tmpJ=FLDkkFROM[tmp1]; tmpI=Kvec[tmp1];
            I=[I;Lvec[Int.(tmpI)]]; J=[J;Lvec[Int.(tmpJ)]]; V=[V;tmpV];
            size(tmpV)
            #3.2) For dry points -- This part reflects the `Neumann` boundary condition:
            dFLDdtDry=write(dFLDdt.*mskDry)[:];
            tmp1=findall( (dFLDdtDry .!= 0.0) .* (!isnan).(dFLDdtDry));
            tmpV=dFLDdtDry[tmp1]; tmpIJ=FLDkkFROM[tmp1];
            I=[I;Lvec[Int.(tmpIJ)]]; J=[J;Lvec[Int.(tmpIJ)]]; V=[V;tmpV];
        end; end;
    end;

    nn=sum((!isnan).(mskWet))
    A=sparse(I,J,V,nn,nn)
end

"""
    ScalarPotential(TrspCon)

Scalar potential inversion.
```
TrspPot=ScalarPotential(TrspCon)
```
"""
function ScalarPotential(TrspCon)
    (TrspCon, mskWet, mskDry)=MaskWetPoints(TrspCon)
    (Kvec,Lvec,Kmap,Lmap)=MapWetPoints(mskWet)

    A=MatrixForPoisson(TrspCon,mskWet,mskDry,Kvec,Lvec,Kmap,Lmap)
    yy=write(TrspCon)[:]
    yy=yy[findall(Kvec .!= 0)]

    xx=A\yy
    xx=xx .- median(xx)

    TrspPot=0.0 * write(TrspCon)[:]
    TrspPot[findall(Kvec .!= 0)]=xx
    TrspPot=read(TrspPot,TrspCon)

    return TrspPot
end

"""
    VectorPotential(TrspX,TrspY,Γ,method::Int=1)

Vector potential inversion.
```
TrspPot=ScalarPotential(TrspCon)
```
"""
function VectorPotential(TrspX::MeshArray,TrspY::MeshArray,Γ::Dict,method::Int=1)

    # 1)  streamfunction face by face:

    (fldU,fldV)=exch_UV(TrspX,TrspY);
    fldU=mask(fldU,0.0); fldV=mask(fldV,0.0);

    psi=similar(fldV)
    for I in eachindex(fldV)
        tmp2=cumsum(fldV[I],dims=1)
        psi[I]=[zeros(1,size(tmp2,2));tmp2]
    end

    # 2a)  reset one land value per face to zero (method 1)

    if method==1
        mskW=mask(1.0 .+ 0.0 * mask(view(Γ["hFacW"],:,1),NaN,0.0),0.0)
        mskS=mask(1.0 .+ 0.0 * mask(view(Γ["hFacS"],:,1),NaN,0.0),0.0)
        (mskW,mskS)=exch_UV(mskW,mskS); mskW=abs.(mskW); mskS=abs.(mskS)

        for iF=1:TrspX.grid.nFaces
            tmp2=mskS[iF]
            tmp3a=[ones(1,size(tmp2,2));tmp2]
            tmp3b=[tmp2;ones(1,size(tmp2,2))]
            tmpB=tmp3a.*tmp3b
            #
            tmpA=psi[iF];
            I=findall(tmpB.==0)
            ii=I[1][1]; jj=I[1][2]
            tmpA[:,jj]=tmpA[:,jj] .- tmpB[ii,jj]
            #
            for kk=1:jj-1
                tmpE=tmpA[:,kk+1]-tmpA[:,kk]
                tmpE+=fldU[iF][:,kk]
                tmpE=median(tmpE[findall((!isnan).(tmpE))])
                tmpA[:,kk]=tmpA[:,kk] .- tmpE
            end
            #
            for kk=jj+1:size(tmpA,2)
                tmpE=tmpA[:,kk]-tmpA[:,kk-1]
                tmpE+=fldU[iF][:,kk-1]
                tmpE=median(tmpE[findall((!isnan).(tmpE))])
                tmpA[:,kk]=tmpA[:,kk] .- tmpE
            end;
            #
            psi[iF]=tmpA
        end
    end

    # 2b)  subtract divergent flow line by line (method 2)

    if method==2
        for I in eachindex(fldV)
            tmp2=diff(psi[I],dims=2)+fldU[I]
            tmp3=cumsum(mean(tmp2,dims=1),dims=2)
            psi[I]=psi[I]-ones(size(psi[I],1),1)*[0.0 tmp3]
        end
    end

    # 3) match edge values:

    if fldU.grid.nFaces>1
        TMP1=similar(psi)
        for I in eachindex(TMP1); TMP1[I] = fill(I,size(psi[I])); end
        TMP2=exchange(TMP1) #this is a trick

    for I in 1:TrspX.grid.nFaces-1
        tmp2=exchange(psi) #this is a trick
        tmp3=tmp2[I+1]; tmp3[3:end-2,3:end-2].=NaN #mask out interior points
        TMP3=TMP2[I+1]; tmp3[findall(TMP3.>I+1)].=NaN #mask out edges points coming from unadjusted faces
        tmp3[findall(TMP3.==0)].=NaN
        tmp3[:,1]=tmp3[:,1]-tmp3[:,2]; tmp3[:,end]=tmp3[:,end]-tmp3[:,end-1] #compare edge points
        tmp3[1,:]=tmp3[1,:]-tmp3[2,:]; tmp3[end,:]=tmp3[end,:]-tmp3[end-1,:] #compare edge points
        tmp3[2:end-1,2:end-1] .= NaN #mask out remaining interior points
        psi[I+1]=psi[I+1] .+ median(tmp3[findall((!isnan).(tmp3))]) #adjust the face data
    end
    end

    # 4) put streamfunction at cell center

    for I in 1:TrspX.grid.nFaces
        tmp3=psi[I]
        tmp3=(tmp3[:,1:end-1]+tmp3[:,2:end])/2
        tmp3=(tmp3[1:end-1,:]+tmp3[2:end,:])/2
        psi[I]=tmp3
    end

    # 5) reset one land point to zero

    #all values:
    tmp1=write(psi)

    #all land points:
    tmp2=write(view(Γ["hFacC"],:,1))
    tmp2=findall( (tmp2 .== 0.0) .& (!isnan).(tmp1))

    #closest to Boston:
    tmp_lon=write(Γ["XC"])[tmp2]
    tmp_lat=write(Γ["YC"])[tmp2]
    tmp_dis=(tmp_lat .- 42.3601).^2 + (tmp_lon .- 71.0589).^2
    tmp2=tmp2[findall(tmp_dis .== minimum(tmp_dis))]

    #set that point to zero
    psi=psi .- tmp1[tmp2[1]];

    return psi
end

#include("Interpolation.jl");
using NearestNeighbors
import NearestNeighbors: knn

"""
    knn(xgrid,ygrid::MeshArray,x,y::Array{T,1},k::Int)

Find k nearest neighbors to each point in x,y on xgrid,ygrid

```
lon=collect(0.1:0.5:2.1); lat=collect(0.1:0.5:2.1);
(f,i,j,c)=knn(Γ["XC"],Γ["YC"],lon,lat)
```
"""
function knn(xgrid::MeshArray,ygrid::MeshArray,
        xvec::Array{T,1},yvec::Array{T,1},k=1::Int) where {T}

        #ancillary variables
        γ=xgrid.grid
        a_f=MeshArray(γ,Int); [a_f[ii][:,:].=ii for ii=1:γ.nFaces]
        a_i=MeshArray(γ,Int); [a_i[ii]=collect(1:γ.fSize[ii][1])*ones(Int,1,γ.fSize[ii][2]) for ii=1:γ.nFaces]
        a_j=MeshArray(γ,Int); [a_j[ii]=ones(Int,γ.fSize[ii][1],1)*collect(1:γ.fSize[ii][2])' for ii=1:γ.nFaces]

        #convert to flat Array format
        a_x=write(xgrid)
        a_y=write(ygrid)
        a_f=write(a_f)
        a_i=write(a_i)
        a_j=write(a_j)

        #vector of grid points in Cartesian, 3D, coordinates
        kk=findall(isfinite.(a_x))
        x=sin.(pi/2 .-a_y[kk]*pi/180).*cos.(a_x[kk]*pi/180)
        y=sin.(pi/2 .-a_y[kk]*pi/180).*sin.(a_x[kk]*pi/180)
        z=cos.(pi/2 .-a_y[kk]*pi/180);

        #vector of target points in Cartesian, 3D, coordinates
        xx=sin.(pi/2 .-yvec*pi/180).*cos.(xvec*pi/180);
        yy=sin.(pi/2 .-yvec*pi/180).*sin.(xvec*pi/180);
        zz=cos.(pi/2 .-yvec*pi/180);

        #define tree
        kdtree = KDTree([x y z]')

        #find nearest neighbors
        idxs, _ = knn(kdtree, [xx yy zz]', k, true)
        idxs=[idxs[i][j] for i=1:length(idxs),j=1:k]

        return a_f[kk[idxs]],a_i[kk[idxs]],a_j[kk[idxs]],kk[idxs]
end

knn(xgrid::MeshArray,ygrid::MeshArray,lon::Number,lat::Number) = knn(xgrid::MeshArray,ygrid::MeshArray,[lon],[lat])

"""
    Interpolate(z_in::MeshArray,f,i,j,w)

```
using MeshArrays
lon=[i for i=-179.5:1.0:179.5, j=-89.5:1.0:89.5]
lat=[j for i=-179.5:1.0:179.5, j=-89.5:1.0:89.5]

Γ=GridLoad(GridSpec("LatLonCap",MeshArrays.GRID_LLC90))
(f,i,j,w,j_f,j_x,j_y)=InterpolationFactors(Γ,vec(lon),vec(lat))
DD=Interpolate(Γ["Depth"],f,i,j,w)

using Plots
contourf(vec(lon[:,1]),vec(lat[1,:]),DD,clims=(0.,6000.))
```
"""
function Interpolate(z_in::MeshArray,f,i,j,w)
    z_out=NaN*similar(f[:,1])
    for jj=1:size(f,1)
        if !isnan(sum(w[jj,:]))
            x=[z_in[f[jj,ii]][i[jj,ii],j[jj,ii]] for ii=1:4]
            kk=findall(isfinite.(x))
            ~isempty(kk) ? z_out[jj]=sum(w[jj,kk].*x[kk])/sum(w[jj,kk]) : nothing
        end
    end
    return z_out
end

"""
    InterpolationFactors(Γ,lon::Array{T,1},lat::Array{T,1})

Compute interpolation coefficients etc from grid `Γ` to `lon,lat`

```jldoctest
using MeshArrays
γ=GridSpec("CubeSphere",MeshArrays.GRID_CS32)
Γ=GridLoad(γ)
lon=collect(45.:0.1:46.); lat=collect(60.:0.1:61.)
(f,i,j,w,j_f,j_x,j_y)=InterpolationFactors(Γ,lon,lat)
YC=Interpolate(Γ["YC"],f,i,j,w)
extrema(i)==(9,10)

# output

true
```
"""
function InterpolationFactors(Γ,lon::Array{T,1},lat::Array{T,1}) where {T}
#to match gcmfaces test case (`interp=gcmfaces_interp_coeffs(0.1,0.1);`)
#set: iiTile=17; XC0=6.5000; YC0=-0.1994

        #main loop
        i_f=fill(0,length(lon),4)
        i_i=fill(0,length(lon),4)
        i_j=fill(0,length(lon),4)
        i_w=fill(NaN,length(lon),4)
        j_f=fill(0,length(lon),1)
        j_x=fill(0.0,length(lon),1)
        j_y=fill(0.0,length(lon),1)

        #ancillary variables
        (f,i,j,c)=knn(Γ["XC"],Γ["YC"],lon,lat)

        #1. tiles and τ        
        fs=Γ["XC"].fSize
        s=fill(0,2*length(fs))
        [s[collect(1:2) .+ (i-1)*2]=collect(fs[i]) for i in 1:length(fs)]
        ni=gcd(s); nj=gcd(s); γ=Γ["XC"].grid
        τ=Tiles(γ,ni,nj); tiles=MeshArray(γ,Int);
        [tiles[τ[ii]["face"]][τ[ii]["i"],τ[ii]["j"]].=ii for ii in 1:length(τ)]

        #2. t_XC, t_XC, t_f, t_i, t_j        
        t=vec(write(tiles)[c])
        t_list=unique(t)
        t_XC=Tiles(τ,exchange(Γ["XC"]))
        t_YC=Tiles(τ,exchange(Γ["YC"]))

        t_f=MeshArray(γ,Int); [t_f[ii][:,:].=ii for ii=1:γ.nFaces]
        t_i=MeshArray(γ,Int); [t_i[ii]=collect(1:γ.fSize[ii][1])*ones(Int,1,γ.fSize[ii][2]) for ii=1:γ.nFaces]
        t_j=MeshArray(γ,Int); [t_j[ii]=ones(Int,γ.fSize[ii][1],1)*collect(1:γ.fSize[ii][2])' for ii=1:γ.nFaces]

        t_f=Tiles(τ,exchange(t_f))
        t_i=Tiles(τ,exchange(t_i))
        t_j=Tiles(τ,exchange(t_j))

        x_q=fill(0.0,1,4)
        y_q=fill(0.0,1,4)
        tmpx=fill(0.0,1,4)
        tmpy=fill(0.0,1,4)
        w=fill(0.0,1,4)

        #main loop
        for ii=1:length(t_list)
                tt=t_list[ii]

                ff=τ[tt]["face"]
                ii0=minimum(τ[tt]["i"])+Int(ni/2)
                jj0=minimum(τ[tt]["j"])+Int(nj/2)
                XC0=Γ["XG"].f[ff][ii0,jj0]
                YC0=Γ["YG"].f[ff][ii0,jj0]
                #
                (x_grid,y_grid)=StereographicProjection(XC0,YC0,t_XC[tt],t_YC[tt])
                (x_quad,y_quad,i_quad,j_quad)=QuadArrays(x_grid,y_grid)
                #
                x=minimum(τ[tt]["i"])-0.5 .+collect(-1:ni)*ones(Int,1,nj+2)
                y=minimum(τ[tt]["j"])-0.5 .+ones(Int,ni+2,1)*collect(-1:nj)'
                #
                angsum=fill(0.0,size(x_quad,1))

                for pp in findall(t.==tt)
                        (x_trgt,y_trgt)=StereographicProjection(XC0,YC0,lon[pp],lat[pp])
                        PolygonAngle(x_quad,y_quad,x_trgt,y_trgt,angsum)
                        #angsum=PolygonAngle_deprecated(x_quad,y_quad,x_trgt,y_trgt)
                        if sum(angsum.>180.)>0
                                kk=findall(angsum.>180.)[end]
                                i_f[pp,:].=[t_f[tt][i_quad[kk,i]+1,j_quad[kk,i]+1] for i=1:4]
                                i_i[pp,:].=[t_i[tt][i_quad[kk,i]+1,j_quad[kk,i]+1] for i=1:4]
                                i_j[pp,:].=[t_j[tt][i_quad[kk,i]+1,j_quad[kk,i]+1] for i=1:4]
                                x_q[:].=x_quad[kk,:]
                                y_q[:].=y_quad[kk,:]
                                QuadCoeffs(x_q,y_q,x_trgt,y_trgt,w)
                                i_w[pp,:].=w[:]
                                #
                                [tmpx[i]=x[i_quad[kk,i]+1,j_quad[kk,i]+1] for i=1:4]
                                [tmpy[i]=y[i_quad[kk,i]+1,j_quad[kk,i]+1] for i=1:4]
                                #
                                j_f[pp]=ff
                                j_x[pp]=sum(tmpx[:].*i_w[pp,:])
                                j_y[pp]=sum(tmpy[:].*i_w[pp,:])
                        end
                end
        end

        return i_f,i_i,i_j,i_w,j_f,j_x,j_y
end

InterpolationFactors(Γ,lon::Number,lat::Number) = InterpolationFactors(Γ,[lon],[lat])

"""
    StereographicProjection(XC0,YC0,XC,YC)

Apply stereographic projection that puts `XC0,YC0` at `0.0,0.0`
to target point(s) `XC,YC`

```
lon=collect(45.:0.1:46.); lat=collect(60.:0.1:61.)
x,y=StereographicProjection(45.,60.,lon,lat)
```
"""
function StereographicProjection(XC0::Number,YC0::Number,XC,YC)
        #compute spherical coordinates:
        phi=XC; theta=90 .-YC;
        phi0=XC0; theta0=90-YC0;

        #compute cartesian coordinates:
        X=sind.(theta).*cosd.(phi);
        Y=sind.(theta).*sind.(phi);
        Z=cosd.(theta);

        x=X; y=Y; z=Z;

        #bring chosen point to the north pole:
        xx=x; yy=y; zz=z;
        x=cosd(phi0).*xx+sind(phi0).*yy;
        y=-sind(phi0).*xx+cosd(phi0).*yy;
        z=zz;

        xx=x; yy=y; zz=z;
        x=cosd(theta0)*xx-sind(theta0)*zz;
        y=yy;
        z=sind(theta0)*xx+cosd(theta0)*zz;

        #stereographic projection from the south pole:
        xx=x./(1 .+z);
        yy=y./(1 .+z);

        #nrm=sqrt(xx.^2+yy.^2);
        #msk=1+0*nrm; msk(nrm>tan(pi/4/2))=NaN;%mask points outside of pi/4 cone

        return xx,yy
end


"""
    PolygonAngle(px,py,x=[],y=[])

Compute sum of interior angles for polygons or points-to-polygons (when
`px,py,x,y` is provided as input). `px,py` are `MxN` matrices where each line
specifies one polygon. (optional) `x,y` are position vectors.

```jldoctest
using MeshArrays
px=[0. 0. 1. 1.]; py=[0. 1. 1. 0.];
x=collect(-1.0:0.25:2.0); y=x;
tmp=MeshArrays.PolygonAngle_deprecated(px,py,x,y)

isa(tmp,Array)

# output

true
```
"""
function PolygonAngle_deprecated(px,py,x=[],y=[])

        M=size(px,1); N=size(px,2); P=1;
        doPointsInPolygon=false
        if length(x)>0;
                doPointsInPolygon=true
                x=reshape(x,1,length(x))
                y=reshape(y,1,length(y))
                P=length(x)
        end;

        angsum=zeros(M,P)
        for ii=0:N-1
                ppx=circshift(px,[0 -ii])
                ppy=circshift(py,[0 -ii])
                if ~doPointsInPolygon
                        #compute sum of corner angles
                        v1x=ppx[:,2]-ppx[:,1]
                        v1y=ppy[:,2]-ppy[:,1]
                        v2x=ppx[:,4]-ppx[:,1]
                        v2y=ppy[:,4]-ppy[:,1]
                else;
                        #compute sum of sector angles
                        v1x=ppx[:,1]*ones(1,P)-ones(M,1)*x
                        v1y=ppy[:,1]*ones(1,P)-ones(M,1)*y
                        v2x=ppx[:,2]*ones(1,P)-ones(M,1)*x
                        v2y=ppy[:,2]*ones(1,P)-ones(M,1)*y
                end
                tmp=( v1x.*v2x+v1y.*v2y )./sqrt.( v1x.*v1x+v1y.*v1y )./sqrt.( v2x.*v2x+v2y.*v2y )
                g_acos=acos.( min.(max.(tmp,-1.0),1.0) )
                g_sin= ( v1x.*v2y-v1y.*v2x )./sqrt.( v1x.*v1x+v1y.*v1y )./sqrt.( v2x.*v2x+v2y.*v2y )
                angsum=angsum+rad2deg.(g_acos.*sign.(g_sin));
        end;

        return angsum
end

function PolygonAngle(px::Array,py::Array,angsum::Array)
        M=size(px,1)
        N=size(px,2)
        angsum .= 0.0
        for ii=0:N-1
                i1=mod1(ii+1,N)
                i2=mod1(ii+2,N)
                i4=mod1(ii+4,N)
                for jj=1:M
                #compute sum of sector angles
                v1x=px[jj,i2]-px[jj,i1]
                v1y=py[jj,i2]-py[jj,i1]
                v2x=px[jj,i4]-px[jj,i1]
                v2y=py[jj,i4]-py[jj,i1]

                tmp=( v1x.*v2x+v1y.*v2y )./sqrt.( v1x.*v1x+v1y.*v1y )./sqrt.( v2x.*v2x+v2y.*v2y )
                g_acos=acos.( min.(max.(tmp,-1.0),1.0) )
                g_sin= ( v1x.*v2y-v1y.*v2x )./sqrt.( v1x.*v1x+v1y.*v1y )./sqrt.( v2x.*v2x+v2y.*v2y )
                angsum[jj] += rad2deg(g_acos*sign(g_sin))
                end
        end
end

function PolygonAngle(px::Array,py::Array,x::Array,y::Array,angsum)
        for ii in 1:length(x)
                PolygonAngle(px,py,x[ii],y[ii],view(angsum,:,ii))
        end        
end

function PolygonAngle(px::Array,py::Array,x::Number,y::Number,angsum)
        M=size(px,1)
        N=size(px,2)
        angsum .= 0.0
        for ii=0:N-1
                i1=mod1(ii+1,N)
                i2=mod1(ii+2,N)
                for jj=1:M
                #compute sum of sector angles
                v1x=px[jj,i1] -x
                v1y=py[jj,i1] -y
                v2x=px[jj,i2] -x
                v2y=py[jj,i2] -y

                tmp=( v1x.*v2x+v1y.*v2y )./sqrt.( v1x.*v1x+v1y.*v1y )./sqrt.( v2x.*v2x+v2y.*v2y )
                g_acos=acos.( min.(max.(tmp,-1.0),1.0) )
                g_sin= ( v1x.*v2y-v1y.*v2x )./sqrt.( v1x.*v1x+v1y.*v1y )./sqrt.( v2x.*v2x+v2y.*v2y )
                angsum[jj] += rad2deg(g_acos*sign(g_sin))
                end
        end
end


"""
    QuadArrays(x_grid,y_grid)

Transform x_grid,y_grid (size ni+2,nj+2) into x_quad,y_quad,i_quad,j_quad
quadrilaterals (size ni+1*nj+1,4) where i_quad,j_quad are point indices
"""
function QuadArrays(x_grid::Array{T,2},y_grid::Array{T,2}) where {T}
        ni,nj=size(x_grid) .-2

        x_quad=Array{Float64,2}(undef,(ni+1)*(nj+1),4)
        y_quad=Array{Float64,2}(undef,(ni+1)*(nj+1),4)
        i_quad=Array{Int64,2}(undef,(ni+1)*(nj+1),4)
        j_quad=Array{Int64,2}(undef,(ni+1)*(nj+1),4)

        didj=[[0 0];[1 0];[1 1];[0 1]]
        for pp=1:4
                di=didj[pp,1]
                dj=didj[pp,2]

                #note the shift in indices due to exchange above
                tmp=x_grid[1+di:ni+1+di,1+dj:nj+1+dj]
                x_quad[:,pp]=vec(tmp)
                tmp=y_grid[1+di:ni+1+di,1+dj:nj+1+dj]
                y_quad[:,pp]=vec(tmp)

                tmp=collect(0+di:ni+di)*ones(1,nj+1)
                i_quad[:,pp]=vec(tmp)
                tmp=ones(ni+1,1)*transpose(collect(0+dj:nj+dj));
                j_quad[:,pp]=vec(tmp)
        end

        return x_quad,y_quad,i_quad,j_quad
end

"""
    QuadCoeffs(px,py,ox=[],oy=[],ow=[])

Compute bilinear interpolation coefficients for `ox,oy` within `px,py`
by remapping these quadrilaterals to the `unit square`.
- `px,py` are `Mx4` matrices where each line specifies one quadrilateral.
- `ox,oy` are `MxP` position matrices
- `ow` (output) are the `MxPx4` bilinear interpolation weights

```
QuadCoeffs([-1., 8., 13., -4.]',[-1., 3., 11., 8.]',0.,6.)
QuadCoeffs([0., 2., 3., 1.]',[0., 0., 1., 1.]',0.1,0.1)
```
"""
function QuadCoeffs(px,py,ox=[],oy=[],ow=[])
        #test case from https://www.particleincell.com/2012/quad-interpolation/
        #  QuadCoeffs([-1., 8., 13., -4.]',[-1., 3., 11., 8.]',0.,6.)
        #However the case of a perfect parallelogram needs special treatment
        #  QuadCoeffs([0., 2., 3., 1.]',[0., 0., 1., 1.]',0.1,0.1)
        #Deals with this situtation by falling back to ParaCoeffs
        #  ParaCoeffs([0., 2., 3., 1.]',[0., 0., 1., 1.]',0.1,0.1)

        #1. solve linear problem (`a,b` vectors from `px,py`)
        #  A=[1 0 0 0;1 1 0 0;1 1 1 1;1 0 1 0]; AI = inv(A);
        #  AI=[1 0 0 0;-1 1 0 0;-1 0 0 1; 1 -1 1 -1];
        #  a = AI*px';
        #  b = AI*py';
        #This defines the mapping from logical `l,m` to physical `x,y` as
        #  x=a(1)+a(2)*l+a(3)*m+a(2)*l*m;
        #  y=b(1)+b(2)*l+b(3)*m+b(2)*l*m;

        a=[px[1] -px[1]+px[2] -px[1]+px[4] px[1]-px[2]+px[3]-px[4]]
        a[findall(abs.(a).<1e-8)].=0.0

        b=[py[1] -py[1]+py[2] -py[1]+py[4] py[1]-py[2]+py[3]-py[4]]
        b[findall(abs.(b).<1e-8)].=0.0

        #2. select between the two solutions (to 2nd order
        #non-linear problem below) using polygon interior angles
        angsum=fill(0.0,1)
        PolygonAngle(px,py,angsum)
        sgn=fill(NaN,1)
        isapprox(angsum[1],360.0,rtol=0.01) ? sgn[1]=1.0 : nothing
        isapprox(angsum[1],-360.0,rtol=0.01) ? sgn[1]=-1.0 : nothing
        #ii=findall(isnan.(angsum))
        #length(ii)>0 ? println("warning: edge point was found") : nothing

        #3. solve non-linear problem for `pl,pm` from `px,py` & `a,b`
        # This defines the mapping from physical `x,y` to logical `l,m`

        a=reshape(a,(size(a,1),1,size(a,2))); 
        b=reshape(b,(size(b,1),1,size(b,2))); 

        # quadratic equation coeffs, `aa*mm^2+bb*m+cc=0`
        if ~isempty(ox);
                x=ox; y=oy;
        else;
                x=px; y=py;
                a=repeat(a,1,size(x,2),1);
                b=repeat(b,1,size(x,2),1);
                sgn=repeat(sgn,1,size(x,2));
        end;
        #
        det=fill(0.0,size(x))
        pm=fill(0.0,size(x))
        pl=fill(0.0,size(x))
        for ii=1:size(x,1), jj=1:size(x,2)
                aa = a[ii,jj,4]*b[ii,jj,3]-a[ii,jj,3]*b[ii,jj,4]
                bb = a[ii,jj,4]*b[ii,jj,1]-a[ii,jj,1]*b[ii,jj,4]+a[ii,jj,2]*b[ii,jj,3]-a[ii,jj,3]*b[ii,jj,2]+x[ii,jj]*b[ii,jj,4]-y[ii,jj].*a[ii,jj,4]
                cc = a[ii,jj,2]*b[ii,jj,1]-a[ii,jj,1]*b[ii,jj,2]+x[ii,jj]*b[ii,jj,2]-y[ii,jj].*a[ii,jj,2]
                #compute `pm = (-b+sqrt(b^2-4ac))/(2a)`
                det[ii,jj] = sqrt(bb*bb - 4.0*aa*cc)
                pm[ii,jj]  = (-bb+sgn[ii,jj]*det[ii,jj])/(2.0*aa)
                #compute `pl` by substitution in equation system
                pl[ii,jj]  = (x[ii,jj]-a[ii,jj,1]-a[ii,jj,3]*pm[ii,jj])/(a[ii,jj,2]+a[ii,jj,4]*pm[ii,jj])
        end

        if ~isempty(ox);
                tmp1=(1 .-pl).*(1 .-pm)
                tmp2=pl.*(1 .-pm)
                tmp3=pl.*pm
                tmp4=(1 .-pl).*pm
                ow[:]=cat(tmp1,tmp2,tmp3,tmp4; dims=3)
                #ow[:].=cat(tmp1,tmp2,tmp3,tmp4; dims=3)[:]
                #treat pathological cases if needed
                tmp=ParaCoeffs(px,py,[ox],[oy])
                for kk=1:length(ow)
                        !isfinite(ow[kk]) ? ow[kk]=tmp[kk] : nothing
                end
        end

end

"""
    ParaCoeffs(px,py,ox=[],oy=[])

Compute bilinear interpolation coefficients for `ox,oy` within `px,py`
by remapping these parallelograms to the `unit square`.
- `px,py` are `Mx4` matrices where each line specifies one quadrilateral.
- `ox,oy` are `MxP` position matrices
- `pw` (output) are the `MxPx4` bilinear interpolation weights

```
x=1.0; y=1.0 #Try send the corners to unit square corners?
println(vec(ParaCoeffs([0., 2., 3., 1.]',[0., 0., 1., 1.]',x,y)))
println(vec(QuadCoeffs([0., 2.01, 3., 1.]',[0., 0., 1., 1.]',x,y)))
```
"""
function ParaCoeffs(px,py,ox=[],oy=[])

        tmp1=px[:,1];
        tmp2=-px[:,1]+px[:,2];
        tmp3=-px[:,2]+px[:,3];
        a=[tmp1 tmp2 tmp3];

        tmp1=py[:,1];
        tmp2=-py[:,1]+py[:,2];
        tmp3=-py[:,2]+py[:,3];
        b=[tmp1 tmp2 tmp3];

#        (m,l)=inv([a[1,2] a[1,3];b[1,2] b[1,3]])*[ox[1]-a[1,1]; oy[1]-b[1,1]]
        m=( b[:,3].*(ox-a[:,1])-a[:,3].*(oy-b[:,1]) ) ./(a[:,2].*b[:,3]-a[:,3].*b[:,2])
        l=( -b[:,2].*(ox-a[:,1])+a[:,2].*(oy-b[:,1]) ) ./(a[:,2].*b[:,3]-a[:,3].*b[:,2])

        ow=[];
        tmp1=(1 .-l).*(1 .-m)
        tmp4=l.*(1 .-m)
        tmp3=l.*m
        tmp2=(1 .-l).*m
        ow=cat(tmp1,tmp2,tmp3,tmp4; dims=3)

        return ow
end

export AbstractMeshArray, MeshArray, InnerArray, OuterArray
export gcmgrid, exchange, gradient, convergence, smooth, mask
export simple_periodic_domain, GridSpec, GridLoad, GridOfOnes, GridAddWS!
export Tiles, Interpolate, InterpolationFactors, knn
export ScalarPotential, VectorPotential, ThroughFlow
export StereographicProjection, LatitudeCircles
export maximum, minimum, mean, std

#The following exch_UV differs from normal exchange; incl. exch_UV_N
export exch_UV

#ggebbie: think about using AbstractMeshArray as input type.
"""
    function maximum
    Author ggebbie
    Compute maximum value of gcmgrid type 
# Input
- `x::MeshArrays.gcmarray{Float32,1,Array{Float32,2}}`: input of gcmarray type
- `dryval::Float32`: land value to be eliminated in calculation
# Output
- `xmax::`: maximum value of 2D field
"""
function maximum(x::gcmarray{Float32,1,Array{Float32,2}},dryval::Float32)

    isdry(z) = (z == dryval)
    
    #  vector list of nonzero elements
    xcount = [sum(count(!isdry,x[i])) for i in eachindex(x)]
    if sum(xcount) > 0
        xmax = maximum([maximum(filter(!isdry,x[i])) for i in eachindex(x) if xcount[i] > 0])
    else
        xmax = NaN32
    end
    return xmax
end

"""
    function minimum
    Author ggebbie
    Compute minimum value of gcmgrid type 
# Input
- `x::MeshArrays.gcmarray{Float32,1,Array{Float32,2}}`: input of gcmarray type
- `dryval::Float32`: land value to be eliminated in calculation
# Output
- `xmin::Float32`: minimum value of 2D field
"""
function minimum(x::gcmarray{Float32,1,Array{Float32,2}},dryval::Float32)
    xmin = -maximum(-x,dryval)
    return xmin
end

"""
    function mean
    Author: ggebbie
    Mean of gcmgrid type 
# Input
- `x::MeshArrays.gcmarray{Float32,1,Array{Float32,2}}`: input of gcmarray type
- `dryval::Float32`: land value (doesn't work for NaN32)
# Output
- `xbar::Float32`: mean value (unweighted)
"""
function mean(x::gcmarray{Float32,1,Array{Float32,2}},dryval::Float32)

    isdry(z) = (z == dryval)
    
    #  vector list of nonzero elements
    xcount = [sum(count(!isdry,x[i])) for i in eachindex(x)]
    if sum(xcount) > 0
        # don't assume 0 on land
        xsum = sum([sum(filter(!isdry,x[i])) for i in eachindex(x) if xcount[i] > 0])
        xbar = xsum/sum(xcount)
    else
        xbar = NaN32
    end
    return xbar
end

"""
    function mean
    Author: ggebbie
    Area-weighted mean of gcmgrid type filtered by dryval
# Input
- `x::MeshArrays.gcmarray{Float32,1,Array{Float32,2}}`: input of gcmarray type
- `weight::MeshArrays.gcmarray{Float32,1,Array{Float32,2}}`: weighting variable of gcmarray type
- `dryval::Float32`: land value (doesn't work for NaN32)
# Output
- `xbar::Float32`: mean value (unweighted)
"""
function mean(x::gcmarray{Float32,1,Array{Float32,2}},weight::gcmarray{Float64,1,Array{Float64,2}},dryval::Float64)

    isdry(z) = (z == dryval)
    
    #  vector list of nonzero elements
    xcount = [sum(count(!isdry,x[i])) for i in eachindex(x)]
    if sum(xcount) > 0
        # don't assume 0 on land
        xsum = sum([sum(filter(!isdry,x[i].*weight[i])) for i in eachindex(x) if xcount[i] > 0])
        xdenom = sum([sum(filter(!isdry,weight[i])) for i in eachindex(weight) if xcount[i] > 0])
        xbar = xsum/xdenom
    else
        xbar = NaN32
    end
    return xbar
end

"""
    function mean
    Author: ggebbie
    Area-weighted mean of gcmgrid type filtered by a function
# Input
- `x::MeshArrays.gcmarray{Float32,1,Array{Float32,2}}`: input of gcmarray type
- `weight::MeshArrays.gcmarray{Float32,1,Array{Float32,2}}`: weighting variable of gcmarray type
- `isgood::Function`: returns true is a value to be used in the mean
# Output
- `xbar::Float32`: mean value (weighted and filtered)
"""
function mean(x::gcmarray{Float32,1,Array{Float32,2}},weight::gcmarray{Float64,1,Array{Float64,2}},isgood)

    #  vector list of nonzero elements
    xcount = [sum(count(isgood,x[i])) for i in eachindex(x)]
    if sum(xcount) > 0
        # don't assume 0 on land
        xsum = sum([sum(filter(isgood,x[i].*weight[i])) for i in eachindex(x) if xcount[i] > 0])
        xdenom = sum([sum(filter(isgood,weight[i])) for i in eachindex(weight) if xcount[i] > 0])
        xbar = xsum/xdenom
    else
        xbar = NaN32
    end
    return xbar
end

"""
    function std
    Author: ggebbie
    Compute standard deviation of gcmgrid type 
# Input
- `x::MeshArrays.gcmarray{Float32,1,Array{Float32,2}}`: input of gcmarray type
- `xbar::Float32`: mean value
- `dryval::Float32`: land value (doesn't work for NaN32)
# Output
- `σx::Float32`: standard deviation 
"""
function std(x::gcmarray{Float32,1,Array{Float32,2}},xbar::Float32,dryval)

    isdry(z) = (z == dryval)

    # x prime = fluctuation
    x′ = x .- xbar
    
    #  vector list of nonzero elements
    xcount = [sum(count(!isdry,x′[i])) for i in eachindex(x′)]
    if sum(xcount) > 0
        # don't assume 0 on land
        x²sum = sum([sum(filter(!isdry,x′[i]).^2) for i in eachindex(x) if xcount[i] > 0])
        σx = sqrt(x²sum/(sum(xcount)-1))
    else
        xbar = NaN32
    end
    return σx
end

end # module
