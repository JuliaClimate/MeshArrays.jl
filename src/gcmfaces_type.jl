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
#  grid=gcmgrid("", "llc", 5, fSize, [90 1170], T, read, write)
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
    similar(A)
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

# basic operations

import Base: +, -, *, /, getindex
import Base: isnan, isinf, isfinite
import Base: maximum, minimum, sum, fill

function +(a::gcmfaces)
  c=similar(a)
  for iFace=1:a.grid.nFaces
    c.f[iFace]=+a.f[iFace];
  end
  return c
end

function +(a::gcmfaces,b::gcmfaces)
  c=similar(a)
  for iFace=1:a.grid.nFaces
    c.f[iFace]=a.f[iFace]+b.f[iFace];
  end
  return c
  #the following modifies a:
  #  c=a
  #  c.f=a.f .+ b.f
  #the following fails immutability:
  #  c=gcmfaces()
  #  c.f=a.f .+ b.f
end

function +(a::Number,b::gcmfaces)
  c=similar(b)
  for iFace=1:b.grid.nFaces
    c.f[iFace]=a.+b.f[iFace];
  end
  return c
  #the following is deprecated synthax as of v0.7:
  #  c=gcmfaces()
  #  c.f=a .+ b.f
end

function +(a::gcmfaces,b::Number)
  c=b+a
  return c
end

function -(a::gcmfaces)
  c=similar(a)
  for iFace=1:a.grid.nFaces
    c.f[iFace]=-a.f[iFace];
  end
  return c
end

function -(a::gcmfaces,b::gcmfaces)
  c=similar(a)
  for iFace=1:a.grid.nFaces
    c.f[iFace]=a.f[iFace]-b.f[iFace];
  end
  return c
end

function -(a::Number,b::gcmfaces)
  c=similar(b)
  for iFace=1:b.grid.nFaces
    c.f[iFace]=a.-b.f[iFace];
  end
  return c
end

function -(a::gcmfaces,b::Number)
  c=similar(a)
  for iFace=1:a.grid.nFaces
    c.f[iFace]=a.f[iFace].-b;
  end
  return c
end

function *(a::gcmfaces,b::gcmfaces)
  c=similar(a);
  for iFace=1:a.grid.nFaces
    c.f[iFace]=a.f[iFace].*b.f[iFace];
  end
  return c
end

function *(a::Number,b::gcmfaces)
  c=similar(b);
  for iFace=1:b.grid.nFaces
    c.f[iFace]=a*b.f[iFace];
  end
  return c
end

function *(a::gcmfaces,b::Number)
  c=b*a
  return c
end

function /(a::gcmfaces,b::gcmfaces)
  c=similar(a);
  for iFace=1:a.grid.nFaces
    c.f[iFace]=a.f[iFace]./b.f[iFace];
  end
  return c
end

function /(a::Number,b::gcmfaces)
  c=similar(b);
  for iFace=1:b.grid.nFaces
    c.f[iFace]=a./b.f[iFace];
  end
  return c
end

function /(a::gcmfaces,b::Number)
  c=similar(a);
  for iFace=1:a.grid.nFaces
    c.f[iFace]=a.f[iFace]./b;
  end
  return c
end

##

function isnan(a::gcmfaces)
    c=similar(a);
    for iFace=1:a.grid.nFaces
      c.f[iFace]=isnan.(a.f[iFace]);
    end
    return c
end

function isinf(a::gcmfaces)
    c=similar(a);
    for iFace=1:a.grid.nFaces
      c.f[iFace]=isinf.(a.f[iFace]);
    end
    return c
end

function isfinite(a::gcmfaces)
    c=similar(a);
    for iFace=1:a.grid.nFaces
      c.f[iFace]=isfinite.(a.f[iFace]);
    end
    return c
end

function sum(a::gcmfaces)
    c=0.0;
    for iFace=1:a.grid.nFaces
      tmp1=a.f[iFace];
      c=c+sum(tmp1);
    end
    return c
end

function maximum(a::gcmfaces)
    c=-Inf;
    for iFace=1:a.grid.nFaces
      tmp1=a.f[iFace];
      c=max(c,maximum(tmp1));
    end
    return c
end

function minimum(a::gcmfaces)
    c=Inf;
    for iFace=1:a.grid.nFaces
      tmp1=a.f[iFace];
      c=min(c,minimum(tmp1));
    end
    return c
end

function fill(val::Any,a::gcmfaces)
    c=similar(a);
    for iFace=1:a.grid.nFaces
      c.f[iFace]=fill(val,fsize(a,iFace));
    end
    return c
end

###

function nFacesEtc(a::gcmfaces)
  nFaces=length(a.f)
  ndims(a.f[1])>2 ? n3=size(a.f[1],3) : n3=1
  return nFaces, n3
end
