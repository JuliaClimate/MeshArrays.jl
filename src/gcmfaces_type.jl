## gcmfaces type definition + methods

## type definition

abstract type AbstractGcmfaces{T, N} <: AbstractArray{T, N} end

"""
    gcmfaces{T, N}

gcmfaces data structure. Available constructors:

```
gcmfaces{T,N}(nFaces::Int,grTopo::String,f::Array{Array{T,N},1},
         fSize::Array{NTuple{N, Int}}, aSize::NTuple{N,Int})
gcmfaces(nFaces::Int,grTopo::String,::Type{T},
         fSize::Array{NTuple{N, Int}}, aSize::NTuple{N,Int}) where {T,N}
gcmfaces(nFaces::Int,grTopo::String,v1::Array{Array{T,N},1}) where {T,N}
gcmfaces(A::AbstractGcmfaces{T, N}) where {T,N}
gcmfaces()
```
"""
struct gcmfaces{T, N} <: AbstractGcmfaces{T, N}
   nFaces::Int
   grTopo::String
   f::Array{Array{T,N},1}
   fSize::Array{NTuple{N, Int}}
   aSize::NTuple{N, Int}
end

"""
    gcmsubset{T, N}

gcmsubset data structure. Available constructors:

```
gcmsubset{T,N}(nFaces::Int,grTopo::String,f::Array{Array{T,N},1},
               fSize::Array{NTuple{N, Int}},aSize::NTuple{N, Int},
               i::Array{Array{T,N},1},iSize::Array{NTuple{N, Int}})
gcmsubset(nFaces::Int,grTopo::String,::Type{T},fSize::Array{NTuple{N, Int}},
          aSize::NTuple{N,Int},dims::NTuple{N,Int}) where {T,N}
```
"""
struct gcmsubset{T, N} <: AbstractGcmfaces{T, N}
   nFaces::Int
   grTopo::String
   f::Array{Array{T,N},1}
   fSize::Array{NTuple{N, Int}}
   aSize::NTuple{N, Int}
   i::Array{Array{T,N},1}
   iSize::Array{NTuple{N, Int}}
end

## additional constructors for gcmfaces

function gcmfaces(nFaces::Int,grTopo::String,::Type{T},
  fSize::Array{NTuple{N, Int}},
  aSize::NTuple{N,Int}) where {T,N}
  f=Array{Array{T,N},1}(undef,nFaces)
  for a=1:nFaces
    f[a]=Array{T}(undef,fSize[a])
  end
  gcmfaces{T,N}(nFaces,grTopo,f,fSize,aSize)
end

function gcmfaces(nFaces::Int,grTopo::String,
  v1::Array{Array{T,N},1}) where {T,N}
  fSize=fsize(v1)
  aSize=fsize(v1,0)
  gcmfaces{T,N}(nFaces,grTopo,deepcopy(v1),fSize,aSize)
#  gcmfaces(nFaces,grTopo,T,fs,as)
end

function gcmfaces(A::AbstractGcmfaces{T, N}) where {T,N}
  #should this be called similar? deepcopy?
  fSize=fsize(A)
  aSize=size(A)
  gcmfaces{T,N}(nFaces,grTopo,deepcopy(A.f),fSize,aSize)
#  gcmfaces(nFaces,grTopo,T,fSize,aSize)
end

function gcmfaces()
  if isdefined(MeshArrays,:nFaces)
    nFaces=MeshArrays.nFaces
    grTopo=MeshArrays.grTopo
    T=MeshArrays.ioPrec
    fSize=MeshArrays.facesSize
    aSize=(prod(MeshArrays.ioSize),1)
  else
    nFaces=5
    grTopo="llc"
    T=Float64
    fSize=[(90, 270), (90, 270), (90, 90), (270, 90), (270, 90)]
    aSize=(105300, 1);
  end
  gcmfaces(nFaces,grTopo,T,fSize,aSize)
end

## additional constructors for gcmsubset

#maybe: replace this constructor with one that gets A and sets f to view(A.f)
function gcmsubset(nFaces::Int,grTopo::String,::Type{T},
  fSize::Array{NTuple{N, Int}},aSize::NTuple{N,Int},
  dims::NTuple{N,Int}) where {T,N}
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
  gcmsubset{T,N}(nFaces,grTopo,f,fSize,aSize,i,iSize)
end

## Convenience functions

"""
    fijind(A::gcmfaces,ij::Int)

Compute face and local indices (f,j,k) from global index (ij).
"""
function fijind(A::gcmfaces,ij::Int)
  f=0
  j=0
  k=0
  tmp1=0
  for iFace=1:A.nFaces
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
    fsize(A::AbstractGcmfaces{T, N}) where {T,N}

Return vector of face array sizes. Other methods:
```
fsize(A::AbstractGcmfaces{T, N},i::Int) where {T,N}
fsize(A::Array{Array{T,N},1}) where {T,N}
fsize(A::Array{Array{T,N},1},i::Int) where {T,N}
```
"""
function fsize(A::AbstractGcmfaces{T, N}) where {T,N}
  fs=Array{NTuple{N, Int}}(undef,A.nFaces)
  for i=1:A.nFaces
    fs[i]=size(A.f[i]);
  end
  return fs
end

function fsize(A::AbstractGcmfaces{T, N},i::Int) where {T,N}
  if i>0
    fs=size(A.f[i])
  else
    tmp1=0
    for i=1:A.nFaces
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

function Base.getindex(A::AbstractGcmfaces{T, N}, I::Vararg{Union{Int,AbstractUnitRange,Colon}, N}) where {T,N}
  if typeof(I[1])<:Int
    (f,i,j)=fijind(A,I[1])
    J=Base.tail(Base.tail(I))
    J=(i,j,J...)
    val=A.f[f][J...]
  elseif typeof(I[1])<:AbstractUnitRange
    val=similar(A,eltype(A),length.(I))
    for iFace=1:A.nFaces
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

function Base.setindex!(A::AbstractGcmfaces{T, N}, v, I::Vararg{Int, N}) where {T,N}
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

function Base.view(a::AbstractGcmfaces{T, N}, I::Vararg{Union{Int,AbstractUnitRange,Colon}, N}) where {T,N}
  nFaces=a.nFaces;
  grTopo=a.grTopo;
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
  c=gcmfaces(nFaces,grTopo,v1);
  return c;
end

# Custom pretty-printing

function Base.show(io::IO, z::AbstractGcmfaces{T, N}) where {T,N}

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
    printstyled(io, "$(z.grTopo)\n",color=:blue)
    printstyled(io, "  # of faces  = ",color=:normal)
    printstyled(io, "$(z.nFaces)\n",color=:blue)
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
      for iFace=2:z.nFaces
        printstyled(io, "                ",color=:normal)
        printstyled(io, "$(fs[iFace])\n",color=:blue)
      end
    end

    return
end

#

function Base.similar(A::gcmfaces, ::Type{T}, dims::Dims) where {T}
  if prod(dims)==length(A)
    B=gcmfaces(A.nFaces,A.grTopo,T,A.fSize,A.aSize)
  else
    B=gcmsubset(A.nFaces,A.grTopo,T,A.fSize,A.aSize,dims)
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
    B=gcmsubset(A.nFaces,A.grTopo,T,A.fSize,A.aSize,dims[1])
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
  for iFace=1:a.nFaces
    c.f[iFace]=+a.f[iFace];
  end
  return c
end

function +(a::gcmfaces,b::gcmfaces)
  c=similar(a)
  for iFace=1:a.nFaces
    c.f[iFace]=a.f[iFace]+b.f[iFace];
  end
  return c
  #the following modifies a:
  #  c=a
  #  c.f=a.f .+ b.f
  #the following fails immutability:
  #  c=gcmfaces(a.nFaces,a.grTopo)
  #  c.f=a.f .+ b.f
end

function +(a::Number,b::gcmfaces)
  c=similar(b)
  for iFace=1:b.nFaces
    c.f[iFace]=a.+b.f[iFace];
  end
  return c
  #the following is deprecated synthax as of v0.7:
  #  c=gcmfaces(b.nFaces,b.grTopo)
  #  c.f=a .+ b.f
end

function +(a::gcmfaces,b::Number)
  c=b+a
  return c
end

function -(a::gcmfaces)
  c=similar(a)
  for iFace=1:a.nFaces
    c.f[iFace]=-a.f[iFace];
  end
  return c
end

function -(a::gcmfaces,b::gcmfaces)
  c=similar(a)
  for iFace=1:a.nFaces
    c.f[iFace]=a.f[iFace]-b.f[iFace];
  end
  return c
end

function -(a::Number,b::gcmfaces)
  c=similar(b)
  for iFace=1:b.nFaces
    c.f[iFace]=a.-b.f[iFace];
  end
  return c
end

function -(a::gcmfaces,b::Number)
  c=similar(a)
  for iFace=1:a.nFaces
    c.f[iFace]=a.f[iFace].-b;
  end
  return c
end

function *(a::gcmfaces,b::gcmfaces)
  c=similar(a);
  for iFace=1:a.nFaces
    c.f[iFace]=a.f[iFace].*b.f[iFace];
  end
  return c
end

function *(a::Number,b::gcmfaces)
  c=similar(b);
  for iFace=1:b.nFaces
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
  for iFace=1:a.nFaces
    c.f[iFace]=a.f[iFace]./b.f[iFace];
  end
  return c
end

function /(a::Number,b::gcmfaces)
  c=similar(b);
  for iFace=1:b.nFaces
    c.f[iFace]=a./b.f[iFace];
  end
  return c
end

function /(a::gcmfaces,b::Number)
  c=similar(a);
  for iFace=1:a.nFaces
    c.f[iFace]=a.f[iFace]./b;
  end
  return c
end

##

function isnan(a::gcmfaces)
    c=similar(a);
    for iFace=1:a.nFaces
      c.f[iFace]=isnan.(a.f[iFace]);
    end
    return c
end

function isinf(a::gcmfaces)
    c=similar(a);
    for iFace=1:a.nFaces
      c.f[iFace]=isinf.(a.f[iFace]);
    end
    return c
end

function isfinite(a::gcmfaces)
    c=similar(a);
    for iFace=1:a.nFaces
      c.f[iFace]=isfinite.(a.f[iFace]);
    end
    return c
end

function sum(a::gcmfaces)
    c=0.0;
    for iFace=1:a.nFaces
      tmp1=a.f[iFace];
      c=c+sum(tmp1);
    end
    return c
end

function maximum(a::gcmfaces)
    c=-Inf;
    for iFace=1:a.nFaces
      tmp1=a.f[iFace];
      c=max(c,maximum(tmp1));
    end
    return c
end

function minimum(a::gcmfaces)
    c=Inf;
    for iFace=1:a.nFaces
      tmp1=a.f[iFace];
      c=min(c,minimum(tmp1));
    end
    return c
end

function fill(val::Any,a::gcmfaces)
    c=similar(a);
    for iFace=1:a.nFaces
      c.f[iFace]=fill(val,fsize(a,iFace));
    end
    return c
end
