
## gcmfaces type definition + methods

## type definition

abstract type AbstractGcmfaces{T, N} <: AbstractArray{T, N} end

struct gcmfaces{T, N} <: AbstractGcmfaces{T, N}
   nFaces::Int
   grTopo::String
   f::Array{Array{T,N},1}
end

# additional constructors

function gcmfaces(nFaces::Int,grTopo::String,v1::AbstractArray)
  tmp_eltype=eltype(v1[1])
  tmp_ndims=ndims(v1[1])
  gcmfaces{tmp_eltype,tmp_ndims}(nFaces,grTopo,v1)
end

function gcmfaces(nFaces::Int,grTopo::String)
  tmp0=Array{Float64}(undef,1,1)
  tmp1=fill(tmp0,nFaces)
  gcmfaces{Float64,2}(nFaces,grTopo,tmp1)
end

gcmfaces() = gcmfaces(MeshArrays.nFaces,MeshArrays.grTopo);

## Convenience functions

function fsize(A::AbstractGcmfaces{T, N}) where {T,N}
  fs=Array{NTuple{N, Int}}(undef,A.nFaces)
  for i=1:A.nFaces
    fs[i]=size(A.f[i]);
  end
  return fs
end

function fsize(A::AbstractGcmfaces{T, N},i::Integer) where {T,N}
  fs=size(A.f[i]);
end

## Interface Methods

function Base.size(A::AbstractGcmfaces{T, N}) where {T,N}
  fs=fsize(A)
  tmp1=0
  for i=1:A.nFaces
    tmp1=tmp1+fs[i][1]*fs[i][2]
  end
  #N>2 ? s=(tmp1,1,fs[1][3]) : s=(tmp1,1)
  if N==2; s=(tmp1,1);
  elseif N==3; s=(tmp1,1,fs[1][3]);
  elseif N==4; s=(tmp1,1,fs[1][3],fs[1][4]);
  else; error("N>4 case not implemented yet");
  end
  return s
end

function Base.size(A::AbstractGcmfaces{T, N},dim::Integer) where {T,N}
  tmp1=size(A)
  s=tmp1[dim]
end

function Base.getindex(A::AbstractGcmfaces{T, N}, i::Vararg{Int, N}) where {T,N}

  f=0
  j=0
  k=0
  tmp1=0
  for iFace=1:A.nFaces
    tmpsize=fsize(A,iFace)
    tmp11=tmpsize[1]*tmpsize[2]
    tmp2=tmp1+tmp11
    if tmp1<i[1]<=tmp2
      f=iFace;
      tmp3=(i[1]-tmp1);
      k=Int(ceil(tmp3/tmpsize[1]))
      j=Int(tmp3-tmpsize[1]*(k-1))
    end;
    tmp1=tmp1+tmp11
  end
  if ndims(A)==2; val=A.f[f][j,k];
  elseif ndims(A)==3; val=A.f[f][j,k,i[3]];
  elseif ndims(A)==4; val=A.f[f][j,k,i[3],i[4]];
  else; error("ndims>4 case not implemented yet")
  end
  return val
end

## view

function Base.view(a::AbstractGcmfaces{T, N}, I::Vararg{Any, N}) where {T,N}
  nFaces=a.nFaces;
  grTopo=a.grTopo;
  v1=Array{Any}(undef,nFaces);
  for iFace=1:nFaces
    if N==2;
      v1[iFace]=view(a.f[iFace],:,:);
    elseif N==3;
      v1[iFace]=view(a.f[iFace],:,:,I[3]);
    elseif N==4;
      v1[iFace]=view(a.f[iFace],:,:,I[3],I[4]);
    else;
      error("N>4 case not implemented yet")
    end
  end
  c=gcmfaces(nFaces,grTopo,v1);
  return c;
end

# Custom pretty-printing

function Base.show(io::IO, z::gcmfaces)

#    @printf io " MeshArrays instance with \n"
    printstyled(io, " gcmfaces array \n",color=:normal)
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
      printstyled(io, "  face sizes  = ",color=:normal)
      printstyled(io, "$(size(z.f[1]))\n",color=:blue)
      for iFace=2:z.nFaces
        printstyled(io, "                ",color=:normal)
        printstyled(io, "$(size(z.f[iFace]))\n",color=:blue)
      end
    end

    return
end

# basic operations

import Base: +, -, *, /, getindex
import Base: isnan, isinf, isfinite
import Base: maximum, minimum, sum, fill

function +(a::gcmfaces)
  c=gcmfaces(a.nFaces,a.grTopo,a.f)
  return c
end

function +(a::gcmfaces,b::gcmfaces)
  cf=a.f+b.f
  c=gcmfaces(a.nFaces,a.grTopo,cf)
  return c
  #the following modifies a:
  #  c=a
  #  c.f=a.f .+ b.f
  #the following fails immutability:
  #  c=gcmfaces(a.nFaces,a.grTopo)
  #  c.f=a.f .+ b.f
end

function +(a::Number,b::gcmfaces)
  nFaces=b.nFaces;
  grTopo=b.grTopo;
  v1=Array{Any}(undef,nFaces);
  for iFace=1:nFaces
    tmp1=b.f[iFace];
    tmp2=a*ones(Float64, size(tmp1));
    v1[iFace]=tmp1+tmp2;
  end
  c=gcmfaces(nFaces,grTopo,v1);
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
  c=gcmfaces(a.nFaces,a.grTopo,-a.f)
  return c
end

function -(a::gcmfaces,b::gcmfaces)
  cf=a.f .- b.f
  c=gcmfaces(a.nFaces,a.grTopo,cf)
  return c
end

function -(a::Number,b::gcmfaces)
  nFaces=b.nFaces;
  grTopo=b.grTopo;
  v1=Array{Any}(undef,nFaces);
  for iFace=1:nFaces
    tmpb=b.f[iFace];
    tmpa=a*ones(Float64,size(tmpb));
    v1[iFace]=tmpa-tmpb;
  end
  c=gcmfaces(b.nFaces,b.grTopo,v1);
  return c
end

function -(a::gcmfaces,b::Number)
  nFaces=a.nFaces;
  grTopo=a.grTopo;
  v1=Array{Any}(undef,nFaces);
  for iFace=1:nFaces
    tmpa=a.f[iFace];
    tmpb=b*ones(Float64,size(tmpa));
    v1[iFace]=tmpa-tmpb;
  end
  c=gcmfaces(a.nFaces,a.grTopo,v1);
  return c
end

function *(a::gcmfaces,b::gcmfaces)
  nFaces=a.nFaces;
  grTopo=a.grTopo;
  v1=Array{Any}(undef,nFaces);
  for iFace=1:nFaces
    v1[iFace]=a.f[iFace] .* b.f[iFace];
  end
  c=gcmfaces(nFaces,grTopo,v1);
  return c
end

function *(a::Number,b::gcmfaces)
  v1=a .* b.f;
  c=gcmfaces(b.nFaces,b.grTopo,v1);
  return c
end

function *(a::gcmfaces,b::Number)
  v1=b .* a.f;
  c=gcmfaces(a.nFaces,a.grTopo,v1);
  return c
end

function /(a::gcmfaces,b::gcmfaces)
  nFaces=a.nFaces;
  grTopo=a.grTopo;
  v1=Array{Any}(undef,nFaces);
  for iFace=1:nFaces
    v1[iFace]=a.f[iFace] ./ b.f[iFace];
  end
  c=gcmfaces(nFaces,grTopo,v1);
  return c
end

function /(a::Number,b::gcmfaces)
  nFaces=b.nFaces;
  grTopo=b.grTopo;
  v1=Array{Any}(undef,nFaces);
  for iFace=1:nFaces
    v1[iFace]=a ./ b.f[iFace];
  end
  c=gcmfaces(nFaces,grTopo,v1);
  return c
end

function /(a::gcmfaces,b::Number)
  nFaces=a.nFaces;
  grTopo=a.grTopo;
  v1=a.f ./ b;
  c=gcmfaces(nFaces,grTopo,v1);
  return c
end

##

function isnan(a::gcmfaces)
    nFaces=a.nFaces;
    grTopo=a.grTopo;
    v1=Array{Any}(undef,nFaces);
    for iFace=1:nFaces
      tmp1=a.f[iFace];
      v1[iFace]=isnan.(tmp1);
    end
    c=gcmfaces(nFaces,grTopo,v1);
    return c
end

function isinf(a::gcmfaces)
    nFaces=a.nFaces;
    grTopo=a.grTopo;
    v1=Array{Any}(undef,nFaces);
    for iFace=1:nFaces
      tmp1=a.f[iFace];
      v1[iFace]=isinf.(tmp1);
    end
    c=gcmfaces(nFaces,grTopo,v1);
    return c
end

function isfinite(a::gcmfaces)
    nFaces=a.nFaces;
    grTopo=a.grTopo;
    v1=Array{Any}(undef,nFaces);
    for iFace=1:nFaces
      tmp1=a.f[iFace];
      v1[iFace]=isfinite.(tmp1);
    end
    c=gcmfaces(nFaces,grTopo,v1);
    return c
end

function sum(a::gcmfaces)
    c=0.;
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
    nFaces=a.nFaces;
    grTopo=a.grTopo;
    v1=Array{Any}(undef,nFaces);
    for iFace=1:nFaces
      tmp1=a.f[iFace];
      v1[iFace]=fill(val,size(tmp1));
    end
    c=gcmfaces(nFaces,grTopo,v1);
    return c
end
