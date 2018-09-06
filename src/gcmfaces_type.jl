
## MeshArrays type definition + methods

abstract type abstractGcmfaces end

struct gcmfaces <: abstractGcmfaces
   nFaces::Int
   grTopo::String
   f
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

# additional constructors

function gcmfaces(nFaces::Int,grTopo::String)
  tmp1=Array{Any}(undef,nFaces)
  gcmfaces(nFaces,grTopo,tmp1);
end

gcmfaces() = gcmfaces(MeshArrays.nFaces,MeshArrays.grTopo);

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

function getindex(a::gcmfaces,I::Vararg{Any, N}) where {N}
  nFaces=a.nFaces;
  grTopo=a.grTopo;
  v1=Array{Any}(undef,nFaces);
  for iFace=1:nFaces
    if N==2;
      v1[iFace]=getindex(a.f[iFace],I[1],I[2]);
    elseif N==3;
      v1[iFace]=getindex(a.f[iFace],I[1],I[2],I[3]);
    elseif N==4;
      v1[iFace]=getindex(a.f[iFace],I[1],I[2],I[3],I[4]);
    else;
      error("N>4 case not implemented yet")
    end
  end
  c=gcmfaces(nFaces,grTopo,v1);
  return c;
end

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
