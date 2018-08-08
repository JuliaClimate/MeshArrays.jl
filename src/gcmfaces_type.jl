# gcmfaces_type.jl
#
#	First Draft Implementation
#
# gaelforget (https://github.com/gaelforget/gcmfaces_jl)
# Julia 0.6.2
# Created: 02.01.18
# Last Edit: 02.01.18

## original gcmfaces type definition + methods ##

abstract type abstractGcmfaces end

mutable struct gcmfaces <: abstractGcmfaces
   nFaces::Int
   grTopo::String
   f
end

# Custom pretty-printing

function Base.show(io::IO, z::gcmfaces)

    @printf io " GCMFaces instance with \n"
    print_with_color(:normal,io, "  grid type  = ")
    print_with_color(:blue, io, "$(z.grTopo)\n")
    print_with_color(:normal,io, "  data type  = ")
    print_with_color(:blue, io, "$(typeof(z.f[1][1]))\n")
    print_with_color(:normal,io, "  # of faces = ")
    print_with_color(:blue, io, "$(z.nFaces)\n")
    print_with_color(:normal,io, "  face sizes = ")
    print_with_color(:blue, io, "$(size(z.f[1]))\n")
    for iFace=2:z.nFaces
    print_with_color(:normal,io, "               ")
    print_with_color(:blue, io, "$(size(z.f[iFace]))\n")
    end

    return
end

# additional constructors

gcmfaces(nFaces::Int,grTopo::String) = gcmfaces(nFaces,grTopo,[NaN,NaN,NaN,NaN,NaN])

gcmfaces() = gcmfaces(5,"llc")

# basic operations

import Base: +, -, *, /, getindex

function +(a::gcmfaces,b::gcmfaces)
  c=a
  c.f=a.f+b.f
  return c
end

function +(a::Any,b::gcmfaces)
  c=b
  c.f=a+b.f
  return c
end

function +(a::gcmfaces,b::Any)
  c=a
  c.f=a.f+b
  return c
end

function -(a::gcmfaces,b::gcmfaces)
  c=a
  c.f=a.f-b.f
  return c
end

function -(a::Any,b::gcmfaces)
  c=b
  c.f=a-b.f
  return c
end

function -(a::gcmfaces,b::Any)
  c=a
  c.f=a.f-b
  return c
end

function *(a::gcmfaces,b::gcmfaces)
  #this defines * (and .*) as .* applied face by face
  nFaces=a.nFaces;
  grTopo=a.grTopo;
  v1=Array{Any}(nFaces);
  for iFace=1:nFaces
    v1[iFace]=a.f[iFace].*b.f[iFace];
  end
  c=gcmfaces(nFaces,grTopo,v1);
  return c
end

function *(a::Any,b::gcmfaces)
  nFaces=b.nFaces;
  grTopo=b.grTopo;
  v1=a*b.f;
  c=gcmfaces(nFaces,grTopo,v1);
  return c
end

function *(a::gcmfaces,b::Any)
  nFaces=a.nFaces;
  grTopo=a.grTopo;
  v1=b*a.f;
  c=gcmfaces(nFaces,grTopo,v1);
  return c
end

function /(a::gcmfaces,b::gcmfaces)
  #this defines / (and ./) as ./ applied face by face
  nFaces=a.nFaces;
  grTopo=a.grTopo;
  v1=Array{Any}(nFaces);
  for iFace=1:nFaces
    v1[iFace]=a.f[iFace]./b.f[iFace];
  end
  c=gcmfaces(nFaces,grTopo,v1);
  return c
end

function /(a::Any,b::gcmfaces)
  #this defines / (and ./) as ./ applied face by face
  nFaces=b.nFaces;
  grTopo=b.grTopo;
  v1=Array{Any}(nFaces);
  for iFace=1:nFaces
    v1[iFace]=a./b.f[iFace];
  end
  c=gcmfaces(nFaces,grTopo,v1);
  return c
end

function /(a::gcmfaces,b::Any)
  nFaces=a.nFaces;
  grTopo=a.grTopo;
  v1=a.f/b;
  c=gcmfaces(nFaces,grTopo,v1);
  return c
end

function getindex(a::gcmfaces,I::Vararg{Any, N}) where {N}
  nFaces=a.nFaces;
  grTopo=a.grTopo;
  v1=Array{Any}(nFaces);
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

