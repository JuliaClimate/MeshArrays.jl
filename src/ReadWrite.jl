
import Base: read, write, read!

"""
    nFacesEtc(a::AbstractMeshArray)

Return nFaces, n3, n4 (1 by default)
"""
function nFacesEtc(a::AbstractMeshArray)
  nFaces=length(a.fIndex)
  ndims(a.f)>1 ? n3=size(a.f,2) : n3=1
  ndims(a.f)>2 ? n4=size(a.f,3) : n4=1
  return nFaces, n3, n4
end

"""
    read(fil::String,x::AbstractMeshArray)

Read array from file and return as a MeshArray. 

_The second argument (MeshArray or gcmgrid) provides the grid specifications (x.grid.ioSize)._
```
"""
function read(fil::String,x::AbstractMeshArray)

  (n1,n2)=x.grid.ioSize
  (nFaces,n3,n4)=nFacesEtc(x)

  fid = open(fil)
  xx = Array{eltype(x),1}(undef,n1*n2*n3*n4)
  read!(fid,xx)
  xx = reshape(hton.(xx),(n1,n2,n3,n4))
  close(fid)

  return x.grid.read(xx,x)
end

"""
    read(xx::Array,γ::gcmgrid)

Reformat Array data into a MeshArray shaped after `γ`.
"""
function read(xx::Array,γ::gcmgrid; verbose=false)
  siz=size.(Ref(xx),[1,2,3,4])
  S=[i for i in γ.ioSize[:]]
  verbose ? println(siz) : nothing
  verbose ? println(S) : nothing
  if siz[1:2]==S||siz[1:2]==transpose(S)
    n3=siz[3]; n4=siz[4]
    yy=reshape(xx,(γ.ioSize...,n3,n4))
  elseif siz[1]==prod(S)
    n3=siz[2]; n4=siz[3]
    yy=reshape(xx,(γ.ioSize...,n3,n4))
  elseif mod(siz[1],prod(S))==0
    n3=Int(siz[1]/prod(γ.ioSize[:]))
    n4=1
    yy=reshape(xx,(γ.ioSize...,n3,n4))
  else
    error("unexpected array size")
  end
  yy

  verbose ? println(size(yy)) : nothing

  if n3==1&&n4==1
    read(yy,MeshArray(γ,γ.ioPrec))
  elseif n4==1
    read(yy,MeshArray(γ,γ.ioPrec,n3))
  else
    read(yy,MeshArray(γ,γ.ioPrec,n3,n4))
  end
end

"""
    read(xx::Array,x::AbstractMeshArray)

Reformat Array data into a MeshArray similar to `x`. 
"""
function read(xx::Array,x::AbstractMeshArray)
  y=similar(x; m=x.meta)
  read!(xx,y)
  return y
end
  
"""
    read!(xx::Array,x::AbstractMeshArray)

Reformat array into MeshArray and write into `x`.
"""
function read!(xx::Array,x::AbstractMeshArray)
  facesSize=x.grid.fSize
  (n1,n2)=x.grid.ioSize
  (nFaces,n3,n4)=nFacesEtc(x)

  tmp=zeros(x.grid)
  for i3 in 1:n3
    for i4 in 1:n4
      read_one!(xx[:,:,i3,i4],tmp)
      for f in 1:nFaces
        if (n3>1)&&(n4>1)
          x[f,i3,i4].=tmp[f]
        elseif n3>1
          x[f,i3].=tmp[f]
        else
          x[f].=tmp[f]
        end
      end
    end
  end
end

"""
    read!(xx::Array,x::AbstractMeshArray)

Reformat one array of size x.grid.ioSize, and write **in-place** into MeshArray `x``.
"""
function read_one!(xx::Array,x::AbstractMeshArray; verbose=false)
  test1=in(x.grid.class,["PeriodicChannel","PeriodicDomain"])
  format=(test1 ? :simple : :compact)

  facesSize=x.grid.fSize
  (n1,n2)=x.grid.ioSize
  (nFaces,n3,n4)=nFacesEtc(x)
  i0=1; i1=0;
  j0=1; j1=0;
  for iFace=1:nFaces
    nn=facesSize[iFace][1]; mm=facesSize[iFace][2];
    if format==:compact
      i0=i1+1
      i1=i1+nn*mm
      x.f[iFace]=reshape(xx[:][i0:i1,:],(nn,mm))
    else
      i0=(mod(i1,n1)==0 ? 1 : i1+1)
      j0=(mod(i1,n1)==0&&iFace!==1 ? j0+mm : j0)
      i1=i0+nn-1; j1=j0+mm-1
      x.f[iFace]=reshape(xx[i0:i1,j0:j1],(nn,mm))
    end
end

end


"""
    write(fil::String,x::AbstractMeshArray)

Write MeshArray to binary file. Other methods:

```
write(xx::Array,x::AbstractMeshArray) #to Array
```
"""
function write(fil::String,x::AbstractMeshArray)
  y=x.grid.write(x)
  fid = open(fil,"w")
  write(fid,ntoh.(y))
  close(fid)
end

function write(x::AbstractMeshArray; verbose=false)
  test1=in(x.grid.class,["PeriodicChannel","PeriodicDomain"])
  format=(test1 ? :simple : :compact)

  facesSize=x.grid.fSize
  (n1,n2)=x.grid.ioSize
  (nFaces,n3,n4)=nFacesEtc(x)

  if format==:compact
    y = Array{eltype(x),3}(undef,(n1*n2,n3,n4))
  else
    y = Array{eltype(x),4}(undef,(n1,n2,n3,n4))
  end

  i0=1; i1=0;
  j0=1; j1=0;
  for iFace=1:nFaces;
    nn=facesSize[iFace][1];
    mm=facesSize[iFace][2];
    if format==:compact
      i0=i1+1
      i1=i1+nn*mm
    else
      i0=(mod(i1,n1)==0 ? 1 : i1+1)
      j0=(mod(i1,n1)==0&&iFace!==1 ? j0+mm : j0)
      i1=i0+nn-1; j1=j0+mm-1
      verbose ? println((i0,i1,j0,j1)) : false
    end
    for i4=1:n4
      for i3=1:n3
        if format==:compact
          y[i0:i1,i3,i4]=reshape(x.f[iFace,i3,i4],(nn*mm,1))
        else
          y[i0:i1,j0:j1,i3,i4]=reshape(x.f[iFace,i3,i4],(nn,mm))
        end
      end
    end
  end

  if format==:compact
    y=reshape(y,(n1,n2,n3,n4))
  else
    y=y
  end

  yy=if n3==1&&n4==1
    dropdims(y,dims=(3,4))
  elseif n4==1
    dropdims(y,dims=4)
  else
    y
  end
end

## 

function read_tiles(fil::String,x::AbstractMeshArray)

  (n1,n2)=x.grid.ioSize
  (nFaces,n3)=nFacesEtc(x)

  fid = open(fil)
  xx = Array{eltype(x),2}(undef,(n1*n2,n3))
  read!(fid,xx)
  xx = reshape(hton.(xx),(n1,n2,n3))
  close(fid)

  return x.grid.read(xx,x)

end

function read_tiles(xx::Array,γ::gcmgrid)
  n3=Int(round(prod(size(xx))/prod(γ.ioSize)))
  read_tiles(xx,MeshArray(γ,γ.ioPrec,n3))
end

function read_tiles(xx::Array,x::AbstractMeshArray)
	tmp=similar(x)
	s=x.grid.ioSize
	(n1,n2)=x.grid.fSize[1]
	ni=Int(s[1]/n1)
	nj=Int(s[2]/n2)
	ii=[0]
	for j in 1:nj, i in 1:ni
		ii[1]+=1
		tmp[ii[1]]=xx[(i-1)*n1 .+ collect(1:n1),(j-1)*n2 .+ collect(1:n2)]
	end
	tmp
end

function write_tiles(x::AbstractMeshArray)
  tmp = Array{eltype(x),2}(undef,x.grid.ioSize...)
	s=x.grid.ioSize
	(n1,n2)=x.grid.fSize[1]
	ni=Int(s[1]/n1)
	nj=Int(s[2]/n2)
	ii=[0]
	for j in 1:nj, i in 1:ni
		ii[1]+=1
		tmp[(i-1)*n1 .+ collect(1:n1),(j-1)*n2 .+ collect(1:n2)].=x[ii[1]]
	end
	tmp
end

function write_tiles(fil::String,x::AbstractMeshArray)
  y=x.grid.write(x)
  fid = open(fil,"w")
  write(fid,ntoh.(y))
  close(fid)
end

##
