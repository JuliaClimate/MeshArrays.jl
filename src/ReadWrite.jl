
import Base: read, write, read!

"""
    read(fil::String,x::MeshArray)

Read array from file and return as a MeshArray. 

_The second argument (MeshArray or gcmgrid) provides the grid specifications (x.grid.ioSize)._
```
"""
function read(fil::String,x::MeshArray)

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
function read(xx::Array,γ::gcmgrid)
  siz=size.(Ref(xx),[1,2,3,4])
  if siz[1:2]==[i for i in γ.ioSize]
    n3=siz[3]; n4=siz[4]
    yy=reshape(xx,(γ.ioSize...,n3,n4))
  elseif siz[1]==prod(γ.ioSize[:])
    n3=siz[2]; n4=siz[3]
    yy=reshape(xx,(γ.ioSize...,n3,n4))
  elseif mod(siz[1],prod(γ.ioSize[:]))==0
    n3=Int(siz[1]/prod(γ.ioSize[:]))
    n4=1
    yy=reshape(xx,(γ.ioSize...,n3,n4))
  else
    error("unexpected array size")
  end
  yy

  if n3==1&&n4==1
    read(yy,MeshArray(γ,γ.ioPrec))
  elseif n4==1
    read(yy,MeshArray(γ,γ.ioPrec,n3))
  else
    read(yy,MeshArray(γ,γ.ioPrec,n3,n4))
  end
end

"""
    read(xx::Array,x::MeshArray)

Reformat Array data into a MeshArray similar to `x`. 
"""
function read(xx::Array,x::MeshArray)
  y=similar(x; m=x.meta)
  read!(xx,y)
  return y
end
  
"""
    read!(xx::Array,x::MeshArray)

Reformat array into MeshArray and write into `x`.
"""
function read!(xx::Array,x::MeshArray)
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
    read!(xx::Array,x::MeshArray)

Reformat one array of size x.grid.ioSize, and write **in-place** into MeshArray `x``.
"""
function read_one!(xx::Array,x::MeshArray)
  facesSize=x.grid.fSize
  (n1,n2)=x.grid.ioSize
  (nFaces,n3,n4)=nFacesEtc(x)
  i0=0; i1=0;
  for iFace=1:nFaces
    i0=i1+1;
    nn=facesSize[iFace][1]; mm=facesSize[iFace][2];
    i1=i1+nn*mm;
    x.f[iFace]=reshape(xx[:][i0:i1,:],(nn,mm))
end

end


"""
    write(fil::String,x::MeshArray)

Write MeshArray to binary file. Other methods:

```
write(xx::Array,x::MeshArray) #to Array
```
"""
function write(fil::String,x::MeshArray)
  y=x.grid.write(x)
  fid = open(fil,"w")
  write(fid,ntoh.(y))
  close(fid)
end

function write(x::MeshArray)
  facesSize=x.grid.fSize
  (n1,n2)=x.grid.ioSize
  (nFaces,n3,n4)=nFacesEtc(x)

  y = Array{eltype(x),3}(undef,(n1*n2,n3,n4))
  i0=0; i1=0;
  for iFace=1:nFaces;
    i0=i1+1;
    nn=facesSize[iFace][1];
    mm=facesSize[iFace][2];
    i1=i1+nn*mm;
    for i4=1:n4
      for i3=1:n3
        y[i0:i1,i3,i4]=reshape(x.f[iFace,i3,i4],(nn*mm,1))
      end
    end
  end

  y=reshape(y,(n1,n2,n3,n4));
  yy=if n3==1&&n4==1
    dropdims(y,dims=(3,4))
  elseif n4==1
    dropdims(y,dims=4)
  else
    y
  end
end

## 

function read_tiles(fil::String,x::MeshArray)

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

function read_tiles(xx::Array,x::MeshArray)
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

function write_tiles(x::MeshArray)
  tmp = Array{eltype(x),2}(undef,x.grid.ioSize)
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

function write_tiles(fil::String,x::MeshArray)
  y=x.grid.write(x)
  fid = open(fil,"w")
  write(fid,ntoh.(y))
  close(fid)
end

##