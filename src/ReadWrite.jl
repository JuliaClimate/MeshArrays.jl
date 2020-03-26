
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

  y=similar(x)
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
