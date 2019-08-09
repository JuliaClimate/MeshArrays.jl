
## read_bin function with full list of argument

"""
    read_bin(fil::String,kt::Union{Int,Missing},kk::Union{Int,Missing},prec::DataType,mygrid::gcmgrid)

Read model output from binary file and convert to gcmfaces structure. Other methods:

```
read_bin(fil::String,prec::DataType,mygrid::gcmgrid)
read_bin(fil::String,mygrid::gcmgrid)
```
"""
function read_bin(fil::String,kt::Union{Int,Missing},kk::Union{Int,Missing},prec::DataType,mygrid::gcmgrid)

  if ~ismissing(kt);
    error("non-empty kt option not implemented yet");
  end;

  if ~ismissing(kk);
    error("non-empty kk option not implemented yet");
  end;

  (n1,n2)=mygrid.ioSize

  if prec==Float64;
    reclen=8;
  else;
    reclen=4;
  end;
  tmp1=stat(fil);
  n3=Int64(tmp1.size/n1/n2/reclen);

  fid = open(fil);
  fld = Array{prec,1}(undef,(n1*n2*n3));
  read!(fid,fld);
  fld = hton.(fld);
  close(fid)

  n3>1 ? s=(n1,n2,n3) : s=(n1,n2)
  v0=reshape(fld,s);

  convert2gcmfaces(v0,mygrid)

end

## read_bin with reduced list of argument

# read_bin(fil::String,prec::DataType,mygrid::gcmgrid)
function read_bin(fil::String,prec::DataType,mygrid::gcmgrid)
  read_bin(fil,missing,missing,prec,mygrid::gcmgrid)
end

# read_bin(fil::String,mygrid::gcmgrid)
function read_bin(fil::String,mygrid::gcmgrid)
  read_bin(fil,missing,missing,mygrid.ioPrec,mygrid::gcmgrid)
end

## read_bin with alternative arguments

# read_bin(fil::String,x::gcmfaces)
function read_bin(fil::String,x::gcmfaces)
  read_bin(fil,missing,missing,x.mygrid.ioPrec,x.mygrid::gcmgrid)
end

# read_bin(tmp::Array,mygrid::gcmgrid)
function read_bin(tmp::Array,mygrid::gcmgrid)
  convert2gcmfaces(tmp,mygrid)
end

# read_bin(tmp::Array,x::gcmfaces)
function read_bin(tmp::Array,x::gcmfaces)
  convert2gcmfaces(tmp,x.mygrid)
end
