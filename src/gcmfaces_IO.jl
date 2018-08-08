# gcmfaces_IO.jl
#
#       First Draft Implementation; contents:
#          read_nctiles, read_bin, grid_load
#
# gaelforget (https://github.com/gaelforget/gcmfaces_jl)
# Julia 0.6.2
# Created: 02.01.18
# Last Edit: 02.01.18

## read_nctiles ##

function read_nctiles(fileName,fldName);
#usage: fld=read_nctiles(fileName,fldName);        reads full field (all depths, all times)

fileIn=@sprintf("%s.%04d.nc",fileName,1);
x = ncread(fileIn,fldName);
ndims=length(size(x));

#initialize f
if ndims==2;
  f0=Array{Float64}(90,0);
  f00=Array{Float64}(0,90);
else;
  nr=size(x,3)
  f0=Array{Float64}(90,0,nr);
  f00=Array{Float64}(0,90,nr);
end;
f=[f0,f0,NaN,f00,f00];

#fill in f
for ff=1:13;
  #read one tile
  fileIn=@sprintf("%s.%04d.nc",fileName,ff);
  x = ncread(fileIn,fldName);
  #combine tiles
  if ff<=3;
    f[1]=cat(2,f[1],x);
  elseif ff<=6;
    f[2]=cat(2,f[2],x);
  elseif ff==7;
    f[3]=x;
  elseif ff<=10;
    f[4]=cat(1,f[4],x);
  elseif ff<=13;
    f[5]=cat(1,f[5],x);
  end;

end;

#return f gcmfaces object fld
fld=gcmfaces(5,"llc",f);

end

## read_bin ##

function read_bin(fil::String,kt,kk,prec::DataType);

  if ~isempty(kt);
    error("non-empty kt option not implemented yet");
  end;

  if ~isempty(kk);
    error("non-empty kk option not implemented yet");
  end;

  nFaces=GCMFaces.nFaces;
  grTopo=GCMFaces.grTopo;
  (n1,n2)=GCMFaces.ioSize;
  facesSize=GCMFaces.facesSize;

  if prec==Float64;
    reclen=8;
  else;
    reclen=4;
  end;

  tmp1=stat(fil);
  n3=Int64(tmp1.size/n1/n2/reclen);

  fid = open(fil);
  fld = read(fid,prec,(n1*n2*n3));
  fld = hton.(fld);

  # the following corresponds to convert2gcmfaces.m, and
  # should later be moved to gcmfaces_convert.jl

  v0=reshape(fld,(n1*n2,n3));
  v1=Array{Any}(nFaces);
  i0=0; i1=0;
  for iFace=1:nFaces
    i0=i1+1;
    nn=facesSize[iFace,1]; mm=facesSize[iFace,2];
    i1=i1+nn*mm;
    if n3>1;
      v1[iFace]=reshape(v0[i0:i1,:],(nn,mm,n3));
    else;
      v1[iFace]=reshape(v0[i0:i1,:],(nn,mm));
    end;
  end

  gcmfaces(nFaces,grTopo,v1);
end

function read_bin(fil::String,prec::DataType);
  read_bin(fil,[],[],prec);
end

function read_bin(fil::String);
  read_bin(fil,[],[],Float32);
end
