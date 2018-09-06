
using NetCDF, Printf;

export read_nctiles

## read_nctiles function
#
#examples:
#  fileName="nctiles_grid/GRID"
#  Depth=read_nctiles(fileName,"Depth")
#  hFacC=read_nctiles(fileName,"hFacC")

function read_nctiles(fileName,fldName);

if ~(MeshArrays.grTopo=="llc");
  error("non-llc cases not implemented yet");
end;

fileIn=@sprintf("%s.%04d.nc",fileName,1);
x = ncread(fileIn,fldName);
ndims=length(size(x));

#initialize f
if ndims==2;
  f0=Array{Float64}(undef,90,0);
  f00=Array{Float64}(undef,0,90);
else;
  nr=size(x,3)
  f0=Array{Float64}(undef,90,0,nr);
  f00=Array{Float64}(undef,0,90,nr);
end;
f=[f0,f0,NaN,f00,f00];

#fill in f
for ff=1:13;
  #read one tile
  fileIn=@sprintf("%s.%04d.nc",fileName,ff);
  x = ncread(fileIn,fldName);
  #combine tiles
  if ff<=3;
    f[1]=cat(f[1],x;dims=2);
  elseif ff<=6;
    f[2]=cat(f[2],x;dims=2);
  elseif ff==7;
    f[3]=x;
  elseif ff<=10;
    f[4]=cat(f[4],x;dims=1);
  elseif ff<=13;
    f[5]=cat(f[5],x;dims=1);
  end;

end;

#return f gcmfaces object fld
fld=gcmfaces(5,"llc",f);

end
