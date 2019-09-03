
using NetCDF, Printf;

export read_nctiles

## read_nctiles function
#
#examples:
#  mygrid=GCMGridSpec("LLC90")
#  fileName="nctiles_grid/GRID"
#  Depth=read_nctiles(fileName,"Depth",mygrid)
#  hFacC=read_nctiles(fileName,"hFacC",mygrid)

"""
    read_nctiles(fileName,fldName,mygrid)

Read model output from Netcdf / NCTiles file and convert to MeshArray instance.
"""
function read_nctiles(fileName::String,fldName::String,mygrid::gcmgrid)

if (mygrid.class!="llc")||(mygrid.ioSize!=[90 1170])
  error("non-llc90 cases not implemented yet");
end;

fileIn=@sprintf("%s.%04d.nc",fileName,1);
x = ncread(fileIn,fldName);
ndims=length(size(x));

#initialize f
if ndims==2;
  f0=Array{Float64}(undef,90,0);
  f00=Array{Float64}(undef,0,90);
elseif ndims==3;
  f0=Array{Float64}(undef,90,0,size(x,3));
  f00=Array{Float64}(undef,0,90,size(x,3));
elseif ndims==4;
  f0=Array{Float64}(undef,90,0,size(x,3),size(x,4));
  f00=Array{Float64}(undef,0,90,size(x,3),size(x,4));
end;
f=[f0,f0,f0,f00,f00];

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

fld=MeshArray(mygrid,f)
return fld

end


function findtiles(ni,nj,grid="llc90")
  mytiles = Dict()
  GCMGridSpec()
  if grid=="llc90"
      GCMGridLoad()
  else
      println("Unsupported grid option")
  end
  mytiles["nFaces"]=MeshArrays.nFaces;
  #mytiles.fileFormat=mygrid.fileFormat;
  mytiles["ioSize"]=MeshArrays.ioSize;
  %
  XC=MeshArrays.XC;
  YC=MeshArrays.YC;
  XC11=copy(XC); YC11=copy(XC);
  XCNINJ=copy(XC); YCNINJ=copy(XC);
  iTile=copy(XC); jTile=copy(XC); tileNo=copy(XC);
  tileCount=0;
  for iF=1:XC11.nFaces;
      #global tileCount,XC,YC,XC11,YC11,iTile,jTile,tileNo
      face_XC=XC.f[iF]; face_YC=YC.f[iF];
  #ordering convention that was used in first generation nctile files:
  #    for ii=1:size(face_XC,1)/ni;
  #        for jj=1:size(face_XC,2)/nj;
  #ordering convention that is consistent with MITgcm/pkg/exch2:
      for jj=Int.(1:size(face_XC,2)/nj);
          for ii=Int.(1:size(face_XC,1)/ni);    
              tileCount=tileCount+1;
              tmp_i=(1:ni).+ni*(ii-1)
              tmp_j=(1:nj).+nj*(jj-1)
              tmp_XC=face_XC[tmp_i,tmp_j]
              tmp_YC=face_YC[tmp_i,tmp_j]
              XC11.f[iF][tmp_i,tmp_j].=tmp_XC[1,1]
              YC11.f[iF][tmp_i,tmp_j].=tmp_YC[1,1]
              XCNINJ.f[iF][tmp_i,tmp_j].=tmp_XC[end,end]
              YCNINJ.f[iF][tmp_i,tmp_j].=tmp_YC[end,end]
              iTile.f[iF][tmp_i,tmp_j]=collect(1:ni)*ones(Int,1,nj)
              jTile.f[iF][tmp_i,tmp_j]=ones(Int,ni,1)*collect(1:nj)'
              tileNo.f[iF][tmp_i,tmp_j]=tileCount*ones(Int,ni,nj)
          end
      end
  end

  mytiles["XC"] = XC;
  mytiles["YC"] = YC;
  mytiles["XC11"] = XC11;
  mytiles["YC11"] = YC11;
  mytiles["XCNINJ"] = XCNINJ;
  mytiles["YCNINJ"] = YCNINJ;
  mytiles["iTile"] = iTile;
  mytiles["jTile"] = jTile;
  mytiles["tileNo"] = tileNo;

  return mytiles

end
