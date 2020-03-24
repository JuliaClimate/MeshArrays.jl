
## GridSpec function with default GridName argument:

GridSpec() = GridSpec("LatLonCap","GRID_LLC90/")

## GridSpec function with GridName argument:

"""
    GridSpec(GridName,GridParentDir="./")

Return a `gmcgrid` specification that provides grid files `path`,
`class`, `nFaces`, `ioSize`, `facesSize`, `ioPrec`, & a `read` function
(not yet) using hard-coded values for `"PeriodicDomain"`, `"PeriodicChanel"`,
`"CubeSphere"`, and `"LatLonCap" for now.
"""
function GridSpec(GridName,GridParentDir="./")

grDir=GridParentDir
if GridName=="LatLonCap"
    nFaces=5
    grTopo="LatLonCap"
    ioSize=[90 1170]
    facesSize=[(90, 270), (90, 270), (90, 90), (270, 90), (270, 90)]
    ioPrec=Float64
elseif GridName=="CubeSphere"
    nFaces=6
    grTopo="CubeSphere"
    ioSize=[32 192]
    facesSize=[(32, 32), (32, 32), (32, 32), (32, 32), (32, 32), (32, 32)]
    ioPrec=Float32
elseif GridName=="PeriodicChannel"
    nFaces=1
    grTopo="PeriodicChannel"
    ioSize=[360 160]
    facesSize=[(360, 160)]
    ioPrec=Float32
elseif GridName=="PeriodicDomain"
    nFaces=4
    grTopo="PeriodicDomain"
    ioSize=[80 42]
    facesSize=[(40, 21), (40, 21), (40, 21), (40, 21)]
    ioPrec=Float32
else
    error("unknown GridName case")
end

mygrid=gcmgrid(grDir,grTopo,nFaces,facesSize, ioSize, ioPrec, read, write)

return mygrid

end

## GridLoad function

"""
    GridLoad(mygrid::gcmgrid)

Return a `Dict` of grid variables read from files located in `mygrid.path` (see `?GridSpec`).

Based on the MITgcm naming convention, grid variables are:

- XC, XG, YC, YG, AngleCS, AngleSN, hFacC, hFacS, hFacW, Depth.
- RAC, RAW, RAS, RAZ, DXC, DXG, DYC, DYG.
- DRC, DRF, RC, RF (one-dimensional)
"""
function GridLoad(mygrid::gcmgrid)

    GridVariables=Dict()

    list0=("XC","XG","YC","YG","AngleCS","AngleSN","RAC","RAW","RAS","RAZ",
    "DXC","DXG","DYC","DYG","Depth")
    for ii=1:length(list0)
        tmp1=mygrid.read(mygrid.path*list0[ii]*".data",MeshArray(mygrid,mygrid.ioPrec))
        tmp2=Symbol(list0[ii])
        @eval (($tmp2) = ($tmp1))
        GridVariables[list0[ii]]=tmp1
    end

    mygrid.ioPrec==Float64 ? reclen=8 : reclen=4

    list0=("DRC","DRF","RC","RF")
    for ii=1:length(list0)
        fil=mygrid.path*list0[ii]*".data"
        tmp1=stat(fil)
        n3=Int64(tmp1.size/reclen)

        fid = open(fil)
        tmp1 = Array{mygrid.ioPrec,1}(undef,n3)
        read!(fid,tmp1)
        tmp1 = hton.(tmp1)

        tmp2=Symbol(list0[ii])
        @eval (($tmp2) = ($tmp1))
        GridVariables[list0[ii]]=tmp1
    end

    list0=("hFacC","hFacS","hFacW")
    n3=length(GridVariables["RC"])
    for ii=1:length(list0)
        tmp1=mygrid.read(mygrid.path*list0[ii]*".data",MeshArray(mygrid,mygrid.ioPrec,n3))
        tmp2=Symbol(list0[ii])
        @eval (($tmp2) = ($tmp1))
        GridVariables[list0[ii]]=tmp1
    end

    return GridVariables

end

"""
    GridOfOnes(grTp,nF,nP)

Define all-ones grid variables instead of using `GridSpec` & `GridLoad`. E.g.

```
γ,Γ=GridOfOnes("CubeSphere",6,20);
```
"""
function GridOfOnes(grTp,nF,nP)

    grDir=""
    grTopo=grTp
    nFaces=nF
    if grTopo=="LatLonCap"
        ioSize=[nP nP*nF]
    elseif grTopo=="CubeSphere"
        ioSize=[nP nP*nF]
    elseif grTopo=="PeriodicChanel"
        ioSize=[nP nP]
    elseif grTopo=="PeriodicDomain"
        nFsqrt=Int(sqrt(nF))
        ioSize=[nP*nFsqrt nP*nFsqrt]
    end
    facesSize=Array{NTuple{2, Int},1}(undef,nFaces)
    facesSize[:].=[(nP,nP)]
    ioPrec=Float32

    mygrid=gcmgrid(grDir,grTopo,nFaces,facesSize, ioSize, ioPrec, read, write)

    GridVariables=Dict()
    list0=("XC","XG","YC","YG","RAC","RAZ","DXC","DXG","DYC","DYG","hFacC","hFacS","hFacW","Depth");
    for ii=1:length(list0);
        tmp1=fill(1.,nP,nP*nF);
        tmp1=mygrid.read(tmp1,MeshArray(mygrid,Float64));
        tmp2=Symbol(list0[ii]);
        @eval (($tmp2) = ($tmp1))
        GridVariables[list0[ii]]=tmp1
    end

    return mygrid, GridVariables

end

"""
    TileMap(ni::Int,nj::Int,mygrid::gcmgrid)

Return a `MeshArray` map of tile indices for tile size `ni,nj`
"""
function TileMap(mygrid::gcmgrid,ni::Int,nj::Int)
    nbr=MeshArray(mygrid)
    #
    cnt=0
    for iF=1:mygrid.nFaces
        for jj=Int.(1:mygrid.fSize[iF][2]/nj)
            for ii=Int.(1:mygrid.fSize[iF][1]/ni)
                cnt=cnt+1
                tmp_i=(1:ni).+ni*(ii-1)
                tmp_j=(1:nj).+nj*(jj-1)
                nbr.f[iF][tmp_i,tmp_j]=cnt*ones(Int,ni,nj)
            end
        end
    end
    #
    return nbr
end

"""
    GridAddWS!(GridVariables::Dict)

Compute XW, YW, XS, and YS (vector field locations) from XC, YC (tracer
field locations) and add them to GridVariables.

```
GridVariables=GridLoad(GridSpec("LatLonCap","GRID_LLC90/"))
GridAddWS!(GridVariables)
```
"""
function GridAddWS!(GridVariables::Dict)

    XC=exchange(GridVariables["XC"])
    YC=exchange(GridVariables["YC"])
    nFaces=XC.grid.nFaces
    XW=NaN .* XC; YW=NaN .* YC; XS=NaN .* XC; YS=NaN .* YC;

    for ff=1:nFaces
        tmp1=XC[ff][1:end-2,2:end-1]
        tmp2=XC[ff][2:end-1,2:end-1]
        tmp2[tmp2.-tmp1.>180]=tmp2[tmp2.-tmp1.>180].-360;
        tmp2[tmp1.-tmp2.>180]=tmp2[tmp1.-tmp2.>180].+360;
        XW[ff]=(tmp1.+tmp2)./2;
       #
        tmp1=XC[ff][2:end-1,1:end-2]
        tmp2=XC[ff][2:end-1,2:end-1]
        tmp2[tmp2.-tmp1.>180]=tmp2[tmp2.-tmp1.>180].-360;
        tmp2[tmp1.-tmp2.>180]=tmp2[tmp1.-tmp2.>180].+360;
        XS[ff]=(tmp1.+tmp2)./2;
       #
        tmp1=YC[ff][1:end-2,2:end-1]
        tmp2=YC[ff][2:end-1,2:end-1]
        YW[ff]=(tmp1.+tmp2)./2;
       #
        tmp1=YC[ff][2:end-1,1:end-2]
        tmp2=YC[ff][2:end-1,2:end-1]
        YS[ff]=(tmp1.+tmp2)./2;
    end;

    Xmax=180; Xmin=-180;
    XS[findall(XS.<Xmin)]=XS[findall(XS.<Xmin)].+360;
    XS[findall(XS.>Xmax)]=XS[findall(XS.>Xmax)].-360;
    XW[findall(XW.<Xmin)]=XW[findall(XW.<Xmin)].+360;
    XW[findall(XW.>Xmax)]=XW[findall(XW.>Xmax)].-360;

    GridVariables["XW"]=XW
    GridVariables["XS"]=XS
    GridVariables["YW"]=YW
    GridVariables["YS"]=YS
    return GridVariables
end
