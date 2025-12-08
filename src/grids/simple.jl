
module Grids_simple

using Unitful
import MeshArrays: gcmgrid, varmeta, MeshArray, Dict_to_NamedTuple
import MeshArrays: read_tiles, write_tiles, GridSpec_default

"""
    GridSpec_ones(grTp,nF,nP)

Define grid where all variables will be set 
to one via `GridLoad_ones`.

```
using MeshArrays
γ=MeshArrays.GridSpec_ones("CubeSphere",6,20);
Γ=MeshArrays.GridLoad_ones(γ;option="full")
```
"""
function GridSpec_ones(grTp,nF,nP)
    grDir="_ones"
    grTopo=grTp
    nFaces=nF
    if grTopo=="LatLonCap"
        ioSize=[nP nP*nF]
    elseif grTopo=="CubeSphere"
        ioSize=[nP nP*nF]
    elseif grTopo=="PeriodicChannel"
        ioSize=[nP nP]
    elseif grTopo=="PeriodicDomain"
        nFsqrt=Int(sqrt(nF))
        ioSize=[nP*nFsqrt nP*nFsqrt]
    else
        error("unknown configuration (grTopo)")
    end
    facesSize=Array{NTuple{2, Int},1}(undef,nFaces)
    facesSize[:].=[(nP,nP)]
    ioPrec=Float32

    gcmgrid(grDir,grTopo,nFaces,facesSize, ioSize, ioPrec, read, write)
end

"""
    GridLoad_ones(γ::gcmgrid;option="minimal")

Generate a unit grid, where every grid spacing and area is 1, according to γ. 
"""

"""
    GridLoad_ones(γ::gcmgrid;option="minimal")

Generate a unit grid, where every grid spacing and area is 1, 
as specified by `γ` obtained from `GridSpec_ones`. 

```
using MeshArrays
γ=MeshArrays.GridSpec_ones("CubeSphere",6,20);
Γ=MeshArrays.GridLoad_ones(γ;option="full")
```
"""
function GridLoad_ones(γ::gcmgrid;option="minimal")
    nFaces=γ.nFaces
    ioSize=(γ.ioSize[1],γ.ioSize[2])

    Γ=Dict()

    pc=fill(0.5,2); pg=fill(0.0,2); pu=[0.,0.5]; pv=[0.5,0.];
    if option=="full"
        list_n=("XC","XG","YC","YG","RAC","RAW","RAS","RAZ","DXC","DXG","DYC","DYG","Depth","hFacC","hFacS","hFacW");
        list_u=(u"m",u"m",u"m",u"m",u"m^2",u"m^2",u"m^2",u"m^2",u"m",u"m",u"m",u"m",u"m",1.0,1.0,1.0)
        list_p=(pc,pg,pc,pg,pc,pu,pv,pg,pu,pv,pv,pu,pc,fill(0.5,3),[0.,0.5,0.5],[0.5,0.,0.5])
    elseif option=="light"
        list_n=("XC","XG","YC","YG","RAC","DXC","DXG","DYC","DYG","Depth");
        list_u=(u"m",u"m",u"m",u"m",u"m^2",u"m",u"m",u"m",u"m",u"m")
        list_p=(pc,pg,pc,pg,pc,pu,pv,pv,pu,pc)
    else
        list_n=("XC","YC");
        list_u=(u"°",u"°")
        list_p=(pc,pc)
    end

    for ii=1:length(list_n);
        tmp1=fill(1.,(ioSize[:]))
        m=varmeta(list_u[ii],list_p[ii],missing,list_n[ii],list_n[ii]);
        tmp1=γ.read(tmp1,MeshArray(γ,Float64;meta=m));
        Γ[list_n[ii]]=tmp1
    end

    for i in 1:nFaces
        (np,nq)=γ.fSize[i]
        Γ["XC"][i]=vec(0.5:1.0:np-0.5)*ones(1,nq)
        option=="full" ? Γ["XG"][i]=vec(0.0:1.0:np-1.0)*ones(1,nq) : nothing
        Γ["YC"][i]=ones(np,1)*transpose(vec(0.5:1.0:nq-0.5))
        option=="full" ? Γ["YG"][i]=ones(np,1)*transpose(vec(0.0:1.0:nq-1.0)) : nothing
    end
    
    Dict_to_NamedTuple(Γ)
end

"""
    UnitGrid(ioSize, tileSize; option="minimal")
  
Generates a unit grid, where every grid spacing and area is 1, according to `ioSize, tileSize`,
and returns `Γ,γ`.

Since `ioSize, tileSize` defines a one to one mapping from global to tiled array, here we use 
`read_tiles, write_tiles` instead of the default `read, write`. And we overwrite local `XC,YC`
etc accordingly with global `XC,YC` etc . 
"""
function UnitGrid(ioSize::NTuple{2, Int},tileSize::NTuple{2, Int}; option="minimal")
    nF=div(prod(ioSize),prod(tileSize))
    fSize=fill(tileSize,nF)

    γ=gcmgrid("","PeriodicDomain",nF,fSize, ioSize, Float32, read_tiles, write_tiles)
    Γ=GridLoad_ones(γ;option=option)

    Γ.XC[:]=γ.read([i-0.5 for i in 1:ioSize[1], j in 1:ioSize[2]],γ)
    Γ.YC[:]=γ.read([j-0.5 for i in 1:ioSize[1], j in 1:ioSize[2]],γ)
    if option=="full"
        Γ.XG[:]=γ.read([i-1.0 for i in 1:ioSize[1], j in 1:ioSize[2]],γ)
        Γ.YG[:]=γ.read([j-1.0 for i in 1:ioSize[1], j in 1:ioSize[2]],γ)
    end

    return Γ,γ
end

#deprecated method
function GridOfOnes(grTp,nF,nP;option="minimal")
    γ=GridSpec_ones(grTp,nF,nP)
    γ,GridLoad_ones(γ;option=option)
end

"""
    periodic_domain(np::Integer,nq=missing)

Set up a simple periodic domain of size np x nq

```jldoctest; output = false
using MeshArrays
np=16 #domain size is np x np
Γ=MeshArrays.Grids_simple.periodic_domain(np)
isa(Γ.XC,MeshArray)

# output

true
```

"""
function periodic_domain(np::Integer,nq=missing)
    ismissing(nq) ? nq=np : nothing

    nFaces=1
    ioSize=[np nq]
    facesSize=Array{NTuple{2, Int},1}(undef,nFaces)
    facesSize[:].=[(np,nq)]
    ioPrec=Float32
    γ=gcmgrid("","PeriodicDomain",1,facesSize, ioSize, ioPrec, read, write)

    return GridLoad_ones(γ)
end

## Grid_latlon with specified vertical grid and land mask

rSphere = 6370.0*1000

"""
    GridLoad_lonlatdep(depth,mask)

!! deprecated, replaced with `grid_add_z` !!
"""
GridLoad_lonlatdep(depth,mask) = grid_add_z(GridLoad_lonlat(),depth,mask) 


"""
    grid_add_z(G,depth,mask)

```
import MeshArrays.Grids_simple: xy_IAP, grid_factors, grid_add_z

xy=xy_IAP()
gr=grid_factors(xy)
dep=[10 100 1000]; msk=ones(gr[:XC].fSize[1]...,3)
gr=grid_add_z(gr,dep,msk)
```
"""
grid_add_z(G,depth,mask) = begin
    RC=depth
    RF=[0;0.5*(RC[1:end-1]+RC[2:end]);RC[end]+0.5*(RC[end]-RC[end-1])]
    g=G.XC.grid
    hFacC=read(mask,g)
    merge(G,(hFacC=hFacC,RC=RC,RF=RF,DRF=diff(RF)))    
end

"""
    GridLoad_lonlat()

!! deprecated, replaced with `grid_factors` !!
"""
GridLoad_lonlat(xy=xy_IAP())=grid_factors(xy)

xy_IAP()=begin
    xg=0.0:1:360
    yg=-90.0:1:90
    xc=0.5:1:359.5
    yc=-89.5:1:89.5
    (xc=xc,yc=yc,xg=xg,yg=yg)
end

xy_Oscar()=begin
    xg=-0.125:0.25:359.875
    yg=-89.875:0.25:89.875
    xc=0.0:0.25:359.75
    yc=-89.75:0.25:89.75
    (xc=xc,yc=yc,xg=xg,yg=yg)
end

xy_OISST()=begin
    xg=0.0:0.25:360.0
    yg=-90.0:0.25:90.0
    xc=0.125:0.25:359.875
    yc=-89.875:0.25:89.875
    (xc=xc,yc=yc,xg=xg,yg=yg)
end

"""
    grid_factors(xy::NamedTuple)

```
import MeshArrays: Grids_simple

xy=Grids_simple.xy_OISST()
gr=Grids_simple.grid_factors(xy)
```
"""
grid_factors(xy::NamedTuple; nFaces=1)=begin
    (; xc, yc, xg, yg) = xy

    dx=diff(xg)[1]
    ni=length(xc)
    nj=length(yc)
    nni=Int(ni/sqrt(nFaces))
    nni=Int(nj/sqrt(nFaces))

    g=GridSpec_default(xy,nFaces)

    dxF = rSphere*deg2rad.(cosd.(yc)*dx)
    dyF = rSphere*deg2rad.(diff(yg))
    dxG = rSphere*deg2rad.(cosd.(yg[1:end-1])*dx)
    dyG = rSphere*deg2rad.(diff(yg))
    dxC = dxF
    dyC = rSphere*deg2rad.(diff(yg))
    x=sind.(yg[2:end])-sind.(yg[1:end-1])
    RAC = rSphere*rSphere*dx*deg2rad.(abs.(x))
  
    (
    XG=read(xg[1:end-1]*ones(1,nj),g),
    YG=read(permutedims(yc*ones(1,ni)),g),
    XC=read(xc*ones(1,nj),g),
    YC=read(permutedims(yc*ones(1,ni)),g),
    dxF=read(permutedims(dxF*ones(1,ni)),g),
    dyF=read(permutedims(dyF*ones(1,ni)),g),
    dxG=read(permutedims(dxG*ones(1,ni)),g),
    dyG=read(permutedims(dyG*ones(1,ni)),g),
    dxC=read(permutedims(dxC*ones(1,ni)),g),
    dyC=read(permutedims(dyC*ones(1,ni)),g),
    RAC=read(permutedims(RAC*ones(1,ni)),g),
    )
end

end
