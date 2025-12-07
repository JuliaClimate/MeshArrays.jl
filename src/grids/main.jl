
"""
    GridSpec(category="default",
        path=tempname(); np=nothing, ID=:unknown)

- Select one of the pre-defined grids 
    - or by `category` parameter
    - either by `ID` keyword (see `MeshArrays.GridSpec_MITgcm`)
- Return the corresponding `gcmgrid`
    - incl. a path to grid files (`path`).

```
using MeshArrays
γ=GridSpec()
```
"""
function GridSpec(category="default", 
        path=tempname(); np=nothing, ID=:unknown)
    if category=="default"&&ID==:unknown
        GridSpec_default(Grids_simple.xy_IAP())
    else
        GridSpec_MITgcm(category, path; np=np, ID=ID)
    end
end

function GridSpec_default(xy::NamedTuple, nFaces=1)
    (; xc, yc, xg, yg) = xy

    dx=diff(xg)[1]
    ni=length(xc)
    nj=length(yc)
    nni=Int(ni/sqrt(nFaces))
    nnj=Int(nj/sqrt(nFaces))

    grTopo="PeriodicChannel"
    ioSize=[ni nj]
    facesSize=[(nni, nnj), (nni, nnj), (nni, nnj), (nni, nnj)]
    ioPrec=Float32
    path=tempname()

    g=gcmgrid(path, grTopo, nFaces, facesSize, ioSize, ioPrec, read, write)
end

Dict_to_NamedTuple(tmp::Dict) = (; zip(Symbol.(keys(tmp)), values(tmp))...)

## GridLoad function

"""
    GridLoad(γ=GridSpec(); ID=:default, option=:minimal)

- Return a `NamedTuple` of grid variables read from files located in `γ.path` (see `?GridSpec`).
- option : 
  - option=:minimal (default) to get only grid cell center positions (XC, YC). 
  - option=:light to get a complete set of 2D grid variables. 
  - option=:full  to get a complete set of 2D & 3D grid variables. 
- if `ID` is specified then `GridSpec(ID=ID)` provides `γ`.

For grid variables, we follow the MITgcm naming convention.
Grid variables thus typically include :

- XC, XG, YC, YG, AngleCS, AngleSN, Depth
- RAC, RAW, RAS, RAZ, DXC, DXG, DYC, DYG
- three-dimensional : hFacC, hFacS, hFacW
- one-dimensional : DRC, DRF, RC, RF

For additional detail please refer to the MITgcm documentation : 

https://mitgcm.readthedocs.io/en/latest/algorithm/algorithm.html#spatial-discretization-of-the-dynamical-equations

```jldoctest; output = false
using MeshArrays
γ = GridSpec(ID=:LLC90)
Γ = GridLoad(γ;option="full")

isa(Γ.XC,MeshArray)

# output

true
```
"""
function GridLoad(γ=GridSpec(); ID=:default, option=:minimal)

    gr = (ID!==:default ? GridSpec(ID=ID) : γ)

    Γ=Dict()

    op=string(option)
    if op=="full"
        list_n=("XC","XG","YC","YG","RAC","RAW","RAS","RAZ","DXC","DXG","DYC","DYG","Depth")
        if (!isempty(filter(x -> occursin("AngleCS",x), readdir(gr.path))))
            list_n=(list_n...,"AngleCS","AngleSN");
        end
        list_n=(list_n...,"DRC","DRF","RC","RF")
        list_n=(list_n...,"hFacC","hFacS","hFacW")
    elseif op=="light"
        list_n=("XC","XG","YC","YG","RAC","DXC","DXG","DYC","DYG","Depth")
        if (!isempty(filter(x -> occursin("AngleCS",x), readdir(gr.path))))
            list_n=(list_n...,"AngleCS","AngleSN")
        end
        list_n=(list_n...,"DRC","DRF","RC","RF")
    elseif op=="minimal"||op=="minimum"
        list_n=("XC","YC")
    else
        error("unknown option")
    end

    [Γ[ii]=GridLoadVar(ii,gr) for ii in list_n]
    op=="full"||op=="light" ? GridAddWS!(Γ) : nothing
    return Dict_to_NamedTuple(Γ)
end

"""
    GridLoadVar(nam::String,γ::gcmgrid)

Return a grid variable read from files located in `γ.path` (see `?GridSpec`, `?GridLoad`).

Based on the MITgcm naming convention, grid variables are:

- XC, XG, YC, YG, AngleCS, AngleSN, hFacC, hFacS, hFacW, Depth.
- RAC, RAW, RAS, RAZ, DXC, DXG, DYC, DYG.
- DRC, DRF, RC, RF (one-dimensional)

MITgcm documentation : 

https://mitgcm.readthedocs.io/en/latest/algorithm/algorithm.html#spatial-discretization-of-the-dynamical-equations

```jldoctest; output = false
using MeshArrays

γ = GridSpec("CubeSphere",MeshArrays.Dataset("GRID_CS32"))
XC = GridLoadVar("XC",γ)

isa(XC,MeshArray)

# output

true
```
"""
function GridLoadVar(nam::String,γ::gcmgrid)
    pc=fill(0.5,2); pg=fill(0.0,2); pu=[0.,0.5]; pv=[0.5,0.];
    list_n=("XC","XG","YC","YG","RAC","RAW","RAS","RAZ","DXC","DXG","DYC","DYG","Depth","AngleCS","AngleSN")
    list_u=(u"°",u"°",u"°",u"°",u"m^2",u"m^2",u"m^2",u"m^2",u"m",u"m",u"m",u"m",u"m",1.0,1.0)
    list_p=(pc,pg,pc,pg,pc,pu,pv,pg,pu,pv,pv,pu,pc,pc,pc)
    #
    list3d_n=("hFacC","hFacS","hFacW");
    list3d_u=(1.0,1.0,1.0)
    list3d_p=(fill(0.5,3),[0.,0.5,0.5],[0.5,0.,0.5])
    #
    list1d_n=("DRC","DRF","RC","RF")

    if sum(nam.==list_n)==1
        ii=findall(nam.==list_n)[1]
        m=varmeta(list_u[ii],list_p[ii],missing,list_n[ii],list_n[ii])
        tmp1=γ.read(joinpath(γ.path,list_n[ii]*".data"),MeshArray(γ,γ.ioPrec;meta=m))
    elseif sum(nam.==list1d_n)==1
        fil=joinpath(γ.path,nam*".data")
        γ.ioPrec==Float64 ? reclen=8 : reclen=4
        n3=Int64(stat(fil).size/reclen)

        fid = open(fil)
        tmp1 = Array{γ.ioPrec,1}(undef,n3)
        read!(fid,tmp1)
        tmp1 = hton.(tmp1)
    elseif sum(nam.==list3d_n)==1
        fil=joinpath(γ.path,"RC.data")
        γ.ioPrec==Float64 ? reclen=8 : reclen=4
        n3=Int64(stat(fil).size/reclen)

        ii=findall(nam.==list3d_n)[1]
        m=varmeta(list3d_u[ii],list3d_p[ii],missing,list3d_n[ii],list3d_n[ii]);
        tmp1=γ.read(joinpath(γ.path,list3d_n[ii]*".data"),MeshArray(γ,γ.ioPrec,n3;meta=m))
    else
        tmp1=missing
    end
end

##

"""
    GridAddWS!(Γ::Dict)

Compute XW, YW, XS, and YS (vector field locations) from XC, YC (tracer
field locations) and add them to Γ.
"""
function GridAddWS!(Γ::Dict)

    XC=exchange(Γ["XC"]).MA
    YC=exchange(Γ["YC"]).MA
    nFaces=XC.grid.nFaces
    uX=XC.meta.unit
    uY=YC.meta.unit

    XW=MeshArrays.gcmarray(XC.grid,eltype(XC);meta=varmeta(uX,[0.,0.5],missing,"XW","XW"))
    YW=MeshArrays.gcmarray(XC.grid,eltype(XC);meta=varmeta(uY,[0.,0.5],missing,"YW","YW"))
    XS=MeshArrays.gcmarray(XC.grid,eltype(XC);meta=varmeta(uX,[0.5,0.],missing,"XS","XS"))
    YS=MeshArrays.gcmarray(XC.grid,eltype(XC);meta=varmeta(uY,[0.5,0.],missing,"YS","YS"))

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

    Γ["XW"]=XW
    Γ["XS"]=XS
    Γ["YW"]=YW
    Γ["YS"]=YS
    return Γ
end

##

include("simple.jl")
include("tiles.jl")
include("MITgcm.jl")
include("NEMO.jl")

##

function GridLoad_default(γ=GridSpec())
    xy=Grids_simple.xy_IAP()
    gr=Grids_simple.grid_factors(xy)

    dep=[10 100 1000]; msk=ones(gr[:XC].fSize[1]...,3)
    gr=Grids_simple.grid_add_z(gr,dep,msk)
end
