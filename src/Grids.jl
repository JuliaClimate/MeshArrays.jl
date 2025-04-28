
GRID_LLC90_hash = artifact_hash("GRID_LLC90", artifact_toml)
GRID_LLC90 = joinpath(artifact_path(GRID_LLC90_hash)*"/","GRID_LLC90-1.1/")
GRID_LLC90_download() = artifact"GRID_LLC90"
GRID_LL360_hash = artifact_hash("GRID_LL360", artifact_toml)
GRID_LL360 = joinpath(artifact_path(GRID_LL360_hash)*"/","GRID_LL360-1.0/")
GRID_LL360_download() = artifact"GRID_LL360"
GRID_CS32_hash = artifact_hash("GRID_CS32", artifact_toml)
GRID_CS32 = joinpath(artifact_path(GRID_CS32_hash)*"/","GRID_CS32-1.1/")
GRID_CS32_download() = artifact"GRID_CS32"

Dict_to_NamedTuple(tmp::Dict) = (; zip(Symbol.(keys(tmp)), values(tmp))...)

"""
    land_mask(m::MeshArray)

Define land mask from `m` (1 if m>0; NaN if otherwise).
"""
function land_mask(m::MeshArray)
    μ=m
    μ[findall(μ.>0.0)].=1.0
    μ[findall(μ.==0.0)].=NaN
    μ
end

land_mask(Γ::NamedTuple)=land_mask(Γ.hFacC[:,1])

## GridSpec function with category argument:

"""
    GridSpec(category="PeriodicDomain",path=tempname(); ID=:unknown)

- Select one of the pre-defined grids either by ID (keyword) or by category.
- Return the corresponding `gcmgrid` specification, including the path where grid files can be accessed (`path`).

1. selection by `ID`

- `:LLC90`
- `:CS32`
- `:onedegree`
- `:default`

Example:

```
using MeshArrays
g = GridSpec(ID=:LLC90)
```

note : the path to these fully supported grids are handled internally in `MeshArrays.jl`.

2. by `category` and `path`

- `"PeriodicDomain"`
- `"PeriodicChannel"`
- `"CubeSphere"`
- `"LatLonCap"``

Examples:

```jldoctest; output = false
using MeshArrays
g = GridSpec()
g = GridSpec("PeriodicChannel",MeshArrays.GRID_LL360)
g = GridSpec("CubeSphere",MeshArrays.GRID_CS32)
g = GridSpec("LatLonCap",MeshArrays.GRID_LLC90)
isa(g,gcmgrid)

# output

true
```
"""
function GridSpec(category="PeriodicDomain",path=tempname(); ID=:unknown)

if category=="LatLonCap"
    nFaces=5
    grTopo="LatLonCap"
    ioSize=[90 1170]
    facesSize=[(90, 270), (90, 270), (90, 90), (270, 90), (270, 90)]
    ioPrec=Float64
elseif category=="CubeSphere"
    nFaces=6
    grTopo="CubeSphere"
    ioSize=[32 192]
    facesSize=[(32, 32), (32, 32), (32, 32), (32, 32), (32, 32), (32, 32)]
    ioPrec=Float32
elseif category=="PeriodicChannel"
    nFaces=1
    grTopo="PeriodicChannel"
    ioSize=[360 160]
    facesSize=[(360, 160)]
    ioPrec=Float32
elseif category=="PeriodicDomain"
    nFaces=4
    grTopo="PeriodicDomain"
    ioSize=[80 42]
    facesSize=[(40, 21), (40, 21), (40, 21), (40, 21)]
    ioPrec=Float32
else
    error("unknown category case")
end

if ID==:unknown
    gcmgrid(path,grTopo,nFaces,facesSize, ioSize, ioPrec, read, write)
elseif ID==:LLC90
    GridSpec("LatLonCap",MeshArrays.GRID_LLC90)
elseif ID==:CS32
    GridSpec("CubeSphere",MeshArrays.GRID_CS32)
elseif ID==:onedegree
    GridSpec("PeriodicChannel",MeshArrays.GRID_LL360)
elseif ID==:default
    GridSpec()
else
    error("unknwown grid")
end

end

## GridLoad function

"""
    GridLoad(γ=GridSpec(); ID=:default, option=:minimal)

- Return a `NamedTuple` of grid variables read from files located in `γ.path` (see `?GridSpec`).
- option : 
  - option=:minimal (default) to get only grid cell center positions (XC, YC). 
  - option=:light to get a complete set of 2D grid variables. 
  - option=:full  to get a complete set of 2D & 3D grid variables. 

Based on the MITgcm naming convention, grid variables are:

- XC, XG, YC, YG, AngleCS, AngleSN, hFacC, hFacS, hFacW, Depth.
- RAC, RAW, RAS, RAZ, DXC, DXG, DYC, DYG.
- DRC, DRF, RC, RF (one-dimensional)

MITgcm documentation : 

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

    gr.path==GRID_CS32 ? GRID_CS32_download() : nothing
    gr.path==GRID_LL360 ? GRID_LL360_download() : nothing
    gr.path==GRID_LLC90 ? GRID_LLC90_download() : nothing

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

γ = GridSpec("CubeSphere",MeshArrays.GRID_CS32)
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
    Tiles(γ::gcmgrid,ni::Int,nj::Int)

Define sudomain `tiles` of size `ni,nj`. Each tile is defined by a `Dict` where
`tile,face,i,j` correspond to tile ID, face ID, index ranges.

```jldoctest; output = false
using MeshArrays
γ=GridSpec("LatLonCap",MeshArrays.GRID_LLC90)
τ=Tiles(γ,30,30)

isa(τ[1],NamedTuple)

# output

true
```
"""
function Tiles(γ::gcmgrid,ni::Int,nj::Int)
    nt=Int(prod(γ.ioSize)/ni/nj)
    τ=Array{Dict,1}(undef,nt)
    #
    cnt=0
    for iF=1:γ.nFaces
        for jj=Int.(1:γ.fSize[iF][2]/nj)
            for ii=Int.(1:γ.fSize[iF][1]/ni)
                cnt=cnt+1
                i=(1:ni).+ni*(ii-1)
                j=(1:nj).+nj*(jj-1)
                τ[cnt]=Dict("tile" => cnt, "face" => iF, "i" => i, "j" => j)
            end
        end
    end
    #
    return Dict_to_NamedTuple.(τ)
end

"""
    Tiles(τ::Array{Dict},x::MeshArray)

Return an `Array` of tiles which cover `x` according to tile partition `τ`.

```jldoctest; output = false
using MeshArrays
γ=GridSpec("LatLonCap",MeshArrays.GRID_LLC90)
d=γ.read(γ.path*"Depth.data",MeshArray(γ,γ.ioPrec))
τ=Tiles(γ,30,30)
td=Tiles(τ,d)

D=similar(d)
Tiles!(τ,td,D)

isa(td[1],Array)

# output

true
```
"""
function Tiles(τ::Array,x::MeshArray)
    nt=length(τ)
    tx=Array{typeof(x[1]),1}(undef,nt)
    dn=size(x[1],1)-x.fSize[1][1]
    for ii=1:nt
        f=τ[ii].face
        i0=minimum(τ[ii].i)
        i1=maximum(τ[ii].i)
        j0=minimum(τ[ii].j)
        j1=maximum(τ[ii].j)
        tx[ii]=view(x[f],i0:i1+dn,j0:j1+dn)
    end
    return tx
end

"""
    Tiles!(τ::Array,tx::Array,x::MeshArrays)

Map tiles in `tx` according to tile partition `τ` into `x`.
"""
function Tiles!(T::Array,tx::Array,x::MeshArray)
    nt=length(T)
    dn=size(x[1],1)-x.fSize[1][1]
    for ii=1:nt
        f=T[ii].face
        i0=minimum(T[ii].i)
        i1=maximum(T[ii].i)
        j0=minimum(T[ii].j)
        j1=maximum(T[ii].j)
        x[f][i0:i1+dn,j0:j1+dn].=tx[ii]
    end
end

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
