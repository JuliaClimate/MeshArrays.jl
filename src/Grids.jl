
Dict_to_NamedTuple(tmp::Dict) = (; zip(Symbol.(keys(tmp)), values(tmp))...)

"""
    UnitGrid(γ::gcmgrid;option="minimal")

Generate a unit grid, where every grid spacing and area is 1, according to γ. 
"""
function UnitGrid(γ::gcmgrid;option="minimal")
    nFaces=γ.nFaces
    ioSize=(γ.ioSize[1],γ.ioSize[2])

    Γ=Dict()

    pc=fill(0.5,2); pg=fill(0.0,2); pu=[0.,0.5]; pv=[0.5,0.];
    if option=="full"
        list_n=("XC","XG","YC","YG","RAC","RAW","RAS","RAZ","DXC","DXG","DYC","DYG","Depth","hFacC","hFacS","hFacW");
        list_u=(u"m",u"m",u"m",u"m",u"m^2",u"m^2",u"m^2",u"m^2",u"m",u"m",u"m",u"m",u"m",1.0,1.0,1.0)
        list_p=(pc,pg,pc,pg,pc,pu,pv,pg,pu,pv,pv,pu,pc,fill(0.5,3),[0.,0.5,0.5],[0.5,0.,0.5])
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

    return Dict_to_NamedTuple(Γ)

end

"""
    simple_periodic_domain(np::Integer,nq=missing)

Set up a simple periodic domain of size np x nq

```jldoctest
using MeshArrays
np=16 #domain size is np x np
Γ=simple_periodic_domain(np)
isa(Γ.XC,MeshArray)

# output

true
```

"""
function simple_periodic_domain(np::Integer,nq=missing)
    ismissing(nq) ? nq=np : nothing

    nFaces=1
    ioSize=[np nq]
    facesSize=Array{NTuple{2, Int},1}(undef,nFaces)
    facesSize[:].=[(np,nq)]
    ioPrec=Float32
    γ=gcmgrid("","PeriodicDomain",1,facesSize, ioSize, ioPrec, read, write)

    return UnitGrid(γ)
end

## GridSpec function with default GridName argument:

GridSpec() = GridSpec("PeriodicDomain","./")

## GridSpec function with GridName argument:

"""
    GridSpec(GridName,GridParentDir="./")

Select one of the pre-defined grids (by `GridName`) and return 
the corresponding `gmcgrid` -- a global grid specification 
which contains the grid files location (`GridParentDir`).
    

Possible choices for `GridName`:

- `"PeriodicDomain"`
- `"PeriodicChannel"`
- `"CubeSphere"`
- `"LatLonCap"``

```jldoctest
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

return gcmgrid(grDir,grTopo,nFaces,facesSize, ioSize, ioPrec, read, write)

end

## GridLoad function

"""
    GridLoad(γ::gcmgrid;option="minimal")

Return a `NamedTuple` of grid variables read from files located in `γ.path` (see `?GridSpec`).

By default, option="minimal" means that only grid cell center positions (XC, YC) are loaded. 

The "full" option provides a complete set of grid variables. 

Based on the MITgcm naming convention, grid variables are:

- XC, XG, YC, YG, AngleCS, AngleSN, hFacC, hFacS, hFacW, Depth.
- RAC, RAW, RAS, RAZ, DXC, DXG, DYC, DYG.
- DRC, DRF, RC, RF (one-dimensional)

```jldoctest
using MeshArrays

γ = GridSpec("CubeSphere",MeshArrays.GRID_CS32)
#γ = GridSpec("LatLonCap",MeshArrays.GRID_LLC90)
#γ = GridSpec("PeriodicChannel",MeshArrays.GRID_LL360)

Γ = GridLoad(γ;option="full")

isa(Γ.XC,MeshArray)

# output

true
```
"""
function GridLoad(γ::gcmgrid;option="minimal")

    Γ=Dict()

    pc=fill(0.5,2); pg=fill(0.0,2); pu=[0.,0.5]; pv=[0.5,0.];
    if option=="full"
        list_n=("XC","XG","YC","YG","RAC","RAW","RAS","RAZ","DXC","DXG","DYC","DYG","Depth");
        list_u=(u"°",u"°",u"°",u"°",u"m^2",u"m^2",u"m^2",u"m^2",u"m",u"m",u"m",u"m",u"m")
        list_p=(pc,pg,pc,pg,pc,pu,pv,pg,pu,pv,pv,pu,pc)
    else
        list_n=("XC","YC");
        list_u=(u"°",u"°")
        list_p=(pc,pc)
    end

    if (option=="full")&&(!isempty(filter(x -> occursin("AngleCS",x), readdir(γ.path))))
        list_n=(list_n...,"AngleCS","AngleSN");
        list_u=(list_u...,1.0,1.0)
        list_p=(list_p...,pc,pc)
    end

    for ii=1:length(list_n)
        m=varmeta(list_u[ii],list_p[ii],missing,list_n[ii],list_n[ii])
        tmp1=γ.read(joinpath(γ.path,list_n[ii]*".data"),MeshArray(γ,γ.ioPrec;meta=m))
        tmp2=Symbol(list_n[ii])
        @eval (($tmp2) = ($tmp1))
        Γ[list_n[ii]]=tmp1
    end

    γ.ioPrec==Float64 ? reclen=8 : reclen=4

    if option=="full"
        list_n=("DRC","DRF","RC","RF")
    else
        list_n=()
    end
    for ii=1:length(list_n)
        fil=joinpath(γ.path,list_n[ii]*".data")
        tmp1=stat(fil)
        n3=Int64(tmp1.size/reclen)

        fid = open(fil)
        tmp1 = Array{γ.ioPrec,1}(undef,n3)
        read!(fid,tmp1)
        tmp1 = hton.(tmp1)

        tmp2=Symbol(list_n[ii])
        @eval (($tmp2) = ($tmp1))
        Γ[list_n[ii]]=tmp1
    end

    f=readdir(γ.path)
    if (option=="full")&&(sum(occursin.("hFacC",f))>0)
        list_n=("hFacC","hFacS","hFacW");
        list_u=(1.0,1.0,1.0)
        list_p=(fill(0.5,3),[0.,0.5,0.5],[0.5,0.,0.5])
        n3=length(Γ["RC"])
        for ii=1:length(list_n)
            m=varmeta(list_u[ii],list_p[ii],missing,list_n[ii],list_n[ii]);
            tmp1=γ.read(joinpath(γ.path,list_n[ii]*".data"),MeshArray(γ,γ.ioPrec,n3;meta=m))
            tmp2=Symbol(list_n[ii])
            @eval (($tmp2) = ($tmp1))
            Γ[list_n[ii]]=tmp1
        end
    end

    option=="full" ? GridAddWS!(Γ) : nothing

    return Dict_to_NamedTuple(Γ)

end

"""
    GridOfOnes(grTp,nF,nP;option="minimal")

Define all-ones grid variables instead of using `GridSpec` & `GridLoad`. E.g.

```
γ,Γ=GridOfOnes("CubeSphere",6,20);
```
"""
function GridOfOnes(grTp,nF,nP;option="minimal")

    grDir=""
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
    end
    facesSize=Array{NTuple{2, Int},1}(undef,nFaces)
    facesSize[:].=[(nP,nP)]
    ioPrec=Float32

    γ=gcmgrid(grDir,grTopo,nFaces,facesSize, ioSize, ioPrec, read, write)

    return γ, UnitGrid(γ;option=option)
end

"""
    Tiles(γ::gcmgrid,ni::Int,nj::Int)

Define sudomain `tiles` of size `ni,nj`. Each tile is defined by a `Dict` where
`tile,face,i,j` correspond to tile ID, face ID, index ranges.

```jldoctest
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

```jldoctest
using MeshArrays
γ=GridSpec("LatLonCap",MeshArrays.GRID_LLC90)
d=γ.read(γ.path*"Depth.data",MeshArray(γ,γ.ioPrec))
τ=Tiles(γ,30,30)
td=Tiles(τ,d)

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
    GridAddWS!(Γ::Dict)

Compute XW, YW, XS, and YS (vector field locations) from XC, YC (tracer
field locations) and add them to Γ.
"""
function GridAddWS!(Γ::Dict)

    XC=exchange(Γ["XC"])
    YC=exchange(Γ["YC"])
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
