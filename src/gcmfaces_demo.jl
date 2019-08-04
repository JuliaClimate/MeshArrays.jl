
## demo functions:

"""
    demo1(gridChoice::String)

Demonstrate basic functionalities (load grid, arithmetic, exchange, gradient,
etc.). Call sequence:

```
!isdir("GRID_LLC90") ? error("missing files") : nothing

(D,Dexch,Darr,DD)=MeshArrays.demo1("LLC90");
```
"""
function demo1(gridChoice::String)

    mygrid=GCMGridSpec(gridChoice)

    D=read_bin(mygrid.path*"Depth.data",mygrid.ioPrec,mygrid)

    1000+D
    D+1000
    D+D
    D-1000
    1000-D
    D-D
    1000*D
    D*1000
    D*D
    D/1000
    1000/D
    D/D

    Dexch=exchange(D,4)
    Darr=convert2array(D)
    DD=convert2array(Darr,mygrid)

    GridVariables=GCMGridLoad(mygrid)

    (dFLDdx, dFLDdy)=gradient(GridVariables["YC"],GridVariables)
    (dFLDdxEx,dFLDdyEx)=exchange(dFLDdx,dFLDdy,4)

    view(GridVariables["hFacC"],:,:,40)
    #show(fsize(GridVariables["hFacC"],1))
    #show(fsize(view(GridVariables["hFacC"],:,:,40),1))

    return (D,Dexch,Darr,DD)

end

##

"""
    demo2()

Demonstrate higher level functions using `smooth`. Call sequence:

```
(Rini,Rend,DXCsm,DYCsm)=MeshArrays.demo2();

include(joinpath(dirname(pathof(MeshArrays)),"gcmfaces_plot.jl"))
qwckplot(Rini)
qwckplot(Rend)
```

"""
function demo2()

    #Pre-requisite: either load predefined grid using `demo1` or call `GCMGridOnes`
    isdir("GRID_LLC90") ? GridVariables=GCMGridLoad(GCMGridSpec("LLC90")) : GridVariables=GCMGridOnes("cs",6,100)

    (Rini,Rend,DXCsm,DYCsm)=demo2(GridVariables)
end

function demo2(GridVariables::Dict)

    mygrid=GridVariables["XC"].grid

    #initialize 2D field of random numbers
    tmp1=convert2gcmfaces(GridVariables["XC"])
    tmp1=randn(Float32,size(tmp1))
    Rini=convert2gcmfaces(tmp1,mygrid)

    #apply land mask
    if ndims(GridVariables["hFacC"].f[1])>2
        tmp1=mask(view(GridVariables["hFacC"],:,:,1),NaN,0)
    else
        tmp1=mask(GridVariables["hFacC"],NaN,0)
    end
    msk=fill(1.,tmp1) + 0. *tmp1;
    Rini=msk*Rini;

    #specify smoothing length scales in x, y directions
    DXCsm=3*GridVariables["DXC"]; DYCsm=3*GridVariables["DYC"];

    #apply smoother
    Rend=smooth(Rini,DXCsm,DYCsm,GridVariables);

    return (Rini,Rend,DXCsm,DYCsm)

end

"""
    demo3()

Demonstrate ocean transport computations. Call sequence:

```
!isdir("GRID_LLC90")||!isdir("nctiles_climatology") ? error("missing files") : nothing
include(joinpath(dirname(pathof(MeshArrays)),"gcmfaces_nctiles.jl"))

(UV,LC,Tr)=MeshArrays.demo3();

using Statistics, Plots
plot(dropdims(mean(sum(Tr,dims=2),dims=3),dims=(2,3))/1e6,title="meridional transport")

include(joinpath(dirname(pathof(MeshArrays)),"gcmfaces_plot.jl"))
qwckplot(UV["U"][:,:,1,1],"U component (note varying face orientations)")
qwckplot(UV["V"][:,:,1,1],"V component (note varying face orientations)")
```
"""
function demo3()

    mygrid=GCMGridSpec("LLC90")
    GridVariables=GCMGridLoad(mygrid)

    fileName="nctiles_climatology/UVELMASS/UVELMASS"
    U=Main.read_nctiles(fileName,"UVELMASS",mygrid);
    fileName="nctiles_climatology/VVELMASS/VVELMASS"
    V=Main.read_nctiles(fileName,"VVELMASS",mygrid);

    (UV, LC, Tr)=demo3(U,V,GridVariables)

end

function demo3(U::gcmfaces,V::gcmfaces,GridVariables::Dict)

    LC=LatCircles(-89.0:89.0,GridVariables)

    U=mask(U,0.0)
    V=mask(V,0.0)

    UV=Dict("U"=>U,"V"=>V,"dimensions"=>["x","y","z","t"],"factors"=>["dxory","dz"])

    n=size(U)
    Tr=Array{Float64}(undef,length(LC),n[3],n[4])
    for i=1:length(LC)
        Tr[i,:,:]=TransportThrough(UV,LC[i],GridVariables)
    end

    return UV, LC, Tr

end
