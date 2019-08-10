
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

    D=mygrid.read(mygrid.path*"Depth.data",gcmfaces(mygrid))

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
    Darr=write(D)
    DD=read(Darr,D)

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
qwckplot(Rini,"raw noise")
qwckplot(Rend,"smoothed noise")
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
    tmp1=randn(Float32,Tuple(mygrid.ioSize))
    Rini=read(tmp1,gcmfaces(mygrid,Float32))

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
(UV,LC,Tr)=MeshArrays.demo3();
using Plots; plot(Tr/1e6,title="meridional transport")

include(joinpath(dirname(pathof(MeshArrays)),"gcmfaces_plot.jl"))
qwckplot(UV["U"][:,:,1,1],"U component (note varying face orientations)")
qwckplot(UV["V"][:,:,1,1],"V component (note varying face orientations)")
```
"""
function demo3()

    mygrid=GCMGridSpec("LLC90")
    GridVariables=GCMGridLoad(mygrid)

    TrspX=mygrid.read(mygrid.path*"TrspX.bin",gcmfaces(mygrid,Float32))
    TrspY=mygrid.read(mygrid.path*"TrspY.bin",gcmfaces(mygrid,Float32))
    TauX=mygrid.read(mygrid.path*"TauX.bin",gcmfaces(mygrid,Float32))
    TauY=mygrid.read(mygrid.path*"TauY.bin",gcmfaces(mygrid,Float32))
    SSH=mygrid.read(mygrid.path*"SSH.bin",gcmfaces(mygrid,Float32))

    (UV, LC, Tr)=demo3(TrspX,TrspY,GridVariables)

end

function demo3(U::gcmfaces,V::gcmfaces,GridVariables::Dict)

    LC=LatitudeCircles(-89.0:89.0,GridVariables)

    #UV=Dict("U"=>U,"V"=>V,"dimensions"=>["x","y","z","t"],"factors"=>["dxory","dz"])
    UV=Dict("U"=>U,"V"=>V,"dimensions"=>["x","y"])

    Tr=Array{Float64,1}(undef,length(LC));
    for i=1:length(LC)
       Tr[i]=ThroughFlow(UV,LC[i],GridVariables)
    end

    return UV, LC, Tr

end
