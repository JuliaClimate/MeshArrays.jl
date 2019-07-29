
## demo functions:

"""
demo1(gridChoice)

Demonstrate basic fucntions (arithmetic, exchange, GCMGridLoad, gradient, etc.). Example call:

```
isdir("GRID_LLC90") ? (D,Dexch,Darr,DD)=demo1("LLC90") : nothing
```
"""
function demo1(gridChoice)

    GCMGridSpec(gridChoice)

    D=read_bin(MeshArrays.grDir*"Depth.data",MeshArrays.ioPrec)

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
    DD=convert2array(Darr)

    GCMGridLoad()

    (dFLDdx, dFLDdy)=gradient(MeshArrays.YC)
    (dFLDdxEx,dFLDdyEx)=exchange(dFLDdx,dFLDdy,4)

    view(MeshArrays.hFacC,:,:,40)
    #show(fsize(MeshArrays.hFacC,1))
    #show(fsize(view(MeshArrays.hFacC,:,:,40),1))

    return (D,Dexch,Darr,DD)

end

##

"""
demo2()

Demonstrate higher level functions using smooth() and

```
isdir("GRID_LLC90") ? demo1("LLC90") : GCMGridOnes("cs",6,100)
(Rini,Rend,DXCsm,DYCsm)=demo2()
@time Rend=smooth(Rini,DXCsm,DYCsm)

include(joinpath(dirname(pathof(MeshArrays)),"gcmfaces_plot.jl"))
qwckplot(Rini)
qwckplot(Rend)
```

"""
function demo2()

    #Pre-requisite: either load predefined grid using `demo1` or call `GCMGridOnes`

    #initialize 2D field of random numbers
    tmp1=convert2gcmfaces(MeshArrays.XC);
    tmp1=randn(Float32,size(tmp1));
    Rini=convert2gcmfaces(tmp1);

    #apply land mask
    if ndims(MeshArrays.hFacC.f[1])>2
        tmp1=mask(view(MeshArrays.hFacC,:,:,1),NaN,0);
    else
        tmp1=mask(MeshArrays.hFacC,NaN,0);
    end
    msk=fill(1.,tmp1) + 0. *tmp1;
    Rini=msk*Rini;

    #specify smoothing length scales in x, y directions
    DXCsm=3*MeshArrays.DXC; DYCsm=3*MeshArrays.DYC;

    #apply smoother
    Rend=smooth(Rini,DXCsm,DYCsm);

    return (Rini,Rend,DXCsm,DYCsm)

end

"""
demo3()

Demonstrate computations of ocean meridional transports. Calling sequence:

```
!isdir("GRID_LLC90")||!isdir("nctiles_climatology") ? error("missing files") : nothing

GCMGridSpec("LLC90")
GCMGridLoad()

include(joinpath(dirname(pathof(MeshArrays)),"gcmfaces_nctiles.jl"))
fileName="nctiles_climatology/UVELMASS/UVELMASS"
U=read_nctiles(fileName,"UVELMASS");
fileName="nctiles_climatology/VVELMASS/VVELMASS"
V=read_nctiles(fileName,"VVELMASS");

(UV, LC, Tr)=demo3(U,V);

using Statistics
include(joinpath(dirname(pathof(MeshArrays)),"gcmfaces_plot.jl"))
qwckplot(UV["U"][:,:,1,1],"U component (note varying face orientations)")
qwckplot(UV["V"][:,:,1,1],"V component (note varying face orientations)")
plot(dropdims(mean(sum(Tr,dims=2),dims=3),dims=(2,3))/1e6,title="meridional transport")
```
"""
function demo3(U,V)

    LC=LatCircles(-89.0:89.0)

    U=mask(U,0.0)
    V=mask(V,0.0)

    UV=Dict("U"=>U,"V"=>V,"dimensions"=>["x","y","z","t"],"factors"=>["dxory","dz"])

    n=size(U)
    Tr=Array{Float64}(undef,length(LC),n[3],n[4])
    for i=1:length(LC)
        Tr[i,:,:]=TransportThrough(UV,LC[i])
    end

    return UV, LC, Tr

end
