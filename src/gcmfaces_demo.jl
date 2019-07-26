
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
include(joinpath(dirname(pathof(MeshArrays)),"gcmfaces_plot.jl"))
include(joinpath(dirname(pathof(MeshArrays)),"gcmfaces_nctiles.jl"))
!isdir("GRID_LLC90")||!isdir("nctiles_climatology") ? error("missing files") : nothing

demo3()
```
"""
function demo3()

    GCMGridSpec("LLC90")
    GCMGridLoad()

    LatValues=-179.5:179.5
    LatCircles=Array{Dict}(undef,length(LatValues))

    for j=1:length(LatValues)
        mskCint=1*(MeshArrays.YC .>= LatValues[j])
        mskC=similar(mskCint)
        mskW=similar(mskCint)
        mskS=similar(mskCint)

        mskCint=exchange(mskCint,1)

        for i=1:mskCint.nFaces
            tmp1=mskCint.f[i]
            # tracer mask:
            tmp2=tmp1[2:end-1,1:end-2]+tmp1[2:end-1,3:end]+
            tmp1[1:end-2,2:end-1]+tmp1[3:end,2:end-1]
            mskC.f[i]=1((tmp2.>0).&(tmp1[2:end-1,2:end-1].==0))
            # velocity masks:
            mskW.f[i]=tmp1[2:end-1,2:end-1] - tmp1[1:end-2,2:end-1]
            mskS.f[i]=tmp1[2:end-1,2:end-1] - tmp1[2:end-1,1:end-2]
        end

        LatCircles[j]=Dict("lat"=>LatValues[j],"mskC"=>mskC,"mskW"=>mskW,"mskS"=>mskS)
    end

    fileName="nctiles_climatology/UVELMASS/UVELMASS"
    U=read_nctiles(fileName,"UVELMASS")
    fileName="nctiles_climatology/VVELMASS/VVELMASS"
    V=read_nctiles(fileName,"VVELMASS")

    return LatCircles

end
