
## demo functions:

"""
    demo1(gridChoice::String)

Demonstrate basic functionalities (load grid, arithmetic, exchange, gradient,
etc.). Call sequence:

```
include("examples/Demos.jl")
(D,Dexch,Darr,DD)=demo1("LatLonCap",MeshArrays.GRID_LLC90);
```
"""
function demo1(gridChoice::String,GridParentDir="./")

    γ=GridSpec(gridChoice,GridParentDir)

    D=γ.read(γ.path*"Depth.data",MeshArray(γ,γ.ioPrec))

    1000 .+ D
    D .+ 1000
    D+D
    D .- 1000
    1000 .- D
    D-D
    1000*D
    D*1000
    D*D
    D/1000
    1000 ./ D
    D/D

    Dexch=exchange(D,4)
    Darr=γ.write(D)
    DD=γ.read(Darr,D)
    DD .== D

    Γ=GridLoad(γ,option="full")

    (dFLDdx, dFLDdy)=gradient(Γ.YC,Γ)
    (dFLDdxEx,dFLDdyEx)=exchange(dFLDdx,dFLDdy,4)

    H=Γ.hFacC
    Ha=γ.write(H)
    Hb=γ.read(Ha,H)

    println(typeof(H))
    println(size(H))
    println(H.fSize[1])
    println(size(H.f[1]))
    if ndims(H)==3&&size(H,3)>40
        println(typeof(view(H,:,:,40)))
    elseif ndims(H.f)==2&&size(H.f,2)>40
        println(typeof(view(H,:,40)))
    end

    return (D,Dexch,Darr,DD)

end

##

"""
    demo2()

Demonstrate higher level functions using `smooth`. Call sequence:

```
include("examples/Demos.jl")
(Rini,Rend,DXCsm,DYCsm)=demo2()

using CairoMakie, MeshArrays
heatmap(Rend,title="smoothed noise",colorrange=(-0.5,0.5))
heatmap(Rini,title="raw noise",colorrange=(-0.5,0.5))
```

"""
function demo2()

    #Pre-requisite: either load predefined grid using `demo1` or call `GridOfOnes`
    isdir(MeshArrays.GRID_LLC90) ? Γ=GridLoad(ID=:LLC90,option="full") : (γ,Γ)=Grids_simple.GridOfOnes("CubeSphere",6,100)

    (Rini,Rend,DXCsm,DYCsm)=demo2(Γ)
end

function demo2(Γ::NamedTuple)

    γ=Γ.XC.grid

    #initialize 2D field of random numbers
    tmp1=randn(Float32,Tuple(γ.ioSize))
    Rini=γ.read(tmp1,MeshArray(γ,Float32))

    #apply land mask
    if ndims(Γ.hFacC.f[1])>2
        tmp1=mask(view(Γ.hFacC,:,:,1),NaN,0)
    elseif ndims(Γ.hFacC.f)>1
        #tmp1=mask(view(Γ.hFacC,:,1),NaN,0)
        tmp1=similar(Rini)
        for i=1:length(tmp1.fIndex); tmp1[i]=Γ.hFacC[i,1]; end;
        tmp1=MeshArrays.mask(tmp1,NaN,0)
    else
        tmp1=MeshArrays.mask(Γ.hFacC,NaN,0)
    end
    msk=fill(1.,tmp1) + 0. *tmp1;
    Rini=msk*Rini;

    #specify smoothing length scales in x, y directions
    DXCsm=3*Γ.DXC; DYCsm=3*Γ.DYC;

    #apply smoother
    Rend=MeshArrays.smooth(Rini,DXCsm,DYCsm,Γ);

    return (Rini,Rend,DXCsm,DYCsm)

end

"""
    demo3()

Demonstrate ocean transport computations. Call sequence:

```
include("examples/Demos.jl")
de=demo3()

using CairoMakie, MeshArrays
lines!(Axis(Figure()[1,1],title="meridional transport"),de.lat,de.Tr/1e6)
fig1=current_figure()
fig2=heatmap(1e-6*de.TrspX,title="U comp. in Sv",colorrange=(-10,10))
fig3=heatmap(1e-6*de.TrspY,title="V comp. in Sv",colorrange=(-10,10))
```
"""
function demo3()

    γ=GridSpec(ID=:LLC90)
    Γ=GridLoad(γ,option=:light)

    TrspX=γ.read(γ.path*"TrspX.bin",MeshArray(γ,Float32))
    TrspY=γ.read(γ.path*"TrspY.bin",MeshArray(γ,Float32))
    TauX=γ.read(γ.path*"TauX.bin",MeshArray(γ,Float32))
    TauY=γ.read(γ.path*"TauY.bin",MeshArray(γ,Float32))
    SSH=γ.read(γ.path*"SSH.bin",MeshArray(γ,Float32))

    (UV, LC, Tr)=demo3(TrspX,TrspY,Γ)

    (TrspX=TrspX,TrspY=TrspY,TauX=TauX,TauY=TauY,SSH=SSH,γ=γ,Γ=Γ,LC=LC,Tr=Tr,lat=-89.0:89.0,UV=UV)
end

function demo3(U::MeshArray,V::MeshArray,Γ::NamedTuple)
    LC=LatitudeCircles(-89.0:89.0,Γ)
    #UV=Dict("U"=>U,"V"=>V,"dimensions"=>["x","y","z","t"],"factors"=>["dxory","dz"])
    UV=Dict("U"=>U,"V"=>V,"dimensions"=>["x","y"])
    Tr=Array{Float64,1}(undef,length(LC));
    for i=1:length(LC)
       Tr[i]=ThroughFlow(UV,LC[i],Γ)
    end
    return UV, LC, Tr
end
