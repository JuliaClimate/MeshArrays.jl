
## demo functions:

"""
    demo1(gridChoice::String)

Demonstrate basic functionalities (load grid, arithmetic, exchange, gradient,
etc.). Call sequence:

```
!isdir("GRID_LLC90") ? error("missing files") : nothing

(D,Dexch,Darr,DD)=demo1("LatLonCap","GRID_LLC90/");
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

    Γ=GridLoad(γ)

    (dFLDdx, dFLDdy)=gradient(Γ["YC"],Γ)
    (dFLDdxEx,dFLDdyEx)=exchange(dFLDdx,dFLDdy,4)

    H=Γ["hFacC"]
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
(Rini,Rend,DXCsm,DYCsm)=demo2();

using Plots
include(joinpath(dirname(pathof(MeshArrays)),"../examples/Plots.jl"))
heatmap(Rend,title="smoothed noise",clims=(-0.5,0.5))
heatmap(Rini,title="raw noise",clims=(-0.5,0.5))

#gr(); contour(Rend,clims=(-0.5,0.5),levels=-0.5:0.1:0.5,fill=true)
```

"""
function demo2()

    #Pre-requisite: either load predefined grid using `demo1` or call `GridOfOnes`
    isdir("GRID_LLC90") ? Γ=GridLoad(GridSpec("LatLonCap","GRID_LLC90/")) : (γ,Γ)=GridOfOnes("CubeSphere",6,100)

    (Rini,Rend,DXCsm,DYCsm)=demo2(Γ)
end

function demo2(Γ::Dict)

    γ=Γ["XC"].grid

    #initialize 2D field of random numbers
    tmp1=randn(Float32,Tuple(γ.ioSize))
    Rini=γ.read(tmp1,MeshArray(γ,Float32))

    #apply land mask
    if ndims(Γ["hFacC"].f[1])>2
        tmp1=mask(view(Γ["hFacC"],:,:,1),NaN,0)
    elseif ndims(Γ["hFacC"].f)>1
        #tmp1=mask(view(Γ["hFacC"],:,1),NaN,0)
        tmp1=similar(Rini)
        for i=1:length(tmp1.fIndex); tmp1[i]=Γ["hFacC"][i,1]; end;
        tmp1=mask(tmp1,NaN,0)
    else
        tmp1=mask(Γ["hFacC"],NaN,0)
    end
    msk=fill(1.,tmp1) + 0. *tmp1;
    Rini=msk*Rini;

    #specify smoothing length scales in x, y directions
    DXCsm=3*Γ["DXC"]; DYCsm=3*Γ["DYC"];

    #apply smoother
    Rend=smooth(Rini,DXCsm,DYCsm,Γ);

    return (Rini,Rend,DXCsm,DYCsm)

end

"""
    demo3()

Demonstrate ocean transport computations. Call sequence:

```
(UV,LC,Tr)=demo3();
using Plots; plot(Tr/1e6,title="meridional transport")

using Plots
include(joinpath(dirname(pathof(MeshArrays)),"../examples/Plots.jl"))
heatmap(1e-6*UV["U"],title="U comp. in Sv",clims=(-10,10))
heatmap(1e-6*UV["V"],title="V comp. in Sv",clims=(-10,10))
```
"""
function demo3()

    γ=GridSpec("LatLonCap","GRID_LLC90/")
    Γ=GridLoad(γ)

    TrspX=γ.read(γ.path*"TrspX.bin",MeshArray(γ,Float32))
    TrspY=γ.read(γ.path*"TrspY.bin",MeshArray(γ,Float32))
    TauX=γ.read(γ.path*"TauX.bin",MeshArray(γ,Float32))
    TauY=γ.read(γ.path*"TauY.bin",MeshArray(γ,Float32))
    SSH=γ.read(γ.path*"SSH.bin",MeshArray(γ,Float32))

    (UV, LC, Tr)=demo3(TrspX,TrspY,Γ)

end

function demo3(U::MeshArray,V::MeshArray,Γ::Dict)

    LC=LatitudeCircles(-89.0:89.0,Γ)

    #UV=Dict("U"=>U,"V"=>V,"dimensions"=>["x","y","z","t"],"factors"=>["dxory","dz"])
    UV=Dict("U"=>U,"V"=>V,"dimensions"=>["x","y"])

    Tr=Array{Float64,1}(undef,length(LC));
    for i=1:length(LC)
       Tr[i]=ThroughFlow(UV,LC[i],Γ)
    end

    return UV, LC, Tr

end
