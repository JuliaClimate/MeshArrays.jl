using Test, Documenter
using MeshArrays

p=dirname(pathof(MeshArrays))
include(joinpath(p,"../examples/Demos.jl"))

@testset "MeshArrays tests:" begin
    for nTopo=1:4
        if nTopo==1; grTopo="CubeSphere"; nFaces=6; N=200;
        elseif nTopo==2; grTopo="LatLonCap"; nFaces=5; N=200;
        elseif nTopo==3; grTopo="PeriodicChannel"; nFaces=1; N=400;
        elseif nTopo==4; grTopo="PeriodicDomain"; nFaces=1; N=400;
        end;
        Npt=nFaces*N*N
        γ,Γ=GridOfOnes(grTopo,nFaces,N)
        @test γ.class == grTopo
        Rini= 0.; Rend= 0.;
        (Rini,Rend,DXCsm,DYCsm)=demo2(Γ);
        @test isa(Rend,MeshArray)
        @test sum(isfinite.(Rend)) == Npt
        Sini=sqrt(sum(Rini*Rini)/(Npt-1.0))
        Send=sqrt(sum(Rend*Rend)/(Npt-1.0))
        #println([Sini Send])
        @test isapprox(Sini,1.000; atol=1e-2)
        @test isapprox(Send,0.093; atol=1e-2)
        (dRdx,dRdy)=gradient(Rend,Γ)
        exchange(Rend,2)
        exchange(dRdx,dRdy,1)        
        (dRdx_e,dRdy_e)=exchange(dRdx,dRdy,1)        
        @test isa(dRdx_e,MeshArray)
    end
end

@testset "Transport computations:" begin
    #Load grid and transport / vector field
    γ=GridSpec("LatLonCap",MeshArrays.GRID_LLC90)
    Tx=γ.read(MeshArrays.GRID_LLC90*"TrspX.bin",MeshArray(γ,Float32))
    Ty=γ.read(MeshArrays.GRID_LLC90*"TrspY.bin",MeshArray(γ,Float32))
    Γ=GridLoad(γ); μ =Γ["hFacC"][:,1]
    μ[findall(μ.>0.0)].=1.0; μ[findall(μ.==0.0)].=NaN

    #Meridional transport integral
    uv=Dict("U"=>Tx,"V"=>Ty,"dimensions"=>["x","y"])
    L=-85.0:5.0:85.0; LC=LatitudeCircles(L,Γ)
    T=Array{Float64,1}(undef,length(LC))
    [T[i]=1e-6*ThroughFlow(uv,LC[i],Γ) for i=1:length(LC)]

    #See: OceanTransports/helper_functions.jl
    #u,v,uC,vC=rotate_uv(uv,Γ);

    #Transport convergence   
    TrspCon=μ.*convergence(Tx,Ty)
    #Scalar potential
    TrspPot=ScalarPotential(TrspCon)
    #Divergent transport component
    (TxD,TyD)=gradient(TrspPot,Γ)
    TxD=TxD.*Γ["DXC"]
    TyD=TyD.*Γ["DYC"]
    #Rotational transport component
    TxR = Tx-TxD
    TyR = Ty-TyD
    #vector Potential
    TrspPsi=VectorPotential(TxR,TyR,Γ);

    tst=extrema(write(TrspPsi))

    @test prod(isapprox.(tst,(-8.08f7,2.19f8)))
end

@testset "gcmfaces type:" begin
    γ=GridSpec("CubeSphere",MeshArrays.GRID_CS32)

    MeshArrays.gcmfaces(γ)
    MeshArrays.gcmfaces(γ,Float32)
    tmp=MeshArrays.gcmfaces(γ,Float32,3)
    tmp[3,1,2]=1.0
    view(tmp,:,1,1:2)
    MeshArrays.gcmfaces(γ,tmp.f)
    MeshArrays.fsize(tmp)
    MeshArrays.fsize(tmp,2)
    MeshArrays.fsize(tmp.f)
    MeshArrays.fsize(tmp.f,2)
    size(tmp)
    size(tmp,3)
    show(tmp)

    x=tmp[1:10,1,1:2]; y=x[2]; x[3]=1.0
    view(x,1:3,:,1)
    MeshArrays.gcmsubset(γ,x.f,x.fSize,x.aSize,x.i,x.iSize)
    MeshArrays.fsize(x.f)
    MeshArrays.fsize(x.f,2)
    size(x)
    size(x,3)

    MeshArrays.fijind(tmp,10)
    MeshArrays.fsize(tmp)
    @test isa(tmp,MeshArrays.gcmfaces)

    tmp=MeshArray(γ)
    tmp1=findall(tmp.>0)
    tmp[tmp1]
    tmp[tmp1].=1.0
    size(tmp1)
    tmp1[2]
    view(tmp1,:)
    show(tmp1)
    similar(tmp1)
    #tmp[tmp1]

    show(tmp)
    MeshArrays.getindexetc(tmp,2)
    MeshArray(γ,tmp.f,meta=tmp.meta)
    MeshArray(γ,meta=tmp.meta)
    MeshArray(γ,Float32,3,4)
    MeshArray(γ,Float32,tmp.fSize,tmp.fIndex,2,3)

    tmp1=MeshArray(γ,Float32,3)
    MeshArray(γ,tmp1.f,meta=tmp1.meta)
    MeshArrays.getindexetc(tmp1,2,1)
end

@testset "doctests" begin
    doctest(MeshArrays; manual = false)
end
