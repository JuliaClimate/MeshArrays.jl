using Test, Documenter
using MeshArrays

p=dirname(pathof(MeshArrays))
include(joinpath(p,"../examples/Demos.jl"))

@testset "MeshArrays tests:" begin
for nTopo=1:3
    if nTopo==1; grTopo="CubeSphere"; nFaces=6; N=200;
    elseif nTopo==2; grTopo="LatLonCap"; nFaces=5; N=200;
    elseif nTopo==3; grTopo="PeriodicChannel"; nFaces=1; N=500;
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
end;
end;

@testset "Helmholtz decomposition test:" begin

    function trsp_read(myspec::String,mypath::String)
        γ=GridSpec(myspec,mypath)
        TrspX=γ.read(mypath*"TrspX.bin",MeshArray(γ,Float32))
        TrspY=γ.read(mypath*"TrspY.bin",MeshArray(γ,Float32))
        TauX=γ.read(mypath*"TauX.bin",MeshArray(γ,Float32))
        TauY=γ.read(mypath*"TauY.bin",MeshArray(γ,Float32))
        SSH=γ.read(mypath*"SSH.bin",MeshArray(γ,Float32))
        return TrspX, TrspY, TauX, TauY, SSH
    end
    
Γ=GridLoad(GridSpec("LatLonCap",MeshArrays.GRID_LLC90))
μ =Γ["hFacC"][:,1]; μ[findall(μ.>0.0)].=1.0; μ[findall(μ.==0.0)].=NaN
(Tx,Ty,τx,τy,η)=trsp_read("LatLonCap",MeshArrays.GRID_LLC90);

TrspCon=μ.*convergence(Tx,Ty)
#scalar potential
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
end;

@testset "doctests" begin
    doctest(MeshArrays; manual = false)
end
