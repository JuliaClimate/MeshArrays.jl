using Test
using MeshArrays

p=dirname(pathof(MeshArrays))
include(joinpath(p,"../examples/Demos.jl"))

@testset "MeshArrays tests:" begin
for nTopo=1:3
    if nTopo==1; grTopo="CubeSphere"; nFaces=6; N=200;
    elseif nTopo==2; grTopo="LatLonCap"; nFaces=5; N=200;
    elseif nTopo==3; grTopo="PeriodicChanel"; nFaces=1; N=500;
    end;
    Npt=nFaces*N*N
    GridVariables=GridOfOnes(grTopo,nFaces,N)
    mygrid=GridVariables["XC"].grid
    @test mygrid.class == grTopo
    Rini= 0.; Rend= 0.;
    (Rini,Rend,DXCsm,DYCsm)=demo2(GridVariables);
    @test isa(Rend,MeshArray)
    @test sum(isfinite.(Rend)) == Npt
    Sini=sqrt(sum(Rini*Rini)/(Npt-1.0))
    Send=sqrt(sum(Rend*Rend)/(Npt-1.0))
    #println([Sini Send])
    @test isapprox(Sini,1.000; atol=1e-2)
    @test isapprox(Send,0.093; atol=1e-2)
end;
end;
