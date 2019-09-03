using Test
using MeshArrays

@testset "MeshArrays tests:" begin
for nTopo=1:3
    if nTopo==1; grTopo="cs"; nFaces=6; N=200;
    elseif nTopo==2; grTopo="llc"; nFaces=5; N=200;
    elseif nTopo==3; grTopo="ll"; nFaces=1; N=500;
    end;
    Npt=nFaces*N*N
    GridVariables=GCMGridOnes(grTopo,nFaces,N)
    mygrid=GridVariables["XC"].grid
    @test mygrid.class == grTopo
    Rini= 0.; Rend= 0.;
    (Rini,Rend,DXCsm,DYCsm)=MeshArrays.demo2(GridVariables);
    @test isa(Rend,MeshArray)
    @test sum(isfinite(Rend)) == Npt
    Sini=sqrt(sum(Rini*Rini)/(Npt-1.0))
    Send=sqrt(sum(Rend*Rend)/(Npt-1.0))
    #println([Sini Send])
    @test isapprox(Sini,1.000; atol=1e-2)
    @test isapprox(Send,0.093; atol=1e-2)
end;
end;
