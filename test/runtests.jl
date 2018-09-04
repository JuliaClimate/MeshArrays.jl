using Test
using GCMFaces

@testset "GCMFaces tests:" begin
    N=200
    Npt=6*N*N
    @test GCMGridOnes("cs",6,N) == "GCMGridOnes: passed"
    Rini= 0.; Rend= 0.;
    (Rini,Rend,DXCsm,DYCsm)=demo2();
    @test isa(Rini,gcmfaces)
    Sini=sqrt(sum(Rini*Rini)/(Npt-1.0))
    Send=sqrt(sum(Rend*Rend)/(Npt-1.0))
    #println([Sini Send])
    @test isapprox(Sini,1.000; atol=1e-2)
    @test isapprox(Send,0.093; atol=1e-2)
end;

