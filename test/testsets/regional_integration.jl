@testset "Regional Integration:" begin
    G,M,files=Integration.example()
    @suppress show(M)

    allones=1.0 .+0*G.hFacC
    vol0=sum(G.RAC*G.DRF*G.hFacC)

    vol=Integration.volumes(M,G)
    @test isapprox(sum(vol),vol0)

    G,M,files=Integration.example(option=:streamlined_loop)
    vol=[b(allones) for b in M.h_sum]
    @test isapprox(sum(vol),vol0)

    rgns=Integration.define_regions(option=:basins,grid=G)
    rgns=Integration.define_regions(option=:dlat_10,grid=G)
    rgns=Integration.define_regions(option=(30,10),grid=G)
    @test isa(rgns,NamedTuple)

    G,M,files=Integration.example()
    files=fill("?",3)
    rd0(F,var,tim,tmp)=tim*ones(tmp)
    H=Integration.loops(M,files=files,rd=rd0)
    @test isa(H,Array)

    G,M,files=Integration.example(option=:streamlined_loop)
    files=fill("?",3)
    H=Integration.streamlined_loop(M,files=files,rd=rd0)
    @test isa(H,Array)
end
