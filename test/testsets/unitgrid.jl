@testset "UnitGrid:" begin
    C=MeshArray(randn(20,10))
    D=MeshArray(randn(20,10,3))

    (Γ,γ)=Grids_simple.UnitGrid( (80,90) , (20,30) ; option="full")
    @test isa(γ,gcmgrid)

    γ=Grids_simple.GridSpec_ones("PeriodicDomain",nP=10)
    Γ=Grids_simple.GridLoad_ones(γ,option="full")
    @test isa(Γ,NamedTuple)

    #various read/write functions
    read(write(Γ.XC),γ)
    tmp=tempname()
    write(tmp,Γ.XC)
    MeshArrays.read_tiles(tmp,Γ.XC)
    MeshArrays.write_tiles(Γ.XC)
    MeshArrays.write_tiles(tmp,Γ.XC)
    @test isfile(tmp)

    xy=Grids_simple.xy_OISST()
    xy=Grids_simple.xy_Oscar()

    xy=Grids_simple.xy_IAP()
    gr=Grids_simple.grid_factors(xy)
    dep=[10 100 1000]; msk=ones(gr[:XC].fSize[1]...,3)
    gr=Grids_simple.grid_add_z(gr,dep,msk)

    @test haskey(gr,:hFacC)
end
