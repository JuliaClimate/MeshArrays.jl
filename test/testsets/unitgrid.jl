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

    # Test read(xx::Array, x::AbstractMeshArray) — untested pathway
    # Create a 4D array with shape (ioSize..., n3, n4) for PeriodicDomain grid
    xx_data = randn(10, 10, 1, 1)  # 10×10 ioSize, no extra dims
    A_template = MeshArray(γ, Float64)
    B = read(xx_data, A_template)
    @test isa(B, MeshArray)
    @test eltype(B) == Float64
    @test size(B) == size(A_template)

    # Test read(xx::Array, x::AbstractMeshArray) with n3 > 1
    xx_3d = randn(10, 10, 2, 1)
    C_template = MeshArray(γ, Float64, 2)
    C = read(xx_3d, C_template)
    @test isa(C, MeshArray)
    @test size(C) == size(C_template)
    @test size(C, 2) == 2

    # Test read(fil::String, x::AbstractMeshArray) with allocate=true
    tmp_file = tempname()
    write(tmp_file, Γ.XC)
    A_template2 = MeshArray(γ, Float64)
    D = read(tmp_file, A_template2)
    @test isa(D, MeshArray)
    @test size(D) == size(A_template2)
    @test eltype(D) == Float64

    # Test similar(x, allocate=true) pathway
    A_1d = MeshArray(γ)
    A_1d_full = similar(A_1d; allocate=true)
    @test isa(A_1d_full, MeshArray)
    @test size(A_1d_full) == size(A_1d)
    # Verify faces are allocated (not 0×0)
    @test size(A_1d_full.f[1]) == size(A_1d.f[1])

    A_2d = MeshArray(γ, Float64, 3)
    A_2d_full = similar(A_2d; allocate=true)
    @test isa(A_2d_full, MeshArray)
    @test size(A_2d_full) == size(A_2d)
    @test size(A_2d_full.f[1]) == size(A_2d.f[1])

    A_3d = MeshArray(γ, Float64, 3, 4)
    A_3d_full = similar(A_3d; allocate=true)
    @test isa(A_3d_full, MeshArray)
    @test size(A_3d_full) == size(A_3d)
    @test size(A_3d_full.f[1]) == size(A_3d.f[1])
end
