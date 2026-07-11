@testset "NEMO_GRID:" begin
    lst=NEMO_GRID.variable_NTA()
    nam=NEMO_GRID.variable_in_NEMO(:XC,lst)
    @test isa(nam,Symbol)

    grid_data=Dict(:glamt=>zeros(1442,1021))
    XC_a=NEMO_GRID.convert_one_grid_variable(grid_data,lst[14])

    grid=(XC=XC_a,YC=XC_a,RAC=XC_a,DXG=XC_a,DYG=XC_a)
    G=NEMO_GRID.grid_to_MeshArrays(grid)

    XC_e=NEMO_GRID.exchange(G.XC)
    @test isa(XC_e,MeshArray_wh)

    G=NEMO_GRID.add_angle_CS_SN(G)
    @test haskey(G,:AngleCS)

    z=[1:10]; G=Dict()
    grid_data=(gdept_0=z,gdepw_0=z,e3t_0=z,e3w_0=z)
    NEMO_GRID.add_one_dim_variables!(G,grid_data)
    @test isa(G,Dict)
end
