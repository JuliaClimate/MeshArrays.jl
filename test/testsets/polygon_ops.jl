@testset "polygon operations" begin
    fil=MeshArrays.Dataset("countries_shp1",do_read=false)
    pol=MeshArrays.read_shp(fil,format=:Shapefile)
    name,rule=MeshArrays.within_pol(pol; ID=11)
    rule_vec = (x,y) -> rule.(x,y)

    np=10000; lo=-180 .+360*rand(np); la=-90 .+180*rand(np);
    np_in=sum(rule_vec(lo,la))
    @test np_in>0

    path_MITgcm=MITgcm.getdata("mitgcmsmallverif")
    path_grid=joinpath(path_MITgcm,"MITgcm","verification","tutorial_held_suarez_cs","input")
    pols,pols3D=MeshArrays.Polygons.polygons_demo(path_grid)
    Depth=GridLoadVar("Depth",GridSpec(ID=:CS32))
    @test isa(pols[1][1,1],GI.Polygon)

    fig=MeshArrays.plot_examples(:polygons_plot,pols,color=Depth)
    MeshArrays.plot_examples(:polygons_plot_dev1,pols,pols3D,sphere_view=true)
    MeshArrays.plot_examples(:polygons_plot_dev1,pols,pols3D,sphere_view=false)
    @test isa(fig,Makie.FigureAxisPlot)
end
