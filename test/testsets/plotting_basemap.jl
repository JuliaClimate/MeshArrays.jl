@testset "plotting" begin
    lon,lat,earth_img=demo.get_basemap()
    plot_examples(:basemap,lon,lat,earth_img)

    lon0=-160
    proj=Proj.Transformation(MA_preset=2,lon0=lon0)
    plot_examples(:baseproj,proj,lon0,pol=pol_shp)
end
