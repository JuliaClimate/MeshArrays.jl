@testset "Plotting:" begin
    γ=GridSpec(ID=:LLC90)
    Γ=GridLoad(γ;option="light")
    D=Γ.Depth
    λ=interpolation_setup(Γ=Γ,
        lon=[i for i=-170.:20.0:170., j=-80.:20.0:80.],
        lat=[j for i=-170.:20.0:170., j=-80.:20.0:80.])
    λ=interpolation_setup()

    lines(pol_json); lines!(pol_json)
    plot(pol_json); plot!(pol_json)

    basins=demo.ocean_basins()
    AtlExt=demo.extended_basin(basins,:Atl)
    sections,path_sec=demo.ocean_sections(Γ)
    my_section=demo.one_section(Γ,[127 127],[-25 -68])

    fig=MeshArrays.plot_examples(:smoothing_demo,D,D)
    (fig1,fig2,fig3)=MeshArrays.plot_examples(:interpolation_demo,Γ)

    fake_ov=40e6*cosd.(360*(1:179)./100)*exp.(-0.1*(-20:29).^2)'
    MeshArrays.plot_examples(:meriodional_overturning,Γ,fake_ov)
    MeshArrays.plot_examples(:northward_transport,rand(179))

    MeshArrays.plot_examples(:gradient_EN,λ,D,D)
    MeshArrays.plot_examples(:gradient_xy,λ,D,D)

    ## more methods

    scatter(Γ.XC,Γ.YC,color=:black)
    heatmap(D,interpolation=λ)
    scatter!(current_axis(),Γ.XC,Γ.YC,color=:red)

    heatmap(D) #will display tile by tile
    heatmap(D,interpolation=λ,title="ocean depth") #same but w title

    lon0=-160
    proj=Proj.Transformation(MA_preset=2,lon0=lon0)
    Dint=reshape(Interpolate(D,λ.f,λ.i,λ.j,λ.w),size(λ.lon))

    Interpolate(D,λ)
    InterpolationFactors(Γ,30.0,30.0)

    ###

    MeshArraysMakieExt = Base.get_extension(MeshArrays, :MeshArraysMakieExt)
    pol=MeshArraysMakieExt.pol_to_Makie(pol_shp)
    dest="+proj=eqearth +lon_0=$(lon0) +lat_1=0.0 +x_0=0.0 +y_0=0.0 +ellps=GRS80"
    MeshArraysMakieExt.split(pol,dest)
    MeshArraysMakieExt.split(pol,Observable(dest))
    MeshArraysMakieExt.split(Observable(pol),Observable(dest))
    MeshArraysMakieExt.split(Observable(pol),dest)

    ###

    f = Figure()
    ax = f[1, 1] = Axis(f, aspect = DataAspect(), title = "Ocean Depth (m)")
	pr_ax=MeshArrays.ProjAxis(ax; proj=proj,lon0=lon0)
    for a in [surface! contourf! contour!]
        surf = a(pr_ax,λ.lon,λ.lat,0*λ.lat; color=Dint,
			colorrange=(0.0,6000.0), colormap=:berlin, shading = NoShading)
    end
	lines!(pr_ax; polygons=pol_shp,color=:black,linewidth=0.5)
	MeshArrays.grid_lines!(pr_ax;color=:lightgreen,linewidth=0.5)
	f

	meta=(colorrange=(0.0,6000.0),cmap=:BrBG_10,ttl="Ocean Depth (m)",lon0=lon0)
	data=(lon=λ.lon,lat=λ.lat,var=Dint,meta=meta) #,polygons=pol_shp)
    plot_examples(:projmap,data,lon0,proj)
    plot_examples(:simple_heatmap,data)

    MeshArraysMakieExt.heatmap_globalmap(D)
    MeshArraysMakieExt.heatmap_interpolation(D,λ)
    MeshArraysMakieExt.heatmap_xy(D,1:10,1:10)

end
