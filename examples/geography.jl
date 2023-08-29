### A Pluto.jl notebook ###
# v0.19.27

using Markdown
using InteractiveUtils

# This Pluto notebook uses @bind for interactivity. When running this notebook outside of Pluto, the following 'mock version' of @bind gives bound variables a default value (instead of an error).
macro bind(def, element)
    quote
        local iv = try Base.loaded_modules[Base.PkgId(Base.UUID("6e696c72-6542-2067-7265-42206c756150"), "AbstractPlutoDingetjes")].Bonds.initial_value catch; b -> missing; end
        local el = $(esc(element))
        global $(esc(def)) = Core.applicable(Base.get, el) ? Base.get(el) : iv(el)
        el
    end
end

# ╔═╡ 0c1448a9-52e7-41c9-b19a-e548d92db948
using Pkg; Pkg.activate(".")

# ╔═╡ d123161e-49f1-11ec-1c1b-51871624545d
begin
	using MeshArrays, OceanStateEstimation
	using PlutoUI, Statistics, JLD2
	using CairoMakie, Proj
	using GeoJSON, Shapefile
	using Downloads, ZipFile

	"Done with packages"
end

# ╔═╡ 39788eec-dd8b-4d65-a0f9-ea73a2d6690c
TableOfContents()

# ╔═╡ cc3e9d0c-8b71-432b-b68d-a8b832ca5f26
md"""# Geography and Visualization

This tutorial focuses on the relation between the gridded space in [MeshArrays.jl](https://juliaclimate.github.io/MeshArrays.jl/dev/) and the geospatial space (longitude, lattitude). A common use case is for plotting in the form of maps. Going back and forth between the two is also often useful for analysis of simulations (transports, averages, etc) and e.g. simulating particle trajectories.
"""

# ╔═╡ ca8a8f1b-a225-46c4-93c1-ce1a4b016461
md"""## Interpolation

Here we interpolate from the global grid (see [MeshArrays.jl](https://github.com/JuliaClimate/MeshArrays.jl)) to a set of arbitary locations (longitude, latitude pairs) as is commonly done, for example, to plot gridded fields in geographic coordinates (example below shows Ocean bottom depth) or to compare climate models to sparse field observations (e.g., [Argo profiles](https://github.com/JuliaOcean/ArgoData.jl)).

Later in this notebook, we break down the `MeshArrays.Interpolate` function into several steps. In brief, the program finds a grid point quadrilateral (4 grid points) that encloses the target location. Computation is chuncked in subdomains (tiles) to allow for parallelism. `MeshArrays.InterpolationFactors` outputs interpolation coefficients, which can then be reused repeatedly -- e.g. to speed up calls to `MeshArrays.Interpolate` in loops that generate animations.
"""

# ╔═╡ 85aac4c5-bf04-4f9a-a233-a3f24231762e
begin
	source_bind = @bind source Select([:shp_example,:json_example])
	list=["wintri","natearth2","longlat"]
	proj_bind = @bind proj Select(list)
	lon0_bind = @bind lon0 Select(collect(-160:40:160),default=-160)
	md"""## Projection
	
	See `simple_heatmap` function (in appendices), which uses `Proj.jl` and `Makie.jl`.
	
	$(proj_bind)
	$(lon0_bind)
	$(source_bind)
	"""
end

# ╔═╡ 6aab5feb-06b6-4fbd-8732-83b70f397f00
begin
	#push!(lonPairs,[-68 -63]); push!(latPairs,[-54 -66]); push!(namPairs,"Drake Passage");

	lon1_slct = @bind lon1 Select(-179.0:180.0;default=20.0)
	lon2_slct = @bind lon2 Select(-179.0:180.0;default=-63.0)
	lat1_slct = @bind lat1 Select(-89.0:89.0;default=-20.0)
	lat2_slct = @bind lat2 Select(-89.0:89.0;default=-66.0)
	md"""## Grid Line Paths

The shortest path between two locations on the surface of a sphere is the one that follows a [great circle](https://en.wikipedia.org/wiki/Great-circle_distance). The method used in `MeshArrays.jl` identifies the path along grid cell edges that corresponds closely to the great circle path.

`MeshArrays.jl` provides the `Transect` function to this end, which relies on the `MeshArrays.edge_mask` function. `MeshArrays.edge_mask` can also be used find paths that approximate circles of constant latitude, or delineate a subdomain for example.

| Location 1         | Location 2         | 
|:------------------:|:------------------:|
| lon = $(lon1_slct) | lon = $(lon2_slct) | 
| lat = $(lat1_slct) | lat = $(lat2_slct) | 
	"""
end

# ╔═╡ c81ad3bb-a924-4393-bc73-0607cfebf75f
begin
	faceID_slct = @bind faceID Select(1:5;default=1)
	
	md"""## Grid Cell Areas

	When the goal is to compute averages or integrals over a subdomain or e.g. the whole Ocean, it is important to account for the corresponding land mask and grid cell areas. Here we visualize variations in grid cell areas for a specific global grid. Here the grid resolution is increased towards the Equator and at high latitudes. Details will differ for other grids, but they should not be overlooked in computing e.g. global means or zonal averages.
	
	- select a grid subdomain : $(faceID_slct)
	"""
end

# ╔═╡ f5d7639e-f524-4134-9188-7bc8586c6883
md"""## More"""

# ╔═╡ 8d2e6d5f-7d89-4a27-82cb-0a587914d717
md"""#### Interpolation Scheme


The Interpolation Code? Let's Break It Down:

- find nearest neighbor (`MeshArray` & `set`)
- define subdomain tiles (`MeshArray` -> `tiles`)
- exchange and start loop (`tile` & `subset`)
    - local stereographic projection
    - define array of quadrilaterals
    - find enclosing quadrilaterals
    - compute interpolation coefficients
"""

# ╔═╡ bf537ab9-b0b0-43f8-a6ad-61798e7b5016
md"""#### Simple Polygon Map"""

# ╔═╡ 25144e1b-21fc-4cc9-b63d-7b26eab1a673
md"""## Appendices

#### Data
"""

# ╔═╡ 39924391-38ce-46a1-877f-80a7975340a0
begin
	pth=MeshArrays.GRID_LLC90
	γ=GridSpec("LatLonCap",pth)
	Γ=GridLoad(γ;option="full")	
	#LC=LatitudeCircles(-89.0:89.0,Γ)
	basins=MeshArrays.ocean_basins()
	"Done with grid"
end

# ╔═╡ be38ff51-3526-44a0-9d8c-9209355e4a4a
begin
	file_int=MeshArrays.interpolation_setup()
	λ=MeshArrays.interpolation_setup(file_int)
	μ=MeshArrays.land_mask(Γ)
	
	Depth_interpolated=Interpolate(λ.μ*Γ.Depth,λ.f,λ.i,λ.j,λ.w)
	Depth_interpolated=reshape(Depth_interpolated,size(λ.lon))
	#Depth_interpolated[findall(Depth_interpolated.==0.0)].=NaN
	
	"Done with interpolating Γ.Depth"
end

# ╔═╡ 303b2f54-a446-4dbe-a9a4-e6f424003722
begin
	nn=30

	mm=Int(round(size(Γ.XC[3],1)/nn))
	oo=Int(floor(mm*mm*13/10))
	#tileID=0*Γ.YC; [tileID[t.face][t.i,t.j].=t.tile for t in τ];	
	tileID_select = @bind ii NumberField(1:mm*mm*13,default=oo)

	τ=Tiles(γ,nn,nn)
	XC=Tiles(τ,Γ.XC)
	YC=Tiles(τ,Γ.YC)
	Depth_tiled=Tiles(τ,Γ.Depth)

	md"""
	!!! note 
	    The rest of this notebook focuses on using fields defined on native ocean model grids --without interpolation to a regular grid. This approach generally allows for improved accuracy calculations in model space.
	
	## Domain Splitting

	Here we decompose the global domain into $(mm*mm) tiles of $(nn) by $(nn).
	
	- tile size : $(nn)
	- tile ID : $(tileID_select)
	"""
end

# ╔═╡ d9a53bc6-19e2-48d9-b9c3-f76c3a533197
MeshArrays.examples_plot(:tiled_example,λ,Depth_interpolated,XC,YC,Depth_tiled,ii)

# ╔═╡ 42562cd8-04c5-4aef-87dc-5d0f87a9d204
begin
		lons=[lon1 lon2]
		lats=[lat1 lat2]
		my_section=MeshArrays.one_section(Γ,lons,lats)
end

# ╔═╡ 794d822e-2101-41ba-8780-fac95fc02075
begin
	fig1=heatmap(λ.μ*Γ.Depth,λ)
	MeshArrays.examples_plot(:one_section,fig1,lons,lats,my_section)
	fig1
end

# ╔═╡ b0d576fc-971a-47c7-9a57-f2c788083bcd
begin

	basin_slct = @bind basin_nam Select(basins.name;default=basins.name[1])

	#basinID_slct = @bind basinID Select(1:maximum(basins);default=1)
	
	md"""## Regional Masks

	There are many ways to split up the global domain into Oceans, Seas, etc. Here is one example.
	
	- select a basin nam : $(basin_slct)
	"""
end

# ╔═╡ 2afc6934-9e07-46bf-a4f1-df33ce5a6fe9
MeshArrays.examples_plot(:ocean_basins,Γ,λ,basins,basin_nam)

# ╔═╡ a8defb50-fce4-4a0b-ac33-deb95f0b826b
MeshArrays.examples_plot(:cell_area,Γ,λ,1)

# ╔═╡ 6b72d272-eefc-45f2-9442-ef38057e4f09
(fig01,fig2,fig3)=MeshArrays.examples_plot(:interpolation_demo,Γ)

# ╔═╡ 6cc62cf0-cb30-4d93-aad6-2ab16f60f95f
begin
	sec,path_sec=MeshArrays.ocean_sections(Γ)
	sec
end

# ╔═╡ 897a49b2-9763-4020-a476-5e0fccda1cfb
begin
	(Γ.XW,Γ.YW)
	LC=LatitudeCircles(79.0,Γ)
	aa=LC[1].tabC
	locClat=( lon=[Γ.XC[aa[i,1]][aa[i,2],aa[i,3]] for i in 1:size(aa,1)],
				lat=[Γ.YC[aa[i,1]][aa[i,2],aa[i,3]] for i in 1:size(aa,1)])

	"Done with latitude line"
end

# ╔═╡ 229d395f-f3b4-40df-a482-264d540be875
begin
	if source==:shp_example
		fil=MeshArrays.download_polygons("ne_110m_admin_0_countries.shp")
	elseif source==:json_example
		fil=MeshArrays.download_polygons("countries.geojson")
	end
	l=MeshArrays.read_polygons(fil)
	"Done With Reading Country Polygons"
end

# ╔═╡ 54f6cb2f-025e-40e4-97a8-0f8ad4b35278
begin
	f1=Figure()
	ax1=Axis(f1[1,1])
	[lines!(ax1,l1,color = :black, linewidth = 0.5) for l1 in l]
	f1
end

# ╔═╡ d448431e-13a9-440c-9463-9174d7400cf1
begin
	proj_ID=findall(list.==proj)[1]
	trans=Proj.Transformation(MA_preset=proj_ID,lon0=lon0)
	ll=Makie.coordinates.(split(l,lon0))
end

# ╔═╡ 41960267-fff9-4bc4-a7bc-aceea2217c63
begin
	meta=(colorrange=(0.0,6000.0),cmap=:BrBG_10,ttl="Ocean Depth (m)",lon0=lon0)
	data=(lon=λ.lon,lat=λ.lat,var=Depth_interpolated,meta=meta,polygons=ll)
	summary(data)
end

# ╔═╡ 9c30b8df-b95d-40e1-b7e3-2dacbdda49ed
#simple_heatmap(data)
MeshArrays.examples_plot(:projmap,data,trans)

# ╔═╡ 7fc02644-b037-425a-b660-cc6904a95037
md"""#### Code"""

# ╔═╡ Cell order:
# ╟─39788eec-dd8b-4d65-a0f9-ea73a2d6690c
# ╟─cc3e9d0c-8b71-432b-b68d-a8b832ca5f26
# ╟─ca8a8f1b-a225-46c4-93c1-ce1a4b016461
# ╠═be38ff51-3526-44a0-9d8c-9209355e4a4a
# ╟─85aac4c5-bf04-4f9a-a233-a3f24231762e
# ╟─9c30b8df-b95d-40e1-b7e3-2dacbdda49ed
# ╟─303b2f54-a446-4dbe-a9a4-e6f424003722
# ╟─d9a53bc6-19e2-48d9-b9c3-f76c3a533197
# ╟─6aab5feb-06b6-4fbd-8732-83b70f397f00
# ╟─794d822e-2101-41ba-8780-fac95fc02075
# ╟─42562cd8-04c5-4aef-87dc-5d0f87a9d204
# ╟─b0d576fc-971a-47c7-9a57-f2c788083bcd
# ╟─2afc6934-9e07-46bf-a4f1-df33ce5a6fe9
# ╟─c81ad3bb-a924-4393-bc73-0607cfebf75f
# ╟─a8defb50-fce4-4a0b-ac33-deb95f0b826b
# ╟─f5d7639e-f524-4134-9188-7bc8586c6883
# ╟─8d2e6d5f-7d89-4a27-82cb-0a587914d717
# ╟─6b72d272-eefc-45f2-9442-ef38057e4f09
# ╟─bf537ab9-b0b0-43f8-a6ad-61798e7b5016
# ╟─54f6cb2f-025e-40e4-97a8-0f8ad4b35278
# ╟─25144e1b-21fc-4cc9-b63d-7b26eab1a673
# ╠═41960267-fff9-4bc4-a7bc-aceea2217c63
# ╠═39924391-38ce-46a1-877f-80a7975340a0
# ╠═6cc62cf0-cb30-4d93-aad6-2ab16f60f95f
# ╠═897a49b2-9763-4020-a476-5e0fccda1cfb
# ╠═229d395f-f3b4-40df-a482-264d540be875
# ╠═d448431e-13a9-440c-9463-9174d7400cf1
# ╟─7fc02644-b037-425a-b660-cc6904a95037
# ╠═0c1448a9-52e7-41c9-b19a-e548d92db948
# ╠═d123161e-49f1-11ec-1c1b-51871624545d
