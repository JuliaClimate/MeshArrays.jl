### A Pluto.jl notebook ###
# v0.19.22

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

# ╔═╡ d123161e-49f1-11ec-1c1b-51871624545d
begin
	using Pkg; Pkg.activate()
	using MeshArrays, OceanStateEstimation
	using PlutoUI, Statistics, JLD2
	using CairoMakie
	using Proj
	"Done with packages"
end

# ╔═╡ 5f1085c7-f78a-433f-853e-0e3505fd99f9
module projections
	using CairoMakie, Proj, MeshArrays
	fil=joinpath(dirname(pathof(MeshArrays)),"..","examples","projections.jl")
	include(fil)
end

# ╔═╡ 0c1448a9-52e7-41c9-b19a-e548d92db948
module polygons
	using Downloads, GeoJSON, GeoInterface, Shapefile
	using GeometryBasics, Observables, MeshArrays
	fil=joinpath(dirname(pathof(MeshArrays)),"..","examples","polygons.jl")
	include(fil)
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
	lon0_bind = @bind lon0 Select(collect(-160:40:160),default=80)
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

# ╔═╡ 6b72d272-eefc-45f2-9442-ef38057e4f09
#(fig1,fig2,fig3)=interp_example()

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
	
	basins=read(joinpath(MeshArrays.GRID_LLC90,"v4_basin.bin"),MeshArray(γ,Float32))
	basin_list=["Pacific","Atlantic","indian","Arctic","Bering Sea","South China Sea","Gulf of Mexico",
		"Okhotsk Sea","Hudson Bay","Mediterranean Sea","Java Sea","North Sea","Japan Sea",
		"Timor Sea","East China Sea","Red Sea","Gulf","Baffin Bay","GIN Seas","Barents Sea"]

	"Done with grid"
end

# ╔═╡ 303b2f54-a446-4dbe-a9a4-e6f424003722
begin
	nn=10

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

# ╔═╡ b0d576fc-971a-47c7-9a57-f2c788083bcd
begin

	basin_slct = @bind basin_nam Select(basin_list;default=basin_list[1])

	#basinID_slct = @bind basinID Select(1:maximum(basins);default=1)
	
	md"""## Regional Masks

	There are many ways to split up the global domain into Oceans, Seas, etc. Here is one example.
	
	- select a basin nam : $(basin_slct)
	"""
end

# ╔═╡ e76c02c3-bae6-4110-b26c-3f6b7547453e
function interp_example()
	lon=collect(0.1:0.5:2.1); lat=collect(0.1:0.5:2.1);
	#lon=[66.75]; lat=[-69.75];

	(f,i,j,w,_,_,_)=InterpolationFactors(Γ,vec(lon),vec(lat))
	Depth_int=Interpolate(Γ.Depth,f,i,j,w)

	#find nearest neighbor (`MeshArray` & `set`)
	(f,i,j,c)=knn(Γ.XC,Γ.YC,lon,lat)
	[write(Γ.XC)[c] write(Γ.YC)[c]]

	#define subdomain tiles (`MeshArray` -> `tiles`)
	ni=30; nj=30;
	τ=Tiles(γ,ni,nj)

	tiles=MeshArray(γ,Int);
	[tiles[τ[i].face][τ[i].i,τ[i].j].=i for i in 1:length(τ)]

	#

	XCtiles=Tiles(τ,exchange(Γ.XC))
	YCtiles=Tiles(τ,exchange(Γ.YC))

	iiTile=tiles[f[1]][i[1],j[1]]; iiFace=τ[iiTile].face

	ii0=minimum(τ[iiTile].i)+Int(ni/2); jj0=minimum(τ[iiTile].j)+Int(nj/2)
	XC0=Γ.XG.f[iiFace][ii0,jj0]; YC0=Γ.YG.f[iiFace][ii0,jj0]
	#XC0=66.5000; YC0=-64.4201

	#
	
	fig1 = Figure(resolution = (900,600), backgroundcolor = :grey95)
	ax1 = Axis(fig1[1,1],xlabel="longitude",ylabel="latitude")

	scatter!(ax1,XCtiles[iiTile][:],YCtiles[iiTile][:],marker=:+,c=:blue)
	scatter!(ax1,[XC0],[YC0],c=:red)
	scatter!(ax1,lon[:],lat[:],c=:green)

	#Local Stereographic Projection
	(x_grid,y_grid)=StereographicProjection(XC0,YC0,XCtiles[iiTile],YCtiles[iiTile])
	(x_trgt,y_trgt)=StereographicProjection(XC0,YC0,lon,lat)
	~isa(x_trgt,Array) ? x_trgt=[x_trgt] : nothing
	~isa(y_trgt,Array) ? y_trgt=[y_trgt] : nothing

	#Identify Enclosing Quadrilaterals
	(x_quad,y_quad,i_quad,j_quad)=MeshArrays.QuadArrays(x_grid,y_grid)
	angsum=zeros(size(x_quad,1),size(x_trgt,1))
	MeshArrays.PolygonAngle(x_quad,y_quad,x_trgt,y_trgt,angsum)
	ii=findall(angsum.>180.)
	ii=[ii[j].I[1] for j in 1:length(ii)]

	#Interpolation Coefficients
	px=permutedims(x_quad[ii[1],:]); py=permutedims(y_quad[ii[1],:])
	ox=x_trgt[1]; oy=y_trgt[1]; ow=zeros(4)
	MeshArrays.QuadCoeffs(px,py,ox,oy,ow);

	#

	fig2 = Figure(resolution = (900,600), backgroundcolor = :grey95)
	ax2 = Axis(fig2[1,1],xlabel="x",ylabel="y")

	scatter!(ax2,x_grid[:],y_grid[:],marker=:+,c=:blue)
	scatter!(ax2,[0.],[0.],c=:red)
	scatter!(ax2,x_quad[ii,:][:],y_quad[ii,:][:],c=:orange)
	scatter!(ax2,x_trgt,y_trgt,c=:green)

	#

	(f,i,j,w)=InterpolationFactors(Γ,lon,lat)

	lon_a=NaN*similar(lon)
	lat_a=NaN*similar(lat)
	for jj=1:length(lon)
	    if !isnan(sum(w[jj,:]))
	        x=[Γ.XC[f[jj,ii]][i[jj,ii],j[jj,ii]] for ii=1:4]
	        y=[Γ.YC[f[jj,ii]][i[jj,ii],j[jj,ii]] for ii=1:4]
	        lon_a[jj]=sum(w[jj,:].*x)
	        lat_a[jj]=sum(w[jj,:].*y)
	    end
	end
	
	#or equivalently:
	lon_b=Interpolate(Γ.XC,f,i,j,w)
	lat_b=Interpolate(Γ.YC,f,i,j,w)

	#
	
	fig3 = Figure(resolution = (900,600), backgroundcolor = :grey95)
	ax3 = Axis(fig3[1,1],xlabel="longitude",ylabel="latitude")

	scatter!(ax3,XCtiles[iiTile][:],YCtiles[iiTile][:],marker=:+,c=:blue)
	scatter!(ax3,[XC0],[YC0],color=:red,marker=:diamond,markersize=24.0)
	scatter!(ax3,lon,lat,color=:orange,markersize=24.0)
	scatter!(ax3,lon_a,lat_a,color=:black,marker=:star4,markersize=24.0)
	#scatter!(ax3,lon_b,lat_b,color=:black,marker=:star4,markersize=24.0)
	#scatter!(ax3,lon_c,lat_c,color=:black,marker=:star4,markersize=24.0)
	
	(fig1,fig2,fig3)
end

# ╔═╡ 963e421c-43fb-43d3-b667-1b9912f940b8
begin
	lons=[lon1 lon2]
	lats=[lat1 lat2]
	name="An Example Grid Path"

	x0,y0,z0,R=MeshArrays.rotate_points(lons,lats)
	x,y,z=MeshArrays.rotate_XCYC(Γ,R)
	mskCint=1.0*(z.>0)
	mskCedge,mskWedge,mskSedge=MeshArrays.edge_mask(mskCint)
	
	#for plotting
	tabC=MeshArrays.MskToTab(mskCedge)
	locCwhole=( lon=[Γ.XC[tabC[i,1]][tabC[i,2],tabC[i,3]] for i in 1:size(tabC,1)],
				lat=[Γ.YC[tabC[i,1]][tabC[i,2],tabC[i,3]] for i in 1:size(tabC,1)])

	mskCedge,mskWedge,mskSedge=MeshArrays.shorter_paths!((x,y,z),(x0,y0,z0),(mskCedge,mskWedge,mskSedge))

	#for plotting
	tabC=MeshArrays.MskToTab(mskCedge)
	locCshort=( lon=[Γ.XC[tabC[i,1]][tabC[i,2],tabC[i,3]] for i in 1:size(tabC,1)],
				lat=[Γ.YC[tabC[i,1]][tabC[i,2],tabC[i,3]] for i in 1:size(tabC,1)])
	
	"Done with transport line masks"
end

# ╔═╡ 2507c335-2886-42b5-b6b8-63279a2d60fe
begin
	lonPairs=[]    
	latPairs=[]    
	namPairs=[]    
	
	push!(lonPairs,[-173 -164]); push!(latPairs,[65.5 65.5]); push!(namPairs,"Bering Strait");
	push!(lonPairs,[-5 -5]); push!(latPairs,[34 40]); push!(namPairs,"Gibraltar");
	push!(lonPairs,[-81 -77]); push!(latPairs,[28 26]); push!(namPairs,"Florida Strait");
	push!(lonPairs,[-81 -79]); push!(latPairs,[28 22]); push!(namPairs,"Florida Strait W1");
	push!(lonPairs,[-76 -76]); push!(latPairs,[21 8]); push!(namPairs,"Florida Strait S1");
	push!(lonPairs,[-77 -77]); push!(latPairs,[26 24]); push!(namPairs,"Florida Strait E1");
	push!(lonPairs,[-77 -77]); push!(latPairs,[24 22]); push!(namPairs,"Florida Strait E2");
	push!(lonPairs,[-65 -50]); push!(latPairs,[66 66]); push!(namPairs,"Davis Strait");
	push!(lonPairs,[-35 -20]); push!(latPairs,[67 65]); push!(namPairs,"Denmark Strait");
	push!(lonPairs,[-16 -7]); push!(latPairs,[65 62.5]); push!(namPairs,"Iceland Faroe");
	push!(lonPairs,[-6.5 -4]); push!(latPairs,[62.5 57]); push!(namPairs,"Faroe Scotland");
	push!(lonPairs,[-4 8]); push!(latPairs,[57 62]); push!(namPairs,"Scotland Norway");
	push!(lonPairs,[-68 -63]); push!(latPairs,[-54 -66]); push!(namPairs,"Drake Passage");
	push!(lonPairs,[103 103]); push!(latPairs,[4 -1]); push!(namPairs,"Indonesia W1");
	push!(lonPairs,[104 109]); push!(latPairs,[-3 -8]); push!(namPairs,"Indonesia W2");
	push!(lonPairs,[113 118]); push!(latPairs,[-8.5 -8.5]); push!(namPairs,"Indonesia W3");
	push!(lonPairs,[118 127 ]); push!(latPairs,[-8.5 -15]); push!(namPairs,"Indonesia W4");
	push!(lonPairs,[127 127]); push!(latPairs,[-25 -68]); push!(namPairs,"Australia Antarctica");
	push!(lonPairs,[38 46]); push!(latPairs,[-10 -22]); push!(namPairs,"Madagascar Channel");
	push!(lonPairs,[46 46]); push!(latPairs,[-22 -69]); push!(namPairs,"Madagascar Antarctica");
	push!(lonPairs,[20 20]); push!(latPairs,[-30 -69.5]); push!(namPairs,"South Africa Antarctica");
	push!(lonPairs,[-76 -72]); push!(latPairs,[21 18.5]); push!(namPairs,"Florida Strait E3");
	push!(lonPairs,[-72 -72]); push!(latPairs,[18.5 10]); push!(namPairs,"Florida Strait E4");
end

# ╔═╡ 6cc62cf0-cb30-4d93-aad6-2ab16f60f95f
if true
	pth_trsp=joinpath(tempdir(),"ECCO_transport_lines")
	!isdir(pth_trsp) ? mkdir(pth_trsp) : nothing
	
	for ii in 1:length(lonPairs)
		lons=Float64.(lonPairs[ii])
		lats=Float64.(latPairs[ii])
		name=namPairs[ii]
		Trsct=Transect(name,lons,lats,Γ)
		fil_Trsct=joinpath(tempdir(),"ECCO_transport_lines","$(Trsct.name).jld2")
		!isfile(fil_Trsct) ? jldsave(fil_Trsct,tabC=Trsct.tabC,tabW=Trsct.tabW,tabS=Trsct.tabS) : nothing
	end

	readdir(pth_trsp)
else
	">> switch false to true inside this cell to compute all transects"
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
		l=polygons.PolygonReading.process_shp()
	elseif source==:json_example
		l=polygons.PolygonReading.process_json()
	end
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
	ll=polygons.LineSplitting.regroup([polygons.LineSplitting.split(LineString(a),lon0) for a in l])
end

# ╔═╡ 7fc02644-b037-425a-b660-cc6904a95037
md"""#### Code"""

# ╔═╡ 75c98143-dd36-4446-adfa-c440fdf8aab1
function setup_interp(Γ)
	μ =Γ.hFacC[:,1]
	μ[findall(μ.>0.0)].=1.0
	μ[findall(μ.==0.0)].=NaN

	fil=joinpath(ScratchSpaces.ECCO,"interp_coeffs_halfdeg.jld2")
	if !isfile(fil)
		OceanStateEstimation.ECCOdiags_add("release2")
		OceanStateEstimation.ECCOdiags_add("interp_coeffs")
	end

	λ = load(fil)
	λ = MeshArrays.Dict_to_NamedTuple(λ)
end

# ╔═╡ be38ff51-3526-44a0-9d8c-9209355e4a4a
begin
	λ=setup_interp(Γ)
	DD=Interpolate(λ.μ*Γ.Depth,λ.f,λ.i,λ.j,λ.w)
	DD=reshape(DD,size(λ.lon))
	#DD[findall(DD.==0.0)].=NaN
	"Done with interpolating Γ.Depth"
end

# ╔═╡ d9a53bc6-19e2-48d9-b9c3-f76c3a533197
let
	fig = Figure(resolution = (900,600), backgroundcolor = :grey95,colormap=:thermal)
	ax = Axis(fig[1,1],xlabel="longitude",ylabel="latitude",title="grid cell area (log10 of m^2)")
	hm1=heatmap!(ax,λ.lon[:,1],λ.lat[1,:],DD,colormap=:grayC,colorrange=(0.0,4000.0))
	sc1=scatter!(ax,XC[ii][:],YC[ii][:],color=Depth_tiled[ii][:],
		markersize=4.0,colormap=:thermal,colorrange=(0.0,6000.0))
	fig
end

# ╔═╡ ed36a2a5-44ea-43a7-a3bd-13f234b6580d
let	
	fig_path = Figure(resolution = (900,600), backgroundcolor = :grey95)
	ax = Axis(fig_path[1,1],xlabel="longitude",ylabel="latitude",title="Ocean Depth Shown in Colors (interpolated Γ.Depth)")
	hm1=heatmap!(ax,λ.lon[:,1],λ.lat[1,:],DD,colormap=:spring,colorrange=(0.0,6000.0))

	scatter!(ax,locCwhole.lon[:],locCwhole.lat[:],color=:blue,markersize=2.0)
	scatter!(ax,locCshort.lon[:],locCshort.lat[:],color=:black,markersize=4.0)
	scatter!(ax,lons[:],lats[:],color=:blue)

	#scatter!(ax,locClat.lon[:],locClat.lat[:],color=:cyan,markersize=4.0)

	Colorbar(fig_path[1,2], hm1, height = Relative(0.65))

	#save("fig_path.png",fig_path)
	fig_path
end

# ╔═╡ 31756b2e-20df-47c9-aaa8-5583e6a81267
let
	fig = Figure(resolution = (900,600), backgroundcolor = :grey95)
	ax = Axis(fig[1,1],xlabel="longitude",ylabel="latitude",title=basin_nam* " (shown in red)")

	basinID=findall(basin_list.==basin_nam)[1]
	
	mx=maximum(basins)
	for ff in 1:length(Γ.RAC)
		col=λ.μ[ff][:].*(basins[ff][:].==basinID)
		kk=findall(col.>0.0)
		!isempty(kk) ? scatter!(ax,Γ.XC[ff][kk],Γ.YC[ff][kk],color=:red,markersize=2.0) : nothing
		kk=findall((col.==0.0).*(!isnan).(λ.μ[ff][:]))
		!isempty(kk) ? scatter!(ax,Γ.XC[ff][kk],Γ.YC[ff][kk],color=basins[ff][kk],
			colorrange=(0.0,mx),markersize=2.0,colormap=:lisbon) : nothing
	end
	Colorbar(fig[1,2], colormap=:lisbon, colorrange=(0.0, mx), height = Relative(0.65))

	fig
end

# ╔═╡ 1477cd5d-7ee3-4af8-95cc-13309db00520
begin
	fig = Figure(resolution = (900,600), backgroundcolor = :grey95,colormap=:thermal)
	ax = Axis(fig[1,1],xlabel="longitude",ylabel="latitude",title="grid cell area (log10 of m^2)")

	rng=(0.0, maximum(Γ.RAC))
	rng=(8.8,10.2)
	for ff in 1:length(Γ.RAC)
		col=log10.(λ.μ[ff][:].*Γ.RAC[ff][:])
		kk=findall((!isnan).(col))
		if ff==faceID
			scatter!(ax,Γ.XC[ff][kk],Γ.YC[ff][kk],color=col[kk],
				colorrange = rng,markersize=2.0,colormap=:thermal)
		else
			scatter!(ax,Γ.XC[ff][kk],Γ.YC[ff][kk],color=:gray,markersize=2.0)
		end
	end
	Colorbar(fig[1,2], colorrange=rng, height = Relative(0.65))

	fig
end

# ╔═╡ 41960267-fff9-4bc4-a7bc-aceea2217c63
begin
	meta=(colorrange=(0.0,6000.0),cmap=:BrBG_10,ttl="Ocean Depth (m)",lon0=lon0)
	#data=(lon=λ.lon,lat=λ.lat,var=DD,meta=meta,polygons=l)
	data=(lon=λ.lon,lat=λ.lat,var=DD,meta=meta,polygons=coordinates.(ll))
	summary(data)
end

# ╔═╡ 9c30b8df-b95d-40e1-b7e3-2dacbdda49ed
#simple_heatmap(data)
projections.ProjMakie.projmap(data,proj_ID)

# ╔═╡ 0b188f96-87b1-41e4-ba66-1fa4d5252bd8
function simple_heatmap(dat)	
	lons = dat.lon[:,1]
	lats = dat.lat[1,:]
	field = dat.var

	fig = Figure(resolution = (1200,800), fontsize = 22)
	ax = Axis(fig[1,1])
	hm1 = heatmap!(ax, lons, lats, field, colorrange=dat.meta.colorrange, colormap=dat.meta.cmap)
	Colorbar(fig[1,2], hm1, height = Relative(0.65))

	fig
end

# ╔═╡ Cell order:
# ╟─39788eec-dd8b-4d65-a0f9-ea73a2d6690c
# ╟─cc3e9d0c-8b71-432b-b68d-a8b832ca5f26
# ╟─ca8a8f1b-a225-46c4-93c1-ce1a4b016461
# ╟─be38ff51-3526-44a0-9d8c-9209355e4a4a
# ╟─85aac4c5-bf04-4f9a-a233-a3f24231762e
# ╟─9c30b8df-b95d-40e1-b7e3-2dacbdda49ed
# ╟─303b2f54-a446-4dbe-a9a4-e6f424003722
# ╟─d9a53bc6-19e2-48d9-b9c3-f76c3a533197
# ╟─6aab5feb-06b6-4fbd-8732-83b70f397f00
# ╟─ed36a2a5-44ea-43a7-a3bd-13f234b6580d
# ╟─b0d576fc-971a-47c7-9a57-f2c788083bcd
# ╟─31756b2e-20df-47c9-aaa8-5583e6a81267
# ╟─c81ad3bb-a924-4393-bc73-0607cfebf75f
# ╟─1477cd5d-7ee3-4af8-95cc-13309db00520
# ╟─f5d7639e-f524-4134-9188-7bc8586c6883
# ╟─8d2e6d5f-7d89-4a27-82cb-0a587914d717
# ╟─e76c02c3-bae6-4110-b26c-3f6b7547453e
# ╠═6b72d272-eefc-45f2-9442-ef38057e4f09
# ╟─bf537ab9-b0b0-43f8-a6ad-61798e7b5016
# ╟─54f6cb2f-025e-40e4-97a8-0f8ad4b35278
# ╟─25144e1b-21fc-4cc9-b63d-7b26eab1a673
# ╠═41960267-fff9-4bc4-a7bc-aceea2217c63
# ╟─39924391-38ce-46a1-877f-80a7975340a0
# ╟─963e421c-43fb-43d3-b667-1b9912f940b8
# ╟─6cc62cf0-cb30-4d93-aad6-2ab16f60f95f
# ╟─2507c335-2886-42b5-b6b8-63279a2d60fe
# ╟─897a49b2-9763-4020-a476-5e0fccda1cfb
# ╟─229d395f-f3b4-40df-a482-264d540be875
# ╟─d448431e-13a9-440c-9463-9174d7400cf1
# ╟─7fc02644-b037-425a-b660-cc6904a95037
# ╟─75c98143-dd36-4446-adfa-c440fdf8aab1
# ╟─5f1085c7-f78a-433f-853e-0e3505fd99f9
# ╟─0c1448a9-52e7-41c9-b19a-e548d92db948
# ╟─0b188f96-87b1-41e4-ba66-1fa4d5252bd8
# ╟─d123161e-49f1-11ec-1c1b-51871624545d
