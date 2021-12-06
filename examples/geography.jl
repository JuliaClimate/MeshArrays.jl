### A Pluto.jl notebook ###
# v0.17.2

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
	using Pkg
	Pkg.activate()
	using MeshArrays, PlutoUI, Statistics, JLD2
	import CairoMakie as Mkie
	"Done with packages"
end

# ╔═╡ cc3e9d0c-8b71-432b-b68d-a8b832ca5f26
md"""# Geography and Visualization

This tutorial focuses on the relation between the gridded space in [MeshArrays.jl](https://juliaclimate.github.io/MeshArrays.jl/dev/) and the geospatial space (longitude, lattitude). A common use case is for plotting in the form of maps. Going back and forth between the two is also often useful for analysis of simulations (transports, averages, etc) and e.g. simulating particle trajectories.
"""

# ╔═╡ 39788eec-dd8b-4d65-a0f9-ea73a2d6690c
TableOfContents()

# ╔═╡ ca8a8f1b-a225-46c4-93c1-ce1a4b016461
md"""## Interpolation scheme

Here we interpolate from the global grid (a `MeshArray`) to a set of arbitary locations (longitude, latitude pairs) as is commonly done e.g. to compare climate models to sparse field observations, or simply to plot gridded fields in geographic coordinates as done below for ocean depth.

Later in this notebook, we break down the `MeshArrays.Interpolate` function. In brief, the program finds a grid point quadrilateral (4 grid points) that encloses the chosen location. Computation is chuncked in subdomains (tiles) to allow for parallelism. `MeshArrays.InterpolationFactors` outputs interpolation coefficients -- reusing those is easy and fast.
"""

# ╔═╡ 2a02b528-5d8e-47a3-9fa9-83872d57ee6f
@doc MeshArrays.Interpolate

# ╔═╡ 6aab5feb-06b6-4fbd-8732-83b70f397f00
begin
	#push!(lonPairs,[-68 -63]); push!(latPairs,[-54 -66]); push!(namPairs,"Drake Passage");

	lon1_slct = @bind lon1 Select(-179.0:180.0;default=20.0)
	lon2_slct = @bind lon2 Select(-179.0:180.0;default=-63.0)
	lat1_slct = @bind lat1 Select(-89.0:89.0;default=-20.0)
	lat2_slct = @bind lat2 Select(-89.0:89.0;default=-66.0)
	md"""## Shortest Path (A to B)

The shortest path between two locations on the surface of a sphere is the one that follows a [great circle](https://en.wikipedia.org/wiki/Great-circle_distance). The method used below identifies the path along grid cell edges that corresponds closely to the great circle path.

In `MeshArrays.jl`, we provide the `Transect` function to this end. It relies on the `MeshArrays.edge_mask` function, which can also be used find paths that approximate circles of constant latitude, or delineate a subdomain for example.

| Location 1         | Location 2         | 
|:------------------:|:------------------:|
| lon = $(lon1_slct) | lon = $(lon2_slct) | 
| lat = $(lat1_slct) | lat = $(lat2_slct) | 
	"""
end

# ╔═╡ 0cfd749e-af65-4716-86c0-0a587b9647f9
@doc MeshArrays.edge_mask

# ╔═╡ c81ad3bb-a924-4393-bc73-0607cfebf75f
begin
	faceID_slct = @bind faceID Select(1:5;default=1)
	
	md"""## Grid Cells and Areas

	When the goal is to compute averages or integrals over a subdomain or e.g. the whole Ocean, it is important to account for the corresponding land mask and grid cell areas. Here we visualize variations in grid cell areas for a specific global grid. Details will differ for other grids, but should not be overlooked.
	
	- select a grid subdomain : $(faceID_slct)
	"""
end

# ╔═╡ 25144e1b-21fc-4cc9-b63d-7b26eab1a673
md"""## Appendices"""

# ╔═╡ 75c98143-dd36-4446-adfa-c440fdf8aab1
function setup_interp(Γ)
	μ =Γ.hFacC[:,1]
	μ[findall(μ.>0.0)].=1.0
	μ[findall(μ.==0.0)].=NaN

	if !isfile(joinpath(tempdir(),"interp_coeffs_halfdeg.jld2"))
		lon=[i for i=-179.75:0.5:179.75, j=-89.75:0.5:89.75]
		lat=[j for i=-179.75:0.5:179.75, j=-89.75:0.5:89.75]
		
		(f,i,j,w)=InterpolationFactors(Γ,vec(lon),vec(lat))
		jldsave(joinpath(tempdir(),"interp_coeffs_halfdeg.jld2"); 
			lon=lon, lat=lat, f=f, i=i, j=j, w=w, μ=μ)
	end

	λ = load(joinpath(tempdir(),"interp_coeffs_halfdeg.jld2"))
	λ = MeshArrays.Dict_to_NamedTuple(λ)
end

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

# ╔═╡ be38ff51-3526-44a0-9d8c-9209355e4a4a
begin
	λ=setup_interp(Γ)
	DD=Interpolate(λ.μ*Γ.Depth,λ.f,λ.i,λ.j,λ.w)
	DD=reshape(DD,size(λ.lon))
	#DD[findall(DD.==0.0)].=NaN
	"Done with interpolating Γ.Depth"
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

# ╔═╡ ed36a2a5-44ea-43a7-a3bd-13f234b6580d
let	
	fig = Mkie.Figure(resolution = (900,600), backgroundcolor = :grey95)
	ax = Mkie.Axis(fig[1,1],xlabel="longitude",ylabel="latitude",title="Ocean Depth Shown in Colors (interpolated Γ.Depth)")
	hm1=Mkie.heatmap!(ax,λ.lon[:,1],λ.lat[1,:],DD,colormap=:spring,colorrange=(0.0,6000.0))

	Mkie.scatter!(ax,locCwhole.lon[:],locCwhole.lat[:],color=:blue,markersize=2.0)
	Mkie.scatter!(ax,locCshort.lon[:],locCshort.lat[:],color=:black,markersize=4.0)
	Mkie.scatter!(ax,lons[:],lats[:],color=:blue)

	#Mkie.scatter!(ax,locClat.lon[:],locClat.lat[:],color=:cyan,markersize=4.0)

	#Mkie.Colorbar(fig[1,2], hm1, height = Mkie.Relative(0.65))

	fig
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

# ╔═╡ 1477cd5d-7ee3-4af8-95cc-13309db00520
begin
	fig = Mkie.Figure(resolution = (900,600), backgroundcolor = :grey95)
	ax = Mkie.Axis(fig[1,1],xlabel="longitude",ylabel="latitude",title="grid cell area")

	mx=maximum(Γ.RAC)
	for ff in 1:length(Γ.RAC)
		col=λ.μ[ff][:].*Γ.RAC[ff][:]
		kk=findall((!isnan).(col))
		if ff==faceID
			Mkie.scatter!(ax,Γ.XC[ff][kk],Γ.YC[ff][kk],color=col[kk],
				colorrange = (0.0, mx),markersize=2.0)
		else
			Mkie.scatter!(ax,Γ.XC[ff][kk],Γ.YC[ff][kk],color=:gray,markersize=2.0)
		end
	end
	Mkie.Colorbar(fig[1,2], colorrange=(0.0, mx), height = Mkie.Relative(0.65))

	fig
end

# ╔═╡ b0d576fc-971a-47c7-9a57-f2c788083bcd
begin

	basin_slct = @bind basin_nam Select(basin_list;default=basin_list[1])

	#basinID_slct = @bind basinID Select(1:maximum(basins);default=1)
	
	md"""## Ocean Basins / Regions

	There are many ways to split up the global domain into Seas and Oceans. Here is one example.
	
	- select a basin nam : $(basin_slct)
	"""
end

# ╔═╡ 31756b2e-20df-47c9-aaa8-5583e6a81267
let
	fig = Mkie.Figure(resolution = (900,600), backgroundcolor = :grey95)
	ax = Mkie.Axis(fig[1,1],xlabel="longitude",ylabel="latitude",title=basin_nam* " (shown in red)")

	basinID=findall(basin_list.==basin_nam)[1]
	
	mx=maximum(basins)
	for ff in 1:length(Γ.RAC)
		col=λ.μ[ff][:].*(basins[ff][:].==basinID)
		kk=findall(col.>0.0)
		!isempty(kk) ? Mkie.scatter!(ax,Γ.XC[ff][kk],Γ.YC[ff][kk],color=:red,markersize=2.0) : nothing
		kk=findall((col.==0.0).*(!isnan).(λ.μ[ff][:]))
		!isempty(kk) ? Mkie.scatter!(ax,Γ.XC[ff][kk],Γ.YC[ff][kk],color=basins[ff][kk],
			colorrange=(0.0,mx),markersize=2.0) : nothing
	end
	Mkie.Colorbar(fig[1,2], colorrange=(0.0, mx), height = Mkie.Relative(0.65))

	fig
end

# ╔═╡ 94b8ba05-dfb8-4075-a260-7968e8fdd78f
md"""#### Compute All Paths"""

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
let
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
end

# ╔═╡ Cell order:
# ╟─cc3e9d0c-8b71-432b-b68d-a8b832ca5f26
# ╟─39788eec-dd8b-4d65-a0f9-ea73a2d6690c
# ╟─ca8a8f1b-a225-46c4-93c1-ce1a4b016461
# ╟─be38ff51-3526-44a0-9d8c-9209355e4a4a
# ╟─2a02b528-5d8e-47a3-9fa9-83872d57ee6f
# ╟─6aab5feb-06b6-4fbd-8732-83b70f397f00
# ╟─ed36a2a5-44ea-43a7-a3bd-13f234b6580d
# ╟─963e421c-43fb-43d3-b667-1b9912f940b8
# ╟─0cfd749e-af65-4716-86c0-0a587b9647f9
# ╟─897a49b2-9763-4020-a476-5e0fccda1cfb
# ╟─c81ad3bb-a924-4393-bc73-0607cfebf75f
# ╟─1477cd5d-7ee3-4af8-95cc-13309db00520
# ╟─b0d576fc-971a-47c7-9a57-f2c788083bcd
# ╟─31756b2e-20df-47c9-aaa8-5583e6a81267
# ╟─25144e1b-21fc-4cc9-b63d-7b26eab1a673
# ╟─75c98143-dd36-4446-adfa-c440fdf8aab1
# ╟─d123161e-49f1-11ec-1c1b-51871624545d
# ╟─39924391-38ce-46a1-877f-80a7975340a0
# ╟─94b8ba05-dfb8-4075-a260-7968e8fdd78f
# ╟─2507c335-2886-42b5-b6b8-63279a2d60fe
# ╟─6cc62cf0-cb30-4d93-aad6-2ab16f60f95f
