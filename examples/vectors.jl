### A Pluto.jl notebook ###
# v0.17.1

using Markdown
using InteractiveUtils

# ╔═╡ 71b1447b-39d9-46e2-966d-a1e6e8dcccc6
begin
	using Pkg
	Pkg.activate()
	
	using MeshArrays, OceanStateEstimation, MITgcmTools
	using JLD2
	import PlutoUI
	import CairoMakie as Mkie
	"Done with packages"
end

# ╔═╡ 7046b9b2-06bf-4ed8-80a3-79d2c52dacaf
md"""# Vector Field Operations

The flow of various quantities within the climate system is often represented as gridded vector fields. On [C-grids](https://en.wikipedia.org/wiki/Arakawa_grids), which are the most commonly used grids in global models, velocity points are staggered (shifted) by half a gridded point compared to tracer points (e.g., temperature fields). 

Vector fields on the C-grid are natively oriented along the grid cell axess. That is they are locally perpendicular to the grid cell edge where they are located. Often for plotting purposes or simply to facilitate analysis of model output by users, it is convenient to rotate these vector fields along Eastward and Northward axes. Doing this requires averaging them to grid cell centers first and then applying a rotation. Function `UVtoUEVN` is provided to this end. 

The result can be easily be interpolated  to any geographic position using `Interpolate`, and `InterpolationFactors` provides a mean to speed-up this step when batch processing a lot of output by pre-computing interpolation factor once and for all.

Another processing step that is often useful is to convert velocities (e.g. in m/s) to transports (e.g. in m^3/s) by multiplying each value by the corresponding grid cell area. That's what function `UVtoTransport` does.

The result can then be passed to functions like `ScalarPotential` and `VectorPotential` that are useful to isolate particular transport components. These higher level functions are illustrated later on in this tutorial.

Transports can also be integrated along model grid line paths that may for examples connect two geographic locations, or track latitude circles, as explained in [Forget et al 2015](https://doi.org/10.5194/gmd-8-3071-2015). Functionalities to this end include `LatitudeCircles` (compute grid edge paths that track latitude circles) and `ThroughFlow` (integrate transports accross such grid edge paths).
"""

# ╔═╡ 3b6fcab2-ae44-41c7-adec-424025092419
PlutoUI.TableOfContents()

# ╔═╡ 1e7616f8-65ec-4cfa-960f-bd26451756bc
md"""## Eastward U , Northward V"""

# ╔═╡ 5d25c473-eb1f-4040-86c7-86dee7add5df
md"""## Meridional Transports

Two widely used ways of presenting global scale transports are demonstrated here :

- the so-called Meridional Overturning Circulation (`MOC`) is often defined by a streamfunction [e.g. Fig 1 in Rousselet et al 2021](https://doi.org/10.1126/sciadv.abf5478). The simplest approach consists in integrating transports across lines of constant latitude (e.g., using `ThroughFlow`) and then vertically starting from the sea floor as shown below.
- a related computation is that of net meridional transports which were recently discussed in detail for heat transport in [Forget and Ferreira 2019](https://doi.org/10.1038/s41561-019-0333-7). The case of sea water volume (~ proportional to mass) is shown below which also corresponds to the meridional streamfunction once intergrated all the way to the sea surface.
- a third related computation, maybe more suprisingly, is the computation of transports convergence ξ. In the time mean, ξ corresponds in first approximation to the map of evaporation-precipitation-runoff plus a globally homogeneous rate of sea level change. Integrating the anomaly in ξ from the South Pole one in turn recovers the total Northward transport (see [Forget and Ferreira 2019](https://doi.org/10.1038/s41561-019-0333-7) for more on this).

!!! note
    Results are here plotted in Sverdrup units defined as 1 Sv = $10^6$ m$^3$/s
"""

# ╔═╡ 1cd94bce-9b10-4957-95d0-f396264df162
md"""## Helmholtz Decomposition

By applying the [Helmholtz Decomposition](https://en.wikipedia.org/wiki/Helmholtz_decomposition) to vertically integrated transports over the Global Ocean, we can decompose this vector field into rotational and divergent components. The divergent component gives the above map of inferred Evaporation-Precipitation-Runoff. The rotational component in turn gives the streamfunction shown below.
"""

# ╔═╡ e59c1819-52fd-49eb-8b0c-63fbbc0baffb
md"""## x-direction U, y-direction V

The maps below reflect the native velocity field orientation. Note how this differ from the Eastward U , Northward V plots shown at the top. The apparent discontinuities in the maps below simply reflect that neighboring subdomains can have different orientaton as is evident in [Forget et al 2015](https://doi.org/10.5194/gmd-8-3071-2015). 

!!! note
    Normally before interpolation one wants to apply rotation to convert to Eastward U , Northward V (see at the top of this notebook). The plot below reflects what happens if one doesn't.
"""

# ╔═╡ 0611986f-3d6e-4236-b68a-747b6018d381
md"""## Gradient And Curl"""

# ╔═╡ 87a4a232-d8eb-4612-8cc5-dc9d1c16dbe8
md"""## Particle Tracking

Materials and particles that tend to follow ocean currents or atmospheric winds can be analyzed in terms of trajectories. These are simply computed by integrating velocities over time within a [Lagrangian framework](https://en.wikipedia.org/wiki/Lagrangian_and_Eulerian_specification_of_the_flow_field).

When vector fields are represented using [MeshArrays.jl](https://github.com/JuliaClimate/MeshArrays.jl), as done in this tutorial, this is easily done via the [IndividualDisplacements.jl](https://github.com/JuliaClimate/IndividualDisplacements.jl) package as shown in [these examples](https://juliaclimate.github.io/IndividualDisplacements.jl/dev/examples/).
"""

# ╔═╡ a8ebbe41-1ac8-44da-b8c4-cbfd4d422227
md"""## Preliminary Steps

Here we load `MeshArrays.jl` and packages that provide a climatological estimate of ocean currents on the global grid defined in [Forget et al 2015](https://doi.org/10.5194/gmd-8-3071-2015) which is then used to demonstrate vector operations listed above.
"""

# ╔═╡ 3c2de9d4-b091-464b-9210-aa84f3d4c5f1
begin
	pth=MeshArrays.GRID_LLC90
	γ=GridSpec("LatLonCap",pth)
	Γ=GridLoad(γ;option="full")
	LC=LatitudeCircles(-89.0:89.0,Γ)
	"Done with grid"
end

# ╔═╡ c2788db6-3ab4-4c22-abf0-ae701a57e94d
begin
	μ =Γ.hFacC[:,1]
	μ[findall(μ.>0.0)].=1.0
	μ[findall(μ.==0.0)].=NaN

	if !isfile(joinpath(tempdir(),"interp_coeffs.jld2"))
		lon=[i for i=-179.:2.0:179., j=-89.:2.0:89.]
		lat=[j for i=-179.:2.0:179., j=-89.:2.0:89.]
		
		(f,i,j,w)=InterpolationFactors(Γ,vec(lon),vec(lat))
		jldsave(joinpath(tempdir(),"interp_coeffs.jld2"); lon=lon, lat=lat, f=f, i=i, j=j, w=w)
	end

	λ = load(joinpath(tempdir(),"interp_coeffs.jld2"))
	λ = MeshArrays.Dict_to_NamedTuple(λ)

	"Done with interpolation coefficients"
end

# ╔═╡ cf99ec64-d142-42ff-9767-3d851229024e
begin
	OceanStateEstimation.get_ecco_velocity_if_needed()
	
	"""
	    read_velocities(γ::gcmgrid,t::Int,pth::String)
	
	Read velocity components `u,v` from files in `pth`for time `t`
	"""
	function read_velocities(γ::gcmgrid,t::Int,pth::String)
	    u=read_nctiles("$pth"*"UVELMASS/UVELMASS","UVELMASS",γ,I=(:,:,:,t))
	    v=read_nctiles("$pth"*"VVELMASS/VVELMASS","VVELMASS",γ,I=(:,:,:,t))
	    return u,v
	end

	"Ready to read velocity fields"
end

# ╔═╡ ea72e6df-68da-4ade-b1d9-b7689e511a40
begin
	U=0.0*Γ.hFacW #varmeta not quite correct
	V=0.0*Γ.hFacS #varmeta not quite correct
	for t=1:12
		(u,v)=read_velocities(Γ.XC.grid,t,ECCOclim_path)
		for i in eachindex(U)
			U[i]=U[i]+u[i]/12.0
			V[i]=V[i]+v[i]/12.0
		end
	end

	(Utr,Vtr)=MeshArrays.UVtoTransport(U,V,Γ)

	"Done reading velocities, and computing transports"
end

# ╔═╡ dda89e9a-4585-49df-82cd-72db00c341fa
begin
	(uE,vN)=UVtoUEVN(U[:,1],V[:,1],Γ)
	uI=reshape(Interpolate(uE,λ.f,λ.i,λ.j,λ.w),size(λ.lon))
	vI=reshape(Interpolate(vN,λ.f,λ.i,λ.j,λ.w),size(λ.lon))
	"Done with Eastward / Westward components"
end

# ╔═╡ 60bda89f-5452-4f89-b4a4-a89568d45a6b
let
		fig1 = Mkie.Figure(resolution = (900,400),markersize=0.1)
		ax1 = Mkie.Axis(fig1[1,1], title="Eastward velocity (in m/s)")
		hm1=Mkie.contourf!(ax1,λ.lon[:,1],λ.lat[1,:],uI,levels=(-0.5:0.1:0.5)./2.0,clims=(-0.5,0.5)./2.0)
		Mkie.Colorbar(fig1[1,2], hm1, height = Mkie.Relative(0.65))
		fig1
end

# ╔═╡ 99e7e24b-578e-47be-818a-623c8e9e4381
let
		fig1 = Mkie.Figure(resolution = (900,400),markersize=0.1)
		ax1 = Mkie.Axis(fig1[1,1], title="Northward velocity (in m/s)")
		hm1=Mkie.contourf!(ax1,λ.lon[:,1],λ.lat[1,:],vI,levels=(-0.5:0.1:0.5)./5.0,clims=(-0.5,0.5)./5.0)
		Mkie.Colorbar(fig1[1,2], hm1, height = Mkie.Relative(0.65))
		fig1
end

# ╔═╡ 09dff33d-1c16-4fa6-a5af-8eb3c15816f3
begin
	nz=size(Γ.hFacC,2); nt=12; nl=length(LC)
	ov=Array{Float64,2}(undef,nl,nz)
	
	#integrate across latitude circles
	for z=1:nz
		UV=Dict("U"=>Utr[:,z],"V"=>Vtr[:,z],"dimensions"=>["x","y"])
		[ov[l,z]=ThroughFlow(UV,LC[l],Γ) for l=1:nl]
	end
	
	#integrate from bottom
	ov=reverse(cumsum(reverse(ov,dims=2),dims=2),dims=2)

	"Done with overturning"
end

# ╔═╡ a6f70839-b79f-42e1-bf4a-8c8978e6618e
let
	x=vec(-89.0:89.0); y=reverse(vec(Γ.RF[1:end-1])); #coordinate variables
	z=reverse(ov,dims=2); z[z.==0.0].=NaN

	fig1 = Mkie.Figure(resolution = (900,400),markersize=0.1)
	ax1 = Mkie.Axis(fig1[1,1], title="Meridional Overturning Streamfunction (in Sv)")
	hm1=Mkie.contourf!(ax1,x,y,1e-6*z,levels=(-40.0:5.0:40.0),clims=(-40,40))
	Mkie.Colorbar(fig1[1,2], hm1, height = Mkie.Relative(0.65))
	fig1
	#savefig("MOC_mean.png")
end

# ╔═╡ 923568dc-85f3-43c1-a631-d408a1662fb3
begin
	MT=fill(0.0,nl)
	
	#integrate across latitude circles and depth
	for z=1:nz
		UV=Dict("U"=>Utr[:,z],"V"=>Vtr[:,z],"dimensions"=>["x","y"])
		[MT[l]=MT[l]+ThroughFlow(UV,LC[l],Γ) for l=1:nl]
	end
	
	"Done with meridional transport"
end

# ╔═╡ 2d6c1787-3c30-4daf-8575-b51211f8e860
let
	x=vec(-89.0:89.0)
	fig1 = Mkie.Figure(resolution = (900,400),markersize=0.1)
	ax1 = Mkie.Axis(fig1[1,1], title="Northward Volume Transport (in Sv)")
	hm1=Mkie.lines!(x,1e-6*MT,xlabel="latitude",ylabel="Transport (in Sv)",label="ECCO estimate")
	fig1
end

# ╔═╡ 3a3ac112-cd2b-4ac9-9029-b79b8b5fa8b4
begin
	Tx=0.0*Utr[:,1]
	Ty=0.0*Vtr[:,1]
	for z=1:nz
		Tx=Tx+Utr[:,z]
		Ty=Ty+Vtr[:,z]
	end

	#convergence & land mask
	TrspCon=μ.*convergence(Tx,Ty)
	
	#scalar potential
	TrspPot=ScalarPotential(TrspCon)
	
	#Divergent transport component
	(TxD,TyD)=gradient(TrspPot,Γ)
	TxD=TxD.*Γ.DXC
	TyD=TyD.*Γ.DYC
	
	#Rotational transport component
	TxR = Tx-TxD
	TyR = Ty-TyD
	
	#vector Potential
	TrspPsi=μ.*VectorPotential(TxR,TyR,Γ)

	"Done with convergence and potentials"
end

# ╔═╡ 9928d59f-0733-4ad9-93e7-f1a8622afd1f
let
		z=reshape(Interpolate(TrspCon,λ.f,λ.i,λ.j,λ.w),size(λ.lon))

		fig1 = Mkie.Figure(resolution = (900,400),markersize=0.1)
		ax1 = Mkie.Axis(fig1[1,1], title="Evaporation - Precipitation - Runoff (in Sv; inferred)")
		hm1=Mkie.contourf!(ax1,λ.lon[:,1],λ.lat[1,:],z,levels=(-1.0:0.25:1.0).*800.0,clims=(-1.0,1.0).*800.0)
		Mkie.Colorbar(fig1[1,2], hm1, height = Mkie.Relative(0.65))
		fig1
end

# ╔═╡ c3b3a104-a1d6-4692-ae0c-2e55821afd03
let
	z=reshape(Interpolate(1e-6*TrspPsi,λ.f,λ.i,λ.j,λ.w),size(λ.lon))

	fig1 = Mkie.Figure(resolution = (900,400),markersize=0.1)
	ax1 = Mkie.Axis(fig1[1,1], title="Streamfunction (in Sv)")
	hm1=Mkie.contourf!(ax1,λ.lon[:,1],λ.lat[1,:],z,levels=(-1.0:0.25:1.0).*40.0,clims=(-1.0,1.0).*40.0)
	Mkie.Colorbar(fig1[1,2], hm1, height = Mkie.Relative(0.65))
	fig1
end

# ╔═╡ 210b5414-ada8-49a1-8ac4-9115ed33e285
let
	(dDdx, dDdy)=gradient(TrspPot,Γ)
	dDdx=reshape(Interpolate(μ.*dDdx,λ.f,λ.i,λ.j,λ.w),size(λ.lon))
	dDdy=reshape(Interpolate(μ.*dDdy,λ.f,λ.i,λ.j,λ.w),size(λ.lon))

	fig1 = Mkie.Figure(resolution = (900,600),markersize=0.1)
	ax1 = Mkie.Axis(fig1[1,1], title="Gradient of scalar potential in x-direction (in 1/s)")
	hm1=Mkie.contourf!(ax1,λ.lon[:,1],λ.lat[1,:],dDdx,levels=(-1.0:0.25:1.0).*0.1)
	ax1 = Mkie.Axis(fig1[2,1], title="Gradient of scalar potential in y-direction (in 1/s)")
	hm1=Mkie.contourf!(ax1,λ.lon[:,1],λ.lat[1,:],dDdy,levels=(-1.0:0.25:1.0).*0.1)
	Mkie.Colorbar(fig1[1:2,2], hm1, height = Mkie.Relative(0.65))
	fig1
end

# ╔═╡ e033825c-019c-44b6-86dd-30fff79f8339
let
	dDdx=U[:,1]
	dDdy=V[:,1]

	DD=Interpolate(dDdx,λ.f,λ.i,λ.j,λ.w)
	DD=reshape(DD,size(λ.lon))
	DD[findall(DD.==0.0)].=NaN
	dDdx_i=DD

	DD=Interpolate(dDdy,λ.f,λ.i,λ.j,λ.w)
	DD=reshape(DD,size(λ.lon))
	DD[findall(DD.==0.0)].=NaN
	dDdy_i=DD

	fig = Mkie.Figure(resolution = (900,600), backgroundcolor = :grey95)
	ax = Mkie.Axis(fig[1,1], title="x-direction velocity (in m/s)",xlabel="longitude",ylabel="latitude")
	hm1=Mkie.heatmap!(ax,λ.lon[:,1],λ.lat[1,:],dDdx_i,colorrange=(-1.0,1.0).*0.2)
	ax = Mkie.Axis(fig[2,1], title="y-direction velocity (in m/s)",xlabel="longitude",ylabel="latitude")
	hm1=Mkie.heatmap!(ax,λ.lon[:,1],λ.lat[1,:],dDdy_i,colorrange=(-1.0,1.0).*0.2)
	Mkie.Colorbar(fig[1:2,2], hm1, height = Mkie.Relative(0.65))
	fig
end

# ╔═╡ 0888db7a-ff14-44c0-9c20-fed722f7e41e
let
	u=U[:,1]
	v=V[:,1]
	tmp=curl(u,v,Γ)

	z=reshape(Interpolate(μ.*tmp,λ.f,λ.i,λ.j,λ.w),size(λ.lon))

	fig1 = Mkie.Figure(resolution = (900,400),markersize=0.1)
	ax1 = Mkie.Axis(fig1[1,1], title="Curl Of velocity Field (in 1/s)")
	hm1=Mkie.contourf!(ax1,λ.lon[:,1],λ.lat[1,:],z,levels=(-1.0:0.25:1.0).*1e-11)
	Mkie.Colorbar(fig1[1,2], hm1, height = Mkie.Relative(0.65))
	fig1
end

# ╔═╡ Cell order:
# ╟─7046b9b2-06bf-4ed8-80a3-79d2c52dacaf
# ╟─3b6fcab2-ae44-41c7-adec-424025092419
# ╟─1e7616f8-65ec-4cfa-960f-bd26451756bc
# ╟─dda89e9a-4585-49df-82cd-72db00c341fa
# ╟─60bda89f-5452-4f89-b4a4-a89568d45a6b
# ╟─99e7e24b-578e-47be-818a-623c8e9e4381
# ╟─5d25c473-eb1f-4040-86c7-86dee7add5df
# ╟─09dff33d-1c16-4fa6-a5af-8eb3c15816f3
# ╟─a6f70839-b79f-42e1-bf4a-8c8978e6618e
# ╟─923568dc-85f3-43c1-a631-d408a1662fb3
# ╟─2d6c1787-3c30-4daf-8575-b51211f8e860
# ╟─3a3ac112-cd2b-4ac9-9029-b79b8b5fa8b4
# ╟─9928d59f-0733-4ad9-93e7-f1a8622afd1f
# ╟─1cd94bce-9b10-4957-95d0-f396264df162
# ╟─c3b3a104-a1d6-4692-ae0c-2e55821afd03
# ╟─e59c1819-52fd-49eb-8b0c-63fbbc0baffb
# ╟─e033825c-019c-44b6-86dd-30fff79f8339
# ╟─0611986f-3d6e-4236-b68a-747b6018d381
# ╟─210b5414-ada8-49a1-8ac4-9115ed33e285
# ╟─0888db7a-ff14-44c0-9c20-fed722f7e41e
# ╟─87a4a232-d8eb-4612-8cc5-dc9d1c16dbe8
# ╟─a8ebbe41-1ac8-44da-b8c4-cbfd4d422227
# ╟─71b1447b-39d9-46e2-966d-a1e6e8dcccc6
# ╟─3c2de9d4-b091-464b-9210-aa84f3d4c5f1
# ╟─c2788db6-3ab4-4c22-abf0-ae701a57e94d
# ╟─cf99ec64-d142-42ff-9767-3d851229024e
# ╟─ea72e6df-68da-4ade-b1d9-b7689e511a40
