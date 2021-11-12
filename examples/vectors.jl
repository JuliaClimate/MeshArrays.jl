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

# ╔═╡ 3b6fcab2-ae44-41c7-adec-424025092419
PlutoUI.TableOfContents()

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

	"Done reading velocities, and transports"
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

	fig1 = Mkie.Figure(resolution = (900,600),markersize=0.1)
	ax1 = Mkie.Axis(fig1[1,1], title="Overturning mean (Eulerian; in Sv)")
	hm1=Mkie.contourf!(ax1,x,y,1e-6*z,levels=(-40.0:5.0:40.0),clims=(-40,40))
	Mkie.Colorbar(fig1[1,2], hm1, height = Mkie.Relative(0.65))
	fig1
	#savefig("MOC_mean.png")
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
		fig1 = Mkie.Figure(resolution = (900,600),markersize=0.1)
		ax1 = Mkie.Axis(fig1[1,1], title="Eastward velocity (in m/s)")
		hm1=Mkie.contourf!(ax1,λ.lon[:,1],λ.lat[1,:],uI,levels=(-0.5:0.1:0.5)./2.0,clims=(-0.5,0.5)./2.0)
		Mkie.Colorbar(fig1[1,2], hm1, height = Mkie.Relative(0.65))
		fig1
end

# ╔═╡ 99e7e24b-578e-47be-818a-623c8e9e4381
let
		fig1 = Mkie.Figure(resolution = (900,600),markersize=0.1)
		ax1 = Mkie.Axis(fig1[1,1], title="Northward velocity (in m/s)")
		hm1=Mkie.contourf!(ax1,λ.lon[:,1],λ.lat[1,:],vI,levels=(-0.5:0.1:0.5)./5.0,clims=(-0.5,0.5)./5.0)
		Mkie.Colorbar(fig1[1,2], hm1, height = Mkie.Relative(0.65))
		fig1
end

# ╔═╡ Cell order:
# ╟─71b1447b-39d9-46e2-966d-a1e6e8dcccc6
# ╟─3b6fcab2-ae44-41c7-adec-424025092419
# ╟─3c2de9d4-b091-464b-9210-aa84f3d4c5f1
# ╟─c2788db6-3ab4-4c22-abf0-ae701a57e94d
# ╟─cf99ec64-d142-42ff-9767-3d851229024e
# ╟─ea72e6df-68da-4ade-b1d9-b7689e511a40
# ╟─09dff33d-1c16-4fa6-a5af-8eb3c15816f3
# ╟─a6f70839-b79f-42e1-bf4a-8c8978e6618e
# ╟─dda89e9a-4585-49df-82cd-72db00c341fa
# ╟─60bda89f-5452-4f89-b4a4-a89568d45a6b
# ╟─99e7e24b-578e-47be-818a-623c8e9e4381
