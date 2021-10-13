### A Pluto.jl notebook ###
# v0.16.1

using Markdown
using InteractiveUtils

# This Pluto notebook uses @bind for interactivity. When running this notebook outside of Pluto, the following 'mock version' of @bind gives bound variables a default value (instead of an error).
macro bind(def, element)
    quote
        local el = $(esc(element))
        global $(esc(def)) = Core.applicable(Base.get, el) ? Base.get(el) : missing
        el
    end
end

# ╔═╡ 1a714fba-2a8e-11ec-182f-8f85cc17b02a
begin
	using Pkg
	Pkg.activate("../docs")
	
	using MeshArrays
		
	import PlutoUI
	import CairoMakie as Mkie
	toc=PlutoUI.TableOfContents()
	
	md"""
	$(toc)
	## Use Julia Packages
	
	"""
end

# ╔═╡ 0358ae4b-f941-4254-aa10-6c37bb09c464
md"""# Get Started

!!! tip
    You can download this notebook or simply read its code in detail [here](https://raw.githubusercontent.com/JuliaClimate/MeshArrays.jl/master/examples/basics.jl). Running notebooks interactively using [Pluto.jl](https://plutojl.org) can be a good way to get familiar with software such as this one.

## Create Your First `MeshArray`

To create your first `MeshArray`, open the `Julia` REPL and type:

```
using MeshArrays
tmp=MeshArray(randn(20,10))
```
"""

# ╔═╡ b1f3ac4b-5f75-40fa-b125-186556ece946
begin
	tmp=MeshArray(randn(20,10))
	PlutoUI.with_terminal() do
		show(tmp)
	end
end

# ╔═╡ a0a7637f-f115-4257-a140-06dc48dd5fba
md"""## Tutorial

For the **first part** of this tutorial, we focus on a use case which proceeds as follows:

1. generate global grid configuration
1. initialize random noise field
1. apply smoothing
1. plot results

Step 3 illustrates how `MeshArrays.jl` computes partial derivatives between neighboring subdomains of the global domain. The chosen smoothing method (`MeshArrays.smooth`) indeed integrates lateral diffusion, a transport process, over time and the global domain. Other processes like e.g. advection by ocean currents would work similarly.

This example is then repeated twice with different grid configurations commonly used in Earth system models (see [Grids](https://juliaclimate.github.io/MeshArrays.jl/dev/start/#Grids)). The repeated code is provided by the `comp` (steps 2 and 3) and `viz` (step 4) functions.

In the **second part** of this tutorial, we further document `Array` operations as performed on a `MeshArray`, and some of the more advanced functionalities that `MeshArray.jl` provides to compute transports within the climate system. 
"""

# ╔═╡ 54030bff-45f7-43b1-b19f-7c0cd05a02ae
begin
	nP_s = @bind nP PlutoUI.NumberField(10:20:30, default=20)
	nQ_s = @bind nQ PlutoUI.NumberField(10:20:30, default=30)
	mP_s = @bind mP PlutoUI.NumberField(1:4, default=4)
	mQ_s = @bind mQ PlutoUI.NumberField(1:4, default=2)
	md"""### 1. Doubly Periodic Domain
	
	Let's start with a doubly periodic domain split into `mP*mQ` subdomains. 
	Each subdomain corresponds to an array of `nP,nQ` grid points. 
	The `UnitGrid()` function is used to generate such a grid configuration.
	
	nP x mP = $(nP_s) x $(mP_s)
		
	nQ x mQ = $(nQ_s) x $(mQ_s)
	
	`(Γ,γ)=UnitGrid( (nP*mP,nQ*mQ) , (nP,nQ) ; option="full")`
	"""
end

# ╔═╡ 1fc99fca-3440-4765-b11d-e10eb560aff7
begin
	(Γ,γ)=UnitGrid( (nP*mP,nQ*mQ) , (nP,nQ) ; option="full")
	
	PlutoUI.with_terminal() do
		show(Γ.XC)
	end
end

# ╔═╡ 2ee0a449-fd98-43e1-90b4-b4e15f237d67
md"""
### 2. Cube Sphere Grid

Now, we instead use a grid that has 6 subdomains, 32x32 points each, covering the six faces of a cube. We note that this _cube sphere_ topology involves connections between subdomain that are slightly more complicated than in the first example. 

Here, grid variables are read from the `MeshArrays.GRID_CS32` folder by the `GridLoad` function. The `GridSpec` function provides the corresponding subdomain sizes. Again we call `comp()` which combines steps 2 and 3, and `viz()` to visualize. 
"""

# ╔═╡ d37ae153-8646-44ae-abd9-7e0e7d479275
md"""### 3. LLC90 grid

The [Lat-Lon-Cap grid](http://www.geosci-model-dev.net/8/3071/2015/) (or LLC) is a global ocean model grid which is widely used in the [MITgcm user community](https://mitgcm.readthedocs.io/en/latest/). It has 5 uneven subdomains, variable grid spacing, and continents [(Forget et al 2015)](http://www.geosci-model-dev.net/8/3071/2015/). LLC90's resolution is one degree albeit with modications in the Arctic and along the Equator.

Here, grid variables are read from the `MeshArrays.GRID_LLC90` folder by the `GridLoad` function. The `GridSpec` function provides the corresponding subdomain sizes. Again we call `comp()` which combines steps 2 and 3, and `viz()` to visualize. 
"""

# ╔═╡ 58f95665-9687-4b4f-af99-1239818f71a3
md"""## Helper Functions

The `comp` and `viz` funtions can be used repeadtedly with different grids.
"""

# ╔═╡ 62f2dc32-456e-4e8a-8c4a-3a0ae236af13
function comp(Γ::NamedTuple)
    γ=Γ.XC.grid
	T=eltype(Γ.XC)

    #initialize 2D field of random numbers
    tmp1=randn(T,Tuple(γ.ioSize))
    Rini=γ.read(tmp1,MeshArray(γ,T))

    #apply land mask
	length(size(Γ.hFacC))>1 ? tmp1=Γ.hFacC[:,1] : tmp=Γ.hFacC
	Rini[findall(tmp1.==0.0)].=NaN

    #specify smoothing length scales in x, y directions
    DXCsm=3*Γ.DXC; DYCsm=3*Γ.DYC;

    #apply smoother
    Rend=smooth(Rini,DXCsm,DYCsm,Γ);

    return Rini,Rend

end


# ╔═╡ cef83984-80dc-4196-8edd-314885fa9448
begin
	Rini_a,Rend_a=comp(Γ)

	md"""Steps 2 (`Rini`) and 3 (`Rend`) are completed in this code cell (using `comp` defined below) and results plotted in the next. 
	"""
end

# ╔═╡ 2573fda7-66cb-484f-9d06-a9bb85a26a9e
let	
	#visualize
	x=γ.write(Γ.XC)[:,1]
	y=γ.write(Γ.YC)[1,:]
	#x=0.5:ioSize[1]-0.5
	#y=0.5:ioSize[2]-0.5
	
	fig = Mkie.Figure(resolution = (900,900), backgroundcolor = :grey95)

	ax1 = Mkie.Axis(fig[1,1])
	hm1=Mkie.heatmap!(ax1,x,y,γ.write(Rini_a),clims=(-0.25,0.25),tickfont = (4, :black))
	ax2 = Mkie.Axis(fig[1,2])
	hm2=Mkie.heatmap!(ax2,x,y,γ.write(Rend_a),clims=(-0.25,0.25),tickfont = (4, :black))
	Mkie.Colorbar(fig[1,3], hm1, height = Mkie.Relative(0.65))
	
	fig
	md"""#### Visualize Global Domain
	
	The `heatmap` display should reflect that `Rend` (R.H.S.) is smoother, and therefore more muted, than is `Rini` (L.H.S.).
	
	!!! note
	    This visualization only uniquely defined for simple grids. Otherwise visualization subdomain by subdomain, or via interpolation are always possible.
	
	$(fig)
	"""
end

# ╔═╡ dacc2261-1d75-4d68-b82c-174e4fae3631
begin
	γ_b=GridSpec("CubeSphere",MeshArrays.GRID_CS32)
	Γ_b=GridLoad(γ_b;option="full")

	Rini_b,Rend_b=comp(Γ_b)
	PlutoUI.with_terminal() do
		show(Rend_b)
	end
end

# ╔═╡ 23198b8d-6bf0-4af9-91d0-7886c9574f81
begin
	γ_c=GridSpec("LatLonCap",MeshArrays.GRID_LLC90)
	Γ_c=GridLoad(γ_c;option="full")
	Rini_c,Rend_c=comp(Γ_c)
	PlutoUI.with_terminal() do
		show(Rend_c)
	end
end

# ╔═╡ 2c29ba59-ffc9-4763-8405-250029016ca5
begin
	let
		γ=γ_c
		Γ=Γ_c
		D=γ.read(γ.path*"Depth.data",MeshArray(γ,γ.ioPrec))

		#Array operations
		1000 .+ D
		D .+ 1000
		D+D
		D .- 1000
		1000 .- D
		D-D
		1000*D
		D*1000
		D*D
		D/1000
		1000 ./ D
		D/D
		
		#Indexing methods
		H=Γ.hFacC
		typeof(H)
	    size(H)
		eltype(H)
		H.f[1,40]
		H.f[:,40]
		view(H,:,40)
		
		#Exchange funtions
		Dexch=exchange(D,4)
		(dDdx, dDdy)=gradient(D,Γ)
		(dDdxEx,dDdyEx)=exchange(dDdx,dDdy,4)
	end
	
	md"""### 4. Array Operations
	
	This cell demonstrates array operations, indexing methods, and exchange functions.
	
	"""
end

# ╔═╡ 5c755d83-6729-4de4-8038-8619b53795ff
begin 
	function demo_Tr(γ,Γ)
		TrspX=γ.read(γ.path*"TrspX.bin",MeshArray(γ,Float32))
		TrspY=γ.read(γ.path*"TrspY.bin",MeshArray(γ,Float32))
		#TauX=γ.read(γ.path*"TauX.bin",MeshArray(γ,Float32))
		#TauY=γ.read(γ.path*"TauY.bin",MeshArray(γ,Float32))
		#SSH=γ.read(γ.path*"SSH.bin",MeshArray(γ,Float32))

		la=-89.0:89.0
	    LC=LatitudeCircles(la,Γ)
	    UV=Dict("U"=>TrspX,"V"=>TrspY,"dimensions"=>["x","y"])

		Tr=Array{Float64,1}(undef,length(LC));
		[Tr[i]=ThroughFlow(UV,LC[i],Γ) for i=1:length(LC)]
		
		fig = Mkie.Figure(resolution = (900,900), backgroundcolor = :grey95)
		ax = Mkie.Axis(fig[1,1], title="Meridional Ocean Transport", 
				ylabel="Sv (1 Sv = 10^6 m^3/s)",xlabel="latitude")
		Mkie.lines!(ax,la,Tr/1e6)
		fig
	end
	
	fig_Tr=demo_Tr(γ_c,Γ_c)
		
	md"""### 4. Transport Computations.
	
	This cell demonstrates methods that compute transports within the Earth System. Accurate transport computations involve factoring in the correct grid scaling factors, and integrating along non-trivial grid line paths [Forget et al 2015](https://doi.org/10.5194/gmd-8-3071-2015).
	
	$(fig_Tr)
	"""
end

# ╔═╡ 4bb0c25c-336c-4438-87df-a93344ea3cfb
function viz(fld)
	fig = Mkie.Figure(resolution = (900,900), backgroundcolor = :grey95)
	nf=length(fld.fSize)
	nn=Int(ceil(nf/2))
	ii=[i for j in 1:2, i in 1:nn]
	jj=[j for j in 1:2, i in 1:nn]
	
	for f in 1:nf
		ax = Mkie.Axis(fig[ii[f],jj[f]], title="face $(f)")

		s=fld.fSize[f]		
		x=collect(0.5:s[1]-0.5)
		y=collect(0.5:s[2]-0.5)
		z=fld[f]

		hm1=Mkie.heatmap!(ax,x,y,z,clims=(-0.25,0.25),tickfont = (4, :black))
	end
	
	Mkie.Colorbar(fig[1:3,3], limits=(-0.25,0.25), height = Mkie.Relative(0.65))
	
	fig
end

# ╔═╡ e34de61b-8630-43f8-84ba-ce7a99d673ab
begin
	f_a=viz(Rend_a)
	md"""#### Visualize Subdomain by Subdomain
	$(f_a)
	"""
end

# ╔═╡ aa95e258-f9d0-4137-8842-e878833740a6
viz(Rend_b)

# ╔═╡ 1ad0fdcc-9c0d-4d17-96d6-37b01610cf85
viz(Rend_c)

# ╔═╡ Cell order:
# ╟─0358ae4b-f941-4254-aa10-6c37bb09c464
# ╟─b1f3ac4b-5f75-40fa-b125-186556ece946
# ╠═1a714fba-2a8e-11ec-182f-8f85cc17b02a
# ╟─a0a7637f-f115-4257-a140-06dc48dd5fba
# ╟─54030bff-45f7-43b1-b19f-7c0cd05a02ae
# ╟─1fc99fca-3440-4765-b11d-e10eb560aff7
# ╟─cef83984-80dc-4196-8edd-314885fa9448
# ╟─2573fda7-66cb-484f-9d06-a9bb85a26a9e
# ╟─e34de61b-8630-43f8-84ba-ce7a99d673ab
# ╟─2ee0a449-fd98-43e1-90b4-b4e15f237d67
# ╟─dacc2261-1d75-4d68-b82c-174e4fae3631
# ╟─aa95e258-f9d0-4137-8842-e878833740a6
# ╟─d37ae153-8646-44ae-abd9-7e0e7d479275
# ╟─23198b8d-6bf0-4af9-91d0-7886c9574f81
# ╟─1ad0fdcc-9c0d-4d17-96d6-37b01610cf85
# ╠═2c29ba59-ffc9-4763-8405-250029016ca5
# ╟─5c755d83-6729-4de4-8038-8619b53795ff
# ╟─58f95665-9687-4b4f-af99-1239818f71a3
# ╠═62f2dc32-456e-4e8a-8c4a-3a0ae236af13
# ╠═4bb0c25c-336c-4438-87df-a93344ea3cfb
