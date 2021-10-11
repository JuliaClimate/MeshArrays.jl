### A Pluto.jl notebook ###
# v0.16.1

using Markdown
using InteractiveUtils

# ╔═╡ 1a714fba-2a8e-11ec-182f-8f85cc17b02a
begin
	using Pkg
	Pkg.activate("../docs")
	
	using MeshArrays, Plots, PlutoUI
	
	p=dirname(pathof(MeshArrays))
	include(joinpath(p,"../examples/Demos.jl"))
	include(joinpath(p,"../examples/Plots.jl"))
	
	"All done with packages."
end

# ╔═╡ 0358ae4b-f941-4254-aa10-6c37bb09c464
md"""# Get Started

To create your first `MeshArray`, open the `Julia` REPL and type:

```
using MeshArrays
tmp=MeshArray(randn(20,10))
```
"""

# ╔═╡ b1f3ac4b-5f75-40fa-b125-186556ece946
begin
	tmp=MeshArray(randn(20,10))
	with_terminal() do
		show(tmp)
	end
end

# ╔═╡ a0a7637f-f115-4257-a140-06dc48dd5fba
md"""## Tutorial

In this tutorial, the same example is repeated three times -- for three different grid configurations commonly used in numerical models (see [Earth Grids](@ref)). The example proceeds as follows in each grid case:

1. generate grid configuration
1. initialize map with random noise
1. apply smoothing across whole domain
1. plot results for the subdomain arrays

Step 3 illustrates how `MeshArrays.jl` computes partial derivatives between neighboring subdomains. The `MeshArrays.smooth()` method is indeed based on lateral diffusion, a transport process, which is integrated through time over the global domain. Other processes like e.g. advection by ocean currents would work similarly.
"""

# ╔═╡ 54030bff-45f7-43b1-b19f-7c0cd05a02ae
md"""

Grid scale noise           |  Smoothed noise
:------------------------------:|:---------------------------------:
![raw](https://user-images.githubusercontent.com/20276764/118325229-2d883d80-b4d1-11eb-953b-ddbb11bcfe1b.png)  |  ![smooth](https://user-images.githubusercontent.com/20276764/118325093-f31ea080-b4d0-11eb-8c6e-8cd0cc2cc255.png)

### 1. Doubly Periodic Domain

Let's start with a doubly periodic domain split into `nF=16` subdomains. Each subdomain corresponds to an array of `nP=20` by `nQ=20` grid points. 
The `UnitGrid()` function is used to generate such a grid configuration.
"""

# ╔═╡ 1fc99fca-3440-4765-b11d-e10eb560aff7
begin
	(nP,nQ,nF)=(20,20,16)
	facesSize=fill((nP,nQ),nF)
	ioSize=[nP nQ*nF]
	
	γ=gcmgrid("","PeriodicDomain",nF,facesSize, ioSize,Float32, read, write)
	Γ=UnitGrid(γ;option="full")
	
	"Done with basic grid"
end

# ╔═╡ cef83984-80dc-4196-8edd-314885fa9448
md"""Then we do steps 2 (`zin`) and 3 (`zout`) as follows. 

The `heatmap` function allows you to visualize that `zout` is indeed smoother (and therefore muted) than is `zin`.
"""

# ╔═╡ 2573fda7-66cb-484f-9d06-a9bb85a26a9e
begin
	#initialize 2D field of random numbers
	tmp1=randn(Float32,Tuple(γ.ioSize))
	zin =γ.read(tmp1,MeshArray(γ,Float32))
	
	#smoothing length scales in x, y directions
	Lx=3*Γ.DXC; Ly=3*Γ.DYC
	
	#apply smoother
	zout=smooth(zin,Lx,Ly,Γ)
	
	#visuzalize
	heatmap(zout,clims=(-0.25,0.25),tickfont = (4, :black))
end

# ╔═╡ 2ee0a449-fd98-43e1-90b4-b4e15f237d67
begin
	let
		(nP,nQ,nF)=(32,32,6)
		facesSize=fill((nP,nQ),nF)
		ioSize=[nP nQ*nF]

		γ=gcmgrid("","CubeSphere",nF,facesSize, ioSize,Float32, read, write)
		Γ=UnitGrid(γ;option="full")
		D=demo2(Γ)
	end
	
	md"""
	### 2. Cube Sphere
	
	In this cell we instead use a grid that has 6 subdomains, 100x100 points each, covering the six faces of a cube. We note that this _cube sphere_ topology involves connections between subdomain that are slightly more complicated than in the first example. And we call `demo2()` which combines steps 2 and 3. 
	"""
end

# ╔═╡ d37ae153-8646-44ae-abd9-7e0e7d479275
begin
	let
		γ=GridSpec("LatLonCap",MeshArrays.GRID_LLC90)
		Γ=GridLoad(γ;option="full")
		D=demo2(Γ)
	end
	
	md"""### 3. LLC90 grid
	
	The [Lat-Lon-Cap grid](http://www.geosci-model-dev.net/8/3071/2015/) (or LLC) is a global ocean model grid which is widely used in the [MITgcm user community](https://mitgcm.readthedocs.io/en/latest/). It has 5 uneven subdomains, variable grid spacing, and continents [(Forget et al 2015)](http://www.geosci-model-dev.net/8/3071/2015/). LLC90's resolution is one degree albeit with modications in the Arctic and along the Equator.
	
	In this case, the grid variables are read from files found in the `MeshArrays.GRID_LLC90` folder by the `GridLoad` function. The `GridSpec` function provides the corresponding domain sizes for this commonly used, `LLC90`, grid. 
	"""
end

# ╔═╡ Cell order:
# ╟─0358ae4b-f941-4254-aa10-6c37bb09c464
# ╟─b1f3ac4b-5f75-40fa-b125-186556ece946
# ╟─1a714fba-2a8e-11ec-182f-8f85cc17b02a
# ╟─a0a7637f-f115-4257-a140-06dc48dd5fba
# ╟─54030bff-45f7-43b1-b19f-7c0cd05a02ae
# ╠═1fc99fca-3440-4765-b11d-e10eb560aff7
# ╟─cef83984-80dc-4196-8edd-314885fa9448
# ╠═2573fda7-66cb-484f-9d06-a9bb85a26a9e
# ╟─2ee0a449-fd98-43e1-90b4-b4e15f237d67
# ╟─d37ae153-8646-44ae-abd9-7e0e7d479275
