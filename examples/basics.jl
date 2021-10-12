### A Pluto.jl notebook ###
# v0.16.1

using Markdown
using InteractiveUtils

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
	## Julia Packages Being Used
	
	"""
end

# ╔═╡ 0dbf6e21-2dd7-497f-9b0b-f67707f2e7e0
module myDemos
	using Plots, MeshArrays
	p=dirname(pathof(MeshArrays))
	include(joinpath(p,"../examples/Demos.jl"))
end

# ╔═╡ 0358ae4b-f941-4254-aa10-6c37bb09c464
md"""# Get Started

!!! tip
    You can download this notebook [here](https://raw.githubusercontent.com/JuliaClimate/MeshArrays.jl/master/examples/basics.jl) and experience it interactively using [Pluto.jl](https://plutojl.org). This can be a good way to get familiar with software such as this one.

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

In this tutorial, the same example is repeated three times -- for three different grid configurations commonly used in numerical models (see [Grids](https://juliaclimate.github.io/MeshArrays.jl/dev/start/#Grids)). The example proceeds as follows in each grid case:

1. generate grid configuration
1. initialize map with random noise
1. apply smoothing across whole domain
1. plot results for the subdomain arrays

Step 3 illustrates how `MeshArrays.jl` computes partial derivatives between neighboring subdomains. The `MeshArrays.smooth()` method is indeed based on lateral diffusion, a transport process, which is integrated through time over the global domain. Other processes like e.g. advection by ocean currents would work similarly.
"""

# ╔═╡ 54030bff-45f7-43b1-b19f-7c0cd05a02ae
md"""### 1. Doubly Periodic Domain

Let's start with a doubly periodic domain split into `nF=16` subdomains. Each subdomain corresponds to an array of `nP=20` by `nQ=20` grid points. 
The `UnitGrid()` function is used to generate such a grid configuration.
"""

# ╔═╡ 1fc99fca-3440-4765-b11d-e10eb560aff7
begin
	(nP,nQ,nF)=(20,20,16)
	facesSize=fill((nP,nQ),nF)
	ioSize=[nP nQ*nF]
	#ioSize=[nP*4 nQ*4]
	
	γ=gcmgrid("","PeriodicDomain",nF,facesSize, ioSize,Float32, read, write)
	Γ=UnitGrid(γ;option="full")
	
	"Done with basic grid"
end

# ╔═╡ cef83984-80dc-4196-8edd-314885fa9448
md"""Then we do steps 2 (`zin`) and 3 (`zout`) as follows. 

!!! note
    The `heatmap` display should reflect that `zout` is smoother (and therefore muted) than is `zin`.
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
	
	#visualize
	x=write(Γ.XC)[:,1]
	y=write(Γ.YC)[1,:]
	
	fig = Mkie.Figure(resolution = (900,900), backgroundcolor = :grey95)

	ax1 = Mkie.Axis(fig[1,1])
	hm1=Mkie.heatmap!(ax1,x,y,write(zin),clims=(-0.25,0.25),tickfont = (4, :black))
	ax2 = Mkie.Axis(fig[1,2])
	hm2=Mkie.heatmap!(ax2,x,y,write(zout),clims=(-0.25,0.25),tickfont = (4, :black))
	Mkie.Colorbar(fig[1,3], hm1, height = Mkie.Relative(0.65))
	
	fig
end

# ╔═╡ 2ee0a449-fd98-43e1-90b4-b4e15f237d67
md"""
### 2. Cube Sphere

In this cell we instead use a grid that has 6 subdomains, 100x100 points each, covering the six faces of a cube. We note that this _cube sphere_ topology involves connections between subdomain that are slightly more complicated than in the first example. And we call `demo2()` which combines steps 2 and 3. 
"""

# ╔═╡ dacc2261-1d75-4d68-b82c-174e4fae3631
let
	(nP,nQ,nF)=(32,32,6)
	facesSize=fill((nP,nQ),nF)
	ioSize=[nP nQ*nF]

	γ=gcmgrid("","CubeSphere",nF,facesSize, ioSize,Float32, read, write)
	Γ=UnitGrid(γ;option="full")
	D=myDemos.demo2(Γ)
	PlutoUI.with_terminal() do
		show(D[1])
	end
end

# ╔═╡ d37ae153-8646-44ae-abd9-7e0e7d479275
md"""### 3. LLC90 grid

The [Lat-Lon-Cap grid](http://www.geosci-model-dev.net/8/3071/2015/) (or LLC) is a global ocean model grid which is widely used in the [MITgcm user community](https://mitgcm.readthedocs.io/en/latest/). It has 5 uneven subdomains, variable grid spacing, and continents [(Forget et al 2015)](http://www.geosci-model-dev.net/8/3071/2015/). LLC90's resolution is one degree albeit with modications in the Arctic and along the Equator.

In this case, the grid variables are read from files found in the `MeshArrays.GRID_LLC90` folder by the `GridLoad` function. The `GridSpec` function provides the corresponding domain sizes for this commonly used, `LLC90`, grid. 
"""

# ╔═╡ 23198b8d-6bf0-4af9-91d0-7886c9574f81
let
	γ=GridSpec("LatLonCap",MeshArrays.GRID_LLC90)
	Γ=GridLoad(γ;option="full")
	D=myDemos.demo2(Γ)
	PlutoUI.with_terminal() do
		show(D[1])
	end
end

# ╔═╡ 58f95665-9687-4b4f-af99-1239818f71a3
md"""### More"""

# ╔═╡ 2cba78ac-d69b-4c47-96d0-7b784c53246b
let
	GridVariables=GridOfOnes("PeriodicDomain",16,20;option="full")
	GridVariables[1].ioSize
	x=write(GridVariables[2].XC)[:,1]
	y=write(GridVariables[2].YC)[1,:]
	Mkie.heatmap(x,y,write(GridVariables[2].YC))
end

# ╔═╡ Cell order:
# ╟─0358ae4b-f941-4254-aa10-6c37bb09c464
# ╟─b1f3ac4b-5f75-40fa-b125-186556ece946
# ╠═1a714fba-2a8e-11ec-182f-8f85cc17b02a
# ╟─0dbf6e21-2dd7-497f-9b0b-f67707f2e7e0
# ╟─a0a7637f-f115-4257-a140-06dc48dd5fba
# ╟─54030bff-45f7-43b1-b19f-7c0cd05a02ae
# ╠═1fc99fca-3440-4765-b11d-e10eb560aff7
# ╟─cef83984-80dc-4196-8edd-314885fa9448
# ╠═2573fda7-66cb-484f-9d06-a9bb85a26a9e
# ╟─2ee0a449-fd98-43e1-90b4-b4e15f237d67
# ╟─dacc2261-1d75-4d68-b82c-174e4fae3631
# ╟─d37ae153-8646-44ae-abd9-7e0e7d479275
# ╠═23198b8d-6bf0-4af9-91d0-7886c9574f81
# ╟─58f95665-9687-4b4f-af99-1239818f71a3
# ╠═2cba78ac-d69b-4c47-96d0-7b784c53246b
