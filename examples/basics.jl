### A Pluto.jl notebook ###
# v0.16.1

using Markdown
using InteractiveUtils

# ╔═╡ 1a714fba-2a8e-11ec-182f-8f85cc17b02a
begin
	using Pkg
	Pkg.activate()
	using MeshArrays, Plots
end

# ╔═╡ 0358ae4b-f941-4254-aa10-6c37bb09c464
md"""# Basics Of MeshArrays.jl"""

# ╔═╡ b1f3ac4b-5f75-40fa-b125-186556ece946
MeshArray(randn(20,10))

# ╔═╡ b72cb2cb-f4c6-4990-9279-710f04881e23
begin
	(nP,nQ,nF)=(20,20,16)
	facesSize=fill((nP,nQ),nF)
	ioSize=[nP nQ*nF]
end

# ╔═╡ ad88ce66-340f-453e-bfb2-14e195d8f82f
begin
	γ=gcmgrid("","PeriodicDomain",
				nF,facesSize, ioSize,
				Float32, read, write)
	Γ=UnitGrid(γ;option="full")
end

# ╔═╡ 1bb800a3-6340-4252-81e4-7fca7ab6d454
begin
	#initialize 2D field of random numbers
	tmp1=randn(Float32,Tuple(γ.ioSize))
	zin =γ.read(tmp1,MeshArray(γ,Float32))
	
	#smoothing length scales in x, y directions
	Lx=3*Γ.DXC; Ly=3*Γ.DYC
	
	#apply smoother
	zout=smooth(zin,Lx,Ly,Γ)
end

# ╔═╡ 19bcbbb6-4055-4725-9202-6cd921509fec
begin
	p=dirname(pathof(MeshArrays))
	include(joinpath(p,"../examples/Plots.jl"))
	heatmap(zout,clims=(-0.25,0.25),tickfont = (4, :black))
end

# ╔═╡ Cell order:
# ╟─0358ae4b-f941-4254-aa10-6c37bb09c464
# ╠═1a714fba-2a8e-11ec-182f-8f85cc17b02a
# ╠═b1f3ac4b-5f75-40fa-b125-186556ece946
# ╠═b72cb2cb-f4c6-4990-9279-710f04881e23
# ╠═ad88ce66-340f-453e-bfb2-14e195d8f82f
# ╠═1bb800a3-6340-4252-81e4-7fca7ab6d454
# ╠═19bcbbb6-4055-4725-9202-6cd921509fec
