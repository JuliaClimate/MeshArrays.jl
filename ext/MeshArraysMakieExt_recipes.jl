
"""
    scatter(XC::AbstractMeshArray,YC::AbstractMeshArray;color=:black,colorrange=[],colorbar=true,title="",kwargs...)

```
scatter(Γ.XC,Γ.YC,color=:black)
MS=log10.(Γ.RAC)*μ
scatter(Γ.XC,Γ.YC,color=MS)
```
"""
function scatter(XC::AbstractMeshArray,YC::AbstractMeshArray;
	color=:black,colorrange=[],colorbar=true,title="",kwargs...)

	if isa(color,MeshArray)
		γ=color.grid
		DD=γ.write(color)	
		!isempty(colorrange) ? cr=colorrange : cr=(nanmin(DD),nanmax(DD))
		cr[1]==cr[2] ? cr=(cr[1]-eps(),cr[2]+eps()) : nothing
	else
		cr=[]
	end

	fig = Figure(size = (900,600), backgroundcolor = :grey95)
	ax = Axis(fig[1,1],xlabel="longitude",ylabel="latitude",title=title)
	scatter!(ax,XC::AbstractMeshArray,YC::AbstractMeshArray; 
	color=color, colorrange=colorrange, colorbar=colorbar, kwargs...)

	fig
end

"""
    scatter!(ax,XC::AbstractMeshArray,YC::AbstractMeshArray;color=:black,colorrange=[],colormap=:viridis)

```
fig=heatmap(Γ.Depth,interpolation=λ)
scatter!(current_axis(),Γ.XC,Γ.YC,color=:red)
fig
```
"""	
function scatter!(ax,XC::AbstractMeshArray,YC::AbstractMeshArray;
	color=:black,colorrange=[],colorbar=true,colormap=:veridis,kwargs...)

	if isa(color,MeshArray)&&isempty(colorrange)
		γ=color.grid
		DD=γ.write(color)	
		!isempty(colorrange) ? cr=colorrange : cr=(nanmin(DD),nanmax(DD))
		cr[1]==cr[2] ? cr=(cr[1]-eps(),cr[2]+eps()) : nothing
	else
		cr=colorrange
	end

	for ff in eachindex(XC)
		if isa(color,Symbol)
			scatter!(ax,XC[ff][:],YC[ff][:];color=color,kwargs...)
		else
			kk=findall((!isnan).(color[ff]))
			scatter!(ax,XC[ff][kk],YC[ff][kk];color=color[ff][kk],
				colorrange = cr,colormap=colormap, kwargs...)
		end
	end

	fig=current_figure()
	colorbar&&isa(color,MeshArray) ? Colorbar(fig[1,2], colorrange=cr, colormap=colormap, height = Relative(0.65)) : nothing

end

"""
    heatmap(MS::AbstractMeshArray; interpolation=nothing,globalmap=false,x=nothing,y=nothing,colorbar=true,title="",kwargs...)

Represent a `MeshArray` as a `heatmap`, or several, depending on keyword parameter choices. 

Additional keyword arguments will passed along to the `Makie.heatmap` call.

```
heatmap(MS) #will display tile by tile
heatmap(MS,interpolation=λ) #will interpolate on the fly
heatmap(MS,interpolation=λ,title="ocean depth") #same but w title
heatmap(MS,x=lon,y=lat) #only for simple domains; will show MS[1]
```
"""
heatmap(MS::AbstractMeshArray;interpolation=nothing,globalmap=false,x=nothing,y=nothing,colorbar=true,title="",kwargs...)=begin
	if (!isnothing(interpolation))||globalmap||(!isnothing(x))
		f=Figure()
		ax=Axis(f[1,1],title=title)
		hm1=heatmap!(ax,MS::AbstractMeshArray;
				interpolation=interpolation,globalmap=globalmap,x=x,y=y,
				kwargs...)
		colorbar ? Colorbar(f[1,2], hm1, height = Relative(0.65)) : nothing
		f
	else
		heatmap_tiled(MS;colorbar=colorbar,title=title,kwargs...)
	end
end

function heatmap!(ax::Axis,MS::AbstractMeshArray;interpolation=nothing,globalmap=false,x=nothing,y=nothing,kwargs...)	
	if !isnothing(interpolation)
		heatmap_interpolation!(ax,MS,interpolation;kwargs...)
	elseif globalmap
		heatmap_globalmap!(ax,MS;kwargs...)
	elseif !isnothing(x)
		heatmap_xy!(ax,MS,x,y;kwargs...)
	else
		heatmap_tiled(MS;kwargs...)
	end
end

function heatmap_globalmap!(ax,MS::AbstractMeshArray;kwargs...)
	γ=MS.grid
	DD=γ.write(MS)	
#	!isempty(colorrange) ? cr=colorrange : cr=(nanmin(DD),nanmax(DD))
	hm1=heatmap!(ax,DD;kwargs...)
end

function heatmap_globalmap(MS::AbstractMeshArray;title="",colorbar=true,kwargs...)
    fig = Figure(size = (900,900), backgroundcolor = :grey95)
    ax = Axis(fig[1,1],xlabel="i index",ylabel="j index",title=title)
	hm1=heatmap_globalmap!(ax,MS;kwargs...)
    colorbar ? Colorbar(fig[1,2], hm1, height = Relative(0.65)) : nothing
    fig
end

heatmap_xy!(ax,MS::AbstractMeshArray,x::Union{UnitRange,Array},y::Union{UnitRange,Array};kwargs...) = heatmap!(ax,x,y,MS[1];kwargs...)
#heatmap_xy!(ax,MS::AbstractMeshArray,x::Union{UnitRange,Array},y::Union{UnitRange,Array};kwargs...) = surface!(ax,x,y,0*x;color=MS[1],shading=NoShading,	kwargs...)

function heatmap_xy(MS::AbstractMeshArray,x::Union{UnitRange,Array},y::Union{UnitRange,Array};title="",colorbar=true,kwargs...)
    fig = Figure(size = (900,400), backgroundcolor = :grey95)
    ax = Axis(fig[1,1],xlabel="longitude",ylabel="latitude",title=title)
	hm1=heatmap_xy!(ax,MS,x,y;kwargs...)
    colorbar ? Colorbar(fig[1,2], hm1, height = Relative(0.65)) : nothing
    fig
end

function heatmap_interpolation!(ax,MS::AbstractMeshArray,λ::NamedTuple;kwargs...)
    DD=Interpolate(MS,λ.f,λ.i,λ.j,λ.w)
	DD=reshape(DD,size(λ.lon))
	hm1=heatmap!(ax,λ.lon[:,1],λ.lat[1,:],DD;kwargs...)
end

function heatmap_interpolation(MS::AbstractMeshArray,λ::NamedTuple;title="",colorbar=true,kwargs...)
    fig = Figure(size = (900,400), backgroundcolor = :grey95)
    ax = Axis(fig[1,1],xlabel="longitude",ylabel="latitude",title=title)
	hm1=heatmap_interpolation!(ax,MS,λ;kwargs...)
    colorbar ? Colorbar(fig[1,2], hm1, height = Relative(0.65)) : nothing
    fig
end

function heatmap_tiled(MS::AbstractMeshArray;title="",
	colorbar=true,colorrange=[],colormap=:viridis,kwargs...)
	fig = Figure(size = (900,900), backgroundcolor = :grey95)
	nf=length(MS.fSize)

	test1=in(MS.grid.class,["PeriodicChannel","PeriodicDomain"])
	test2=in(MS.grid.class,["CubeSphere"])
	test3=in(MS.grid.class,["LatLonCap"])
	if test1
		(np,nq)=Int.(MS.grid.ioSize[:]./MS.fSize[1][:])
		ii=[i for j in 1:np, i in nq:-1:1]
		jj=[j for j in 1:np, i in nq:-1:1]
		kk=(1:nq,np+1)
	elseif test2
		(np,nq)=Int.(MS.grid.ioSize[:]./MS.fSize[3][:])
		ii=[3, 3, 2, 2, 1, 1]
		jj=[1, 2, 2, 3, 3, 4]
		kk=(1:3,5)
	elseif test3
		(np,nq)=Int.(MS.grid.ioSize[:]./MS.fSize[3][:])
		ii=[3, 3, 2, 2, 1]
		jj=[1, 2, 2, 3, 3]
		kk=(1:3,4)
	else
		error("unknown grid class")
	end
	
	!isempty(colorrange) ? cr=colorrange : cr=(nanmin(write(MS)),nanmax(write(MS)))
	cr[1]==cr[2] ? cr=(cr[1]-eps(),cr[2]+eps()) : nothing

	for f in 1:nf
		nf > 1 ? ax = Axis(fig[ii[f],jj[f]], title=title*" face $(f)") : ax = Axis(fig[1,1], title=title)

		s=MS.fSize[f]		
		x=collect(0.5:s[1]-0.5)
		y=collect(0.5:s[2]-0.5)
		z=MS[f]

		hm1=heatmap!(ax,x,y,z;colorrange=cr,colormap=colormap,kwargs...)
	end

	colorbar ? Colorbar(fig[kk...], limits=cr, colormap=colormap, height = Relative(0.65)) : nothing
	
	fig
end

## convert polygon data to use with Makie

"""
	pol_to_Makie(tmp2::Vector)

Convert output of `read_json` or `read_shp`` to a vector of `LineString`.
"""
function pol_to_Makie(pol) #::polyarray)
	if isa(pol,polyarray)
		tmp2=[[x.geometry] for x in pol.data]
	else
		tmp2=[GI.coordinates(a.geometry) for a in pol]
	end	
#	tmp2=pol.data[1].geometry
	tmp22=Vector{Point2{Float64}}[]
	for l1 in tmp2
		if isa(l1[1][1][1],Number)
			push!(tmp22,geo2basic(tuple2vec.(l1[1])))
		else
			for l2 in l1
				push!(tmp22,geo2basic(tuple2vec.(l2[1])))
			end
		end
	end
	
	LineString.(tmp22)
end

tuple2vec(x)=[y for y in x]

## convert to GeometryBasics

to_point2(a::Vector{<: T}) where T = Point2{T}(a[1], a[2])
to_point2(a::AbstractVector{T}) where T <: Number = Point2{T}(a[1], a[2])

"""
	function geo2basic(vector::AbstractVector{<:AbstractVector})

Source : @SimonDanisch , https://github.com/MakieOrg/GeoMakie.jl/pull/125
"""
function geo2basic(vector::AbstractVector{<:AbstractVector})
	if isempty(vector)
		return Point{2, Float64}[]
	else
		# GeoJSON strips the eltype so we need to inspect the elements
		x = first(vector)
		if x isa AbstractVector && length(x) == 2 && x[1] isa Real
			return to_point2.(vector)
		elseif x isa AbstractVector && eltype(x) <: AbstractVector
			linestrings = map(x-> to_point2.(x), vector)
			return Makie.Polygon(linestrings[1], linestrings[2:end])
		else
			error("Unsupported eltype: $(x)")
		end
	end
end

## plot methods

plot(x::AbstractMeshArray; kwargs...) = begin
	if ndims(x) == 1
		heatmap(x; kwargs...)
	else
		heatmap(x[:,1]; kwargs...)
	end
end

function plot(x::Union{gridpath,Vector{gridpath}}; kwargs...)
	fig=Figure(); ax=Axis(fig[1,1],limits=(-180.0,180.0,-90.0,90.0))
	if isa(x,gridpath)
		plot!(x; kwargs...)
	else
		[plot!(y; kwargs...) for y in x]
	end
	fig
end

function plot!(x::gridpath; kwargs...)
	np=size(x.C,1)
	lon=zeros(np)
	lat=zeros(np)
	for p in 1:np
		f=x.C[p,1]; i=x.C[p,2]; j=x.C[p,3]; q=x.C[p,4]
		lon[p]=x.grid.XC[f][i,j]
		lat[p]=x.grid.YC[f][i,j]
	end
	scatter!(lon,lat; kwargs...)
end

"""
    polygons_plot(pols; color=:black, outer_edge=false)

```
polygons_plot(pols,color=Depth,outer_edge=true)
polygons_plot(pols,color=Depth)
polygons_plot(pols)
```
"""
function polygons_plot(pols; color=:black, outer_edge=false)
	pols_Makie=vcat([[Point2.(GI.coordinates(q)[1]) for q in p][:] for p in pols]...)
	if isa(color,MeshArray)
		cols_Makie=vcat([D[:] for D in color]...)
	else
		cols_Makie=color
	end
	if (!outer_edge)&&isa(color,MeshArray)
		poly(pols_Makie, color = cols_Makie)
	else
		poly(pols_Makie,color=:white,strokecolor=cols_Makie,strokewidth=2)
	end
end

lines(pa::polyarray,stuff...;kwargs...) = 
	lines(pol_to_Makie(pa),stuff...;kwargs...)
lines!(pa::polyarray,stuff...;kwargs...) = 
	lines!(pol_to_Makie(pa),stuff...;kwargs...)
plot(pa::polyarray,stuff...;kwargs...) = 
	plot(pol_to_Makie(pa),stuff...;kwargs...)
plot!(pa::polyarray,stuff...;kwargs...) = 
	plot!(pol_to_Makie(pa),stuff...;kwargs...)
