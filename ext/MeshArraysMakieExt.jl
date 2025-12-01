
module MeshArraysMakieExt

using MeshArrays, Makie

import MeshArrays: GI
import MeshArrays: land_mask
import MeshArrays: plot_examples
import MeshArrays: ProjAxis
import MeshArrays: grid_lines!

import MeshArrays: gridpath

import Makie: plot, plot!, heatmap, scatter, scatter!, surface!
import Makie: lines!, heatmap!, contour!, contourf!

LineString=Makie.LineString
Observable=Makie.Observable
GeometryBasics=Makie.GeometryBasics

function plot_examples(ID=Symbol,stuff...;kwargs...)
	if ID==:interpolation_demo
		interpolation_demo(stuff...)
	elseif ID==:projmap
		projmap(stuff...)
	elseif ID==:simple_heatmap
		simple_heatmap(stuff...)
	elseif ID==:smoothing_demo
		smoothing_demo(stuff...)
	elseif ID==:northward_transport
		northward_transport(stuff...)
	elseif ID==:meriodional_overturning
		meriodional_overturning(stuff...)
	elseif ID==:gradient_EN
		gradient_EN(stuff...)
	elseif ID==:gradient_xy
		gradient_xy(stuff...)
	elseif ID==:basemap
		basemap(stuff...)
	elseif ID==:baseproj
		baseproj(stuff...;kwargs...)
	else
		println("unknown plot ID")
	end
end

##

"""
    scatter(XC::MeshArray,YC::MeshArray;color=:black,colorrange=[],colorbar=true,title="",kwargs...)

```
scatter(Γ.XC,Γ.YC,color=:black)
MS=log10.(Γ.RAC)*μ
scatter(Γ.XC,Γ.YC,color=MS)
```
"""
function scatter(XC::MeshArray,YC::MeshArray;
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
	scatter!(ax,XC::MeshArray,YC::MeshArray; 
	color=color, colorrange=colorrange, colorbar=colorbar, kwargs...)

	fig
end

"""
    scatter!(ax,XC::MeshArray,YC::MeshArray;color=:black,colorrange=[],colormap=:viridis)

```
fig=heatmap(Γ.Depth,interpolation=λ)
scatter!(current_axis(),Γ.XC,Γ.YC,color=:red)
fig
```
"""	
function scatter!(ax,XC::MeshArray,YC::MeshArray;
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
    heatmap(MS::MeshArray; interpolation=nothing,globalmap=false,x=nothing,y=nothing,colorbar=true,title="",kwargs...)

Represent a `MeshArray` as a `heatmap`, or several, depending on keyword parameter choices. 

Additional keyword arguments will passed along to the `Makie.heatmap` call.

```
heatmap(MS) #will display tile by tile
heatmap(MS,interpolation=λ) #will interpolate on the fly
heatmap(MS,interpolation=λ,title="ocean depth") #same but w title
heatmap(MS,x=lon,y=lat) #only for simple domains; will show MS[1]
```
"""
heatmap(MS::MeshArray;interpolation=nothing,globalmap=false,x=nothing,y=nothing,colorbar=true,title="",kwargs...)=begin
	if (!isnothing(interpolation))||globalmap||(!isnothing(x))
		f=Figure()
		ax=Axis(f[1,1],title=title)
		hm1=heatmap!(ax,MS::MeshArray;
				interpolation=interpolation,globalmap=globalmap,x=x,y=y,
				kwargs...)
		colorbar ? Colorbar(f[1,2], hm1, height = Relative(0.65)) : nothing
		f
	else
		heatmap_tiled(MS;colorbar=colorbar,title=title,kwargs...)
	end
end

function heatmap!(ax::Axis,MS::MeshArray;interpolation=nothing,globalmap=false,x=nothing,y=nothing,kwargs...)	
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

function heatmap_globalmap!(ax,MS::MeshArray;kwargs...)
	γ=MS.grid
	DD=γ.write(MS)	
#	!isempty(colorrange) ? cr=colorrange : cr=(nanmin(DD),nanmax(DD))
	hm1=heatmap!(ax,DD;kwargs...)
end

function heatmap_globalmap(MS::MeshArray;title="",colorbar=true,kwargs...)
    fig = Figure(size = (900,900), backgroundcolor = :grey95)
    ax = Axis(fig[1,1],xlabel="i index",ylabel="j index",title=title)
	hm1=heatmap_globalmap!(ax,MS;kwargs...)
    colorbar ? Colorbar(fig[1,2], hm1, height = Relative(0.65)) : nothing
    fig
end

heatmap_xy!(ax,MS::MeshArray,x::Union{UnitRange,Array},y::Union{UnitRange,Array};kwargs...) = heatmap!(ax,x,y,MS[1];kwargs...)
#heatmap_xy!(ax,MS::MeshArray,x::Union{UnitRange,Array},y::Union{UnitRange,Array};kwargs...) = surface!(ax,x,y,0*x;color=MS[1],shading=NoShading,	kwargs...)

function heatmap_xy(MS::MeshArray,x::Union{UnitRange,Array},y::Union{UnitRange,Array};title="",colorbar=true,kwargs...)
    fig = Figure(size = (900,400), backgroundcolor = :grey95)
    ax = Axis(fig[1,1],xlabel="longitude",ylabel="latitude",title=title)
	hm1=heatmap_xy!(ax,MS,x,y;kwargs...)
    colorbar ? Colorbar(fig[1,2], hm1, height = Relative(0.65)) : nothing
    fig
end

function heatmap_interpolation!(ax,MS::MeshArray,λ::NamedTuple;kwargs...)
    DD=Interpolate(MS,λ.f,λ.i,λ.j,λ.w)
	DD=reshape(DD,size(λ.lon))
	hm1=heatmap!(ax,λ.lon[:,1],λ.lat[1,:],DD;kwargs...)
end

function heatmap_interpolation(MS::MeshArray,λ::NamedTuple;title="",colorbar=true,kwargs...)
    fig = Figure(size = (900,400), backgroundcolor = :grey95)
    ax = Axis(fig[1,1],xlabel="longitude",ylabel="latitude",title=title)
	hm1=heatmap_interpolation!(ax,MS,λ;kwargs...)
    colorbar ? Colorbar(fig[1,2], hm1, height = Relative(0.65)) : nothing
    fig
end

function heatmap_tiled(MS::MeshArray;title="",
	colorbar=true,colorrange=[],colormap=:viridis,kwargs...)
	fig = Figure(size = (900,900), backgroundcolor = :grey95)
	nf=length(MS.fSize)
	nn=Int(ceil(nf/2))
	ii=[i for j in 1:2, i in 1:nn]
	jj=[j for j in 1:2, i in 1:nn]
	
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

	nf > 1 ? jj=(1:nn,3) : jj=(1,2)
	colorbar ? Colorbar(fig[jj...], limits=cr, colormap=colormap, height = Relative(0.65)) : nothing
	
	fig
end

##

function northward_transport(MT)
	x=vec(-89.0:89.0)
	fig1 = Figure(size = (900,400),markersize=0.1)
	ax1 = Axis(fig1[1,1],xlabel="latitude", ylabel="Transport (in Sv)", 
		title="Northward Volume Transport (in Sv)")
	hm1=lines!(x,1e-6*MT,label="ECCO estimate")
	fig1
end

function meriodional_overturning(Γ,ov)
	x=vec(-89.0:89.0); y=reverse(vec(Γ.RF[1:end-1])); #coordinate variables
	z=reverse(ov,dims=2); z[z.==0.0].=NaN

	fig1 = Figure(size = (900,400),markersize=0.1)
	ax1 = Axis(fig1[1,1], title="Meridional Overturning Streamfunction (in Sv)")
	hm1=contourf!(ax1,x,y,1e-6*z,levels=(-40.0:5.0:40.0))
	Colorbar(fig1[1,2], hm1, height = Relative(0.65))
	fig1
	#savefig("MOC_mean.png")
end

##

function gradient_EN(λ,dDdx,dDdy)
	fig1 = Figure(size = (900,600),markersize=0.1)
	ax1 = Axis(fig1[1,1], title="Gradient of scalar potential in Eastward direction (in 1/s)")
	hm1=heatmap_interpolation!(ax1,dDdx,λ,colorrange=(-0.1,0.1))
	ax1 = Axis(fig1[2,1], title="Gradient of scalar potential in Northward direction (in 1/s)")
	hm1=heatmap_interpolation!(ax1,dDdy,λ,colorrange=(-0.1,0.1))
	Colorbar(fig1[1:2,2], hm1, height = Relative(0.65))
	fig1
end

function gradient_xy(λ,dDdx,dDdy)
	fig = Figure(size = (900,600), backgroundcolor = :grey95)
	ax = Axis(fig[1,1], title="x-direction velocity (in m/s)",xlabel="longitude",ylabel="latitude")
	hm1=heatmap_interpolation!(ax,dDdx,λ,colorrange=(-0.2,0.2))
	ax = Axis(fig[2,1], title="y-direction velocity (in m/s)",xlabel="longitude",ylabel="latitude")
	hm1=heatmap_interpolation!(ax,dDdy,λ,colorrange=(-0.2,0.2))
	Colorbar(fig[1:2,2], hm1, height = Relative(0.65))
	fig
end

##

function simple_heatmap(dat)	
	lons = dat.lon[:,1]
	lats = dat.lat[1,:]
	field = dat.var

	fig = Figure(size = (1200,800), fontsize = 22)
	ax = Axis(fig[1,1])
	hm1 = heatmap!(ax, lons, lats, field, colorrange=dat.meta.colorrange, colormap=dat.meta.cmap)
	Colorbar(fig[1,2], hm1, height = Relative(0.65))

	fig
end

"""
    projmap(data,lon0=lon0,proj=proj)

Inputs:

- data is a `NamedTuple`
- lon0 is the central longitude
- proj is a `Proj.Transformation`

Use examples:

- `MeshArrays.jl/examples/geography.jl`
- `ClimateModels.jl/examples/IPCC.jl`
"""
function projmap(data,lon0=lon0,proj=proj)
	
	f = Figure()
    ax = f[1, 1] = Axis(f, aspect = DataAspect(), title = data.meta.ttl)
	pr_ax=MeshArrays.ProjAxis(ax; proj=proj,lon0=lon0)

    surf = surface!(pr_ax,data.lon,data.lat,0*data.lat.-1; color=data.var, 
	colorrange=data.meta.colorrange, colormap=data.meta.cmap,
        shading = NoShading)

	haskey(data.meta,:gridcolor) ? grid_lines!(pr_ax,color=data.meta.gridcolor,linewidth=0.5) : nothing
	haskey(data,:polygons) ? lines!(pr_ax; polygons=data.polygons,color=:black,linewidth=0.5) : nothing

    Colorbar(f[1,2], surf, height = Relative(0.5))

	f
end

function calc_shift(lo,lo0)
	lon_size=length(lo)
	lon_range=(lo[1,1],lo[end,1])
	Int(lon_size/360)*(360+lo0)
end

##

function interpolation_demo(Γ)
	lon=collect(0.1:0.5:2.1); lat=collect(0.1:0.5:2.1);
	#lon=[66.75]; lat=[-69.75];

	(f,i,j,w,_,_,_)=InterpolationFactors(Γ,vec(lon),vec(lat))
	Depth_int=Interpolate(Γ.Depth,f,i,j,w)

	#find nearest neighbor (`MeshArray` & `set`)
	(f,i,j,c)=knn(Γ.XC,Γ.YC,lon,lat)
	[write(Γ.XC)[c] write(Γ.YC)[c]]

	γ=Γ.XC.grid
	
	#define subdomain tiles (`MeshArray` -> `tiles`)
	ni=30; nj=30;
	τ=Tiles(γ,ni,nj)

	tiles=MeshArray(γ,Int);
	[tiles[τ[i].face][τ[i].i,τ[i].j].=i for i in 1:length(τ)]

	#

	XCtiles=Tiles(τ,exchange(Γ.XC).MA)
	YCtiles=Tiles(τ,exchange(Γ.YC).MA)

	iiTile=tiles[f[1]][i[1],j[1]]; iiFace=τ[iiTile].face

	ii0=minimum(τ[iiTile].i)+Int(ni/2); jj0=minimum(τ[iiTile].j)+Int(nj/2)
	XC0=Γ.XG.f[iiFace][ii0,jj0]; YC0=Γ.YG.f[iiFace][ii0,jj0]
	#XC0=66.5000; YC0=-64.4201

	#
	
	fig1 = Figure(size = (900,600), backgroundcolor = :grey95)
	ax1 = Axis(fig1[1,1],xlabel="longitude",ylabel="latitude")

	scatter!(ax1,XCtiles[iiTile][:],YCtiles[iiTile][:],marker=:+,color=:blue)
	scatter!(ax1,[XC0],[YC0],color=:red)
	scatter!(ax1,lon[:],lat[:],color=:green)

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

	fig2 = Figure(size = (900,600), backgroundcolor = :grey95)
	ax2 = Axis(fig2[1,1],xlabel="x",ylabel="y")

	scatter!(ax2,x_grid[:],y_grid[:],marker=:+,color=:blue)
	scatter!(ax2,[0.],[0.],color=:red)
	scatter!(ax2,x_quad[ii,:][:],y_quad[ii,:][:],color=:orange)
	scatter!(ax2,x_trgt,y_trgt,color=:green)

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
	
	fig3 = Figure(size = (900,600), backgroundcolor = :grey95)
	ax3 = Axis(fig3[1,1],xlabel="longitude",ylabel="latitude")

	scatter!(ax3,XCtiles[iiTile][:],YCtiles[iiTile][:],marker=:+,color=:blue)
	scatter!(ax3,[XC0],[YC0],color=:red,marker=:diamond,markersize=24.0)
	scatter!(ax3,lon,lat,color=:orange,markersize=24.0)
	scatter!(ax3,lon_a,lat_a,color=:black,marker=:star4,markersize=24.0)
	#scatter!(ax3,lon_b,lat_b,color=:black,marker=:star4,markersize=24.0)
	#scatter!(ax3,lon_c,lat_c,color=:black,marker=:star4,markersize=24.0)
	
	(fig1,fig2,fig3)
end

##

function smoothing_demo(Rini_a,Rend_a)	
	fig = Figure(size = (600,600), backgroundcolor = :grey95)

	ax1 = Axis(fig[1,1])
	hm1=heatmap_globalmap!(ax1,Rini_a,colorrange=(-0.25,0.25))
	ax2 = Axis(fig[1,2])
	hm2=heatmap_globalmap!(ax2,Rend_a,colorrange=(-0.25,0.25))
	Colorbar(fig[1,3], hm2, height = Relative(0.65))
	
	fig
end


## convert polygon data to use with Makie

"""
	pol_to_Makie(tmp2::Vector)

Convert output of `read_json` or `read_shp`` to a vector of `LineString`.
"""
function pol_to_Makie(pol::Vector)
	tmp2=[GI.coordinates(a.geometry) for a in pol]
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

##

"""
    basemap(lon,lat,basemap)

```
using MeshArrays, CairoMakie, DataDeps, JLD2
lon,lat,earth_img=demo.get_basemap()
plot_examples(:basemap,lon,lat,earth_img)
```
"""
function basemap(lon,lat,basemap)
    #fig = Figure(size = (1200, 800)) #, backgroundcolor = :grey80)
	fig=with_theme(Figure,theme_light())
    ax = Axis(fig[1, 1])
	#im=image!(ax,lon[:,1],lat[1,:],0.5 .+0.5*Gray.(basemap))
	#im=image!(ax,lon[:,1],lat[1,:],basemap)
	im=image!(ax,lon[1,1]..lon[end,1],lat[1,1]..lat[1,end],basemap)

	hidedecorations!(ax)

	#fig,ax,im
    fig
end

"""
```
import MeshArrays, DataDeps, Shapefile
pol=MeshArrays.Dataset("countries_shp1")

import CairoMakie, Proj
lon0=-160
proj=Proj.Transformation(MA_preset=2,lon0=lon0)
MeshArrays.plot_examples(:baseproj,proj,lon0,pol=pol)
```
"""
function baseproj(proj,lon0; pol=[])
	fi0=Figure(size=(900,600),fontsize=24)
	ax0=Axis(fi0[1,1],backgroundcolor = :gray20)
	pr_ax=ProjAxis(ax0; proj=proj,lon0=lon0)
	lines!(pr_ax,polygons=pol_to_Makie(pol);color=:white, linewidth = 0.5)
	grid_lines!(pr_ax;color=:yellow,linewidth=0.5)
	fi0
end

############################################################
#                                                          #
#     Splitting Of LineString at chosen longitude          #
#                                                          #
############################################################

#This is needed to fix e.g. coast line displays when lon_0 is not 0 but cutting polygons at lon_0+-180.

import Base: split

getlon(p::GeometryBasics.Point) = p[1]

###
function split(tmp::GeometryBasics.LineString, lon0::Real)
	lon0<0.0 ? lon1=lon0+180 : lon1=lon0-180 

	linenodes = GeometryBasics.coordinates(tmp)  # get coordinates of line nodes
	# Find nodes that are on either side of lon0
	cond = getlon.(linenodes) .>= lon1
	# Find interval starts and ends
	end_cond = diff(cond)  # nonzero values denote ends of intervals
	end_inds = findall(!=(0), end_cond)
	start_inds = [firstindex(linenodes);  end_inds .+ 1]  # starts of intervals
	end_inds = [end_inds; lastindex(linenodes)]  # ends of intervals
	# do the splitting
	split_coords = view.(Ref(linenodes), UnitRange.(start_inds, end_inds))  # For each start-end pair, get those coords
	# reconstruct lines from points
	split_lines = if isa(split_coords[1][1],Point)
		[GeometryBasics.LineString(a[:]) for a in split_coords]
	else
		GeometryBasics.LineString.(split_coords)
	end
end

function split(tmp::Vector{<:GeometryBasics.LineString}, lon0::Real)
	[split(a,lon0) for a in tmp]
end

###
split(tmp::GeometryBasics.LineString,dest::Observable) = @lift(split(tmp, $(dest)))

function split(tmp::Vector{<:GeometryBasics.LineString},dest::Observable)
	@lift([split(a,$(dest)) for a in tmp])
end

###
function split(tmp::GeometryBasics.LineString,dest::String)
	if occursin("+lon_0",dest)
		tmp1=split(dest)
		tmp2=findall(occursin.(Ref("+lon_0"),tmp1))[1]
		lon_0=parse(Float64,split(tmp1[tmp2],"=")[2])
		split(tmp,lon_0)
	else
		tmp
	end
end

function split(tmp::Vector{<:GeometryBasics.LineString},dest::String)
	[split(a,dest) for a in tmp]
end

###

split(tmp::Observable,dest::Observable) = @lift(split($(tmp), $(dest)))
split(tmp::Observable,dest::String) = @lift(split($(tmp), (dest)))

##

abstract type AbstractPrAxis <: Makie.AbstractAxis end

struct PrAxis <: AbstractPrAxis
	ax::Axis
	proj::Any
	lon0::Float64
end
 
"""
    ProjAxis(;proj=(x->x),lon0=0.0)

Presets for Proj in MeshArrays.    
"""
function ProjAxis(ax;proj=(x->x),lon0=0.0,omit_grid_lines=true,polygons=Any[])

#    hidespines!(ax)
    hidedecorations!.(ax)
	pr_ax=PrAxis(ax,proj,lon0)
	!omit_grid_lines ? grid_lines!(pr_ax,color=:black,linewidth=0.5) : nothing
	!isempty(polygons) ? lines!(pr_ax; polygons=pol_to_Makie(polygons),color=:black,linewidth=0.5) : nothing
	
	pr_ax
end

function grid_lines!(pr_ax::PrAxis;kwargs...)
	ii=[i for i in -180:30:180, j in -90:1.0:90]';
    jj=[j for i in -180:30:180, j in -90:1.0:90]';
    xl=vcat([[ii[:,i]; NaN] for i in 1:size(ii,2)]...)
    yl=vcat([[jj[:,i]; NaN] for i in 1:size(ii,2)]...)
    tmp=pr_ax.proj.(xl[:],yl[:])
	xl=[a[1] for a in tmp]
	yl=[a[2] for a in tmp]
    lines!(xl,yl; kwargs...)

    tmp=circshift(-179.5:1.0:179.5,-pr_ax.lon0)
    ii=[i for i in tmp, j in -90:30:90];
    jj=[j for i in tmp, j in -90:30:90];
    xl=vcat([[ii[:,i]; NaN] for i in 1:size(ii,2)]...)
    yl=vcat([[jj[:,i]; NaN] for i in 1:size(ii,2)]...)
    tmp=pr_ax.proj.(xl[:],yl[:])
	xl=[a[1] for a in tmp]
	yl=[a[2] for a in tmp]
    lines!(xl,yl; kwargs...)
end

circshift_etc(pr_ax::PrAxis,lon,lat,field) = begin
	if length(size(lon))>1
		lo=lon[:,1]; la=lat[1,:]
	else
		lo=lon[:]; la=lat[:]
	end

	dx=-calc_shift(lo,pr_ax.lon0)

	shlo = circshift(lo,dx)
	csfield = circshift(field,(dx,0))
	
	lon=[i for i in shlo, j in la]
	lat=[j for i in shlo, j in la]
	
	tmp=pr_ax.proj.(lon[:],lat[:])
	x=[a[1] for a in tmp]
	y=[a[2] for a in tmp]
	x=reshape(x,size(lon))
	y=reshape(y,size(lon))

	return x,y,csfield
end	

function surface!(pr_ax::PrAxis,lon,lat,field; color=Float64[], kwargs...)
	x,y,csfield=circshift_etc(pr_ax::PrAxis,lon,lat,field)
	_,_,cscolor=circshift_etc(pr_ax::PrAxis,lon,lat,color)
	surface!(pr_ax.ax,x,y,csfield; color=cscolor, kwargs...)
end

function heatmap!(pr_ax::PrAxis,lon,lat,field; kwargs...)
	x,y,csfield=circshift_etc(pr_ax::PrAxis,lon,lat,field)
	surface!(pr_ax.ax,x,y,0*csfield; color=csfield, kwargs...)
end

function contourf!(pr_ax::PrAxis,lon,lat,field; kwargs...)
	x,y,csfield=circshift_etc(pr_ax::PrAxis,lon,lat,field)
	surface!(pr_ax.ax,x,y,0*csfield; color=csfield, kwargs...)
end

function contour!(pr_ax::PrAxis,lon,lat,field; kwargs...)
	x,y,csfield=circshift_etc(pr_ax::PrAxis,lon,lat,field)
	surface!(pr_ax.ax,x,y,0*csfield; color=csfield, kwargs...)
end

function lines!(pr_ax::PrAxis;polygons=Any[], kwargs...)
	if eltype(polygons)<:Makie.LineString
		po=[Makie.coordinates.(p) for p in split(polygons,pr_ax.lon0)]
		po=[[[Point2(pr_ax.proj(p[1],p[2])) for p in pp] for pp in ppp] for ppp in po]
		[[lines!(pr_ax.ax,pp;kwargs...) for pp in ppp] for ppp in po]
	else
		@warn "untested input type"
	end
end

function scatter!(pr_ax::PrAxis,lon,lat,kargs...; kwargs...)
	tmp=pr_ax.proj.(collect(lon)[:],collect(lat)[:])
	x=[a[1] for a in tmp]
	y=[a[2] for a in tmp]
	x=reshape(x,size(lon))
	y=reshape(y,size(lon))
	scatter!(pr_ax.ax,x,y,kargs...; kwargs...)
end

plot(x::MeshArray; kwargs...) = begin
	if ndims(x) == 1
		heatmap(x; kwargs...)
	else
		heatmap(x[:,1]; kwargs...)
	end
end

##

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

end # module



