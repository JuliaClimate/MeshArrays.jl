
module MeshArraysMakieExt

using MeshArrays, Makie

import MeshArrays: plot_ocean_basins, plot_one_section, land_mask, plot_cell_area
import MeshArrays: projmap, interpolation_demo, read_polygons

import Makie: heatmap
LineString=Makie.LineString
Observable=Makie.Observable

"""
    MeshArrays.plot(c::MeshArray)

This is a test.

```
using MeshArrays, GLMakie

Γ=GridLoad(GridSpec("LatLonCap",MeshArrays.GRID_LLC90))

plot(Γ.YC,facet=3)
"""
function Makie.plot(c::MeshArray;facet=1)
    heatmap(c[facet])
end

##

function heatmap(MS::MeshArray,λ;
    title="",colorrange=[],colormap=:spring,colorbar=true)
    DD=Interpolate(MS,λ.f,λ.i,λ.j,λ.w)
	DD=reshape(DD,size(λ.lon))

    fig = Figure(resolution = (900,600), backgroundcolor = :grey95)
    ax = Axis(fig[1,1],xlabel="longitude",ylabel="latitude",title=title)
    if !isempty(colorrange) 
        hm1=heatmap!(ax,λ.lon[:,1],λ.lat[1,:],DD,colormap=colormap,colorrange=cr)
    else
        hm1=heatmap!(ax,λ.lon[:,1],λ.lat[1,:],DD,colormap=colormap)
    end
    colorbar ? Colorbar(fig[1,2], hm1, height = Relative(0.65)) : nothing
    fig
end

function plot_one_section(fig,lons,lats,my_section)	
    ax=current_axis(fig)
	scatter!(ax,my_section.lon[:],my_section.lat[:],color=:blue,markersize=2.0)
	scatter!(ax,my_section.lon[:],my_section.lat[:],color=:black,markersize=4.0)
	scatter!(ax,lons[:],lats[:],color=:blue)
	fig
end

##

function plot_ocean_basins(Γ,λ,basins,basin_nam)
#    fig=heatmap(λ.μ*Γ.Depth,λ,colormap=:grays,colorbar=false,colorrange=(0,10000))
#    ax=current_axis(fig)
    fig = Figure(resolution = (900,600), backgroundcolor = :grey95,colormap=:thermal)
    ax = Axis(fig[1,1],xlabel="longitude",ylabel="latitude",title="ocean basin IDs")

	basinID=findall(basins.name.==basin_nam)[1]
	
	mx=maximum(basins.map)
	for ff in 1:length(Γ.RAC)
		col=λ.μ[ff][:].*(basins.map[ff][:].==basinID)
		kk=findall(col.>0.0)
		!isempty(kk) ? scatter!(ax,Γ.XC[ff][kk],Γ.YC[ff][kk],color=:red,markersize=2.0) : nothing
		kk=findall((col.==0.0).*(!isnan).(λ.μ[ff][:]))
		!isempty(kk) ? scatter!(ax,Γ.XC[ff][kk],Γ.YC[ff][kk],color=basins.map[ff][kk],
			colorrange=(0.0,mx),markersize=2.0,colormap=:lisbon) : nothing
	end
	Colorbar(fig[1,2], colormap=:lisbon, colorrange=(0.0, mx), height = Relative(0.65))

	fig
end

function land_mask(Γ)
    μ =Γ.hFacC[:,1]
    μ[findall(μ.>0.0)].=1.0
    μ[findall(μ.==0.0)].=NaN
    μ
end


function plot_cell_area(Γ,λ,faceID)
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

##


"""
    projmap(data,trans; omit_lines=false)

Inputs:

- data is a `NamedTuple`
- trans is a `Proj.Transformation`

Use examples:

- `MeshArrays.jl/examples/geography.jl`
- `ClimateModels.jl/examples/IPCC.jl`
"""
function projmap(data,trans; omit_lines=false)
	
    lon0=data.meta.lon0
	dx=-calc_shift(data.lon[:,1],data.meta.lon0)

    lons = circshift(data.lon[:,1],dx)
	lats = data.lat[1,:]
	field = circshift(data.var,(dx,0))

	lon=[i for i in lons, j in lats]
    lat=[j for i in lons, j in lats]

    tmp=trans.(lon[:],lat[:])
	x=[a[1] for a in tmp]
	y=[a[2] for a in tmp]
    x=reshape(x,size(lon))
    y=reshape(y,size(lon))

	f = Figure()
    ax = f[1, 1] = Axis(f, aspect = DataAspect(), title = data.meta.ttl)
    
    surf = surface!(ax,x,y,0*x; color=field, 
	colorrange=data.meta.colorrange, colormap=data.meta.cmap,
        shading = false)

	ii=[i for i in -180:45:180, j in -78.5:1.0:78.5]';
    jj=[j for i in -180:45:180, j in -78.5:1.0:78.5]';
    xl=vcat([[ii[:,i]; NaN] for i in 1:size(ii,2)]...)
    yl=vcat([[jj[:,i]; NaN] for i in 1:size(ii,2)]...)
    tmp=trans.(xl[:],yl[:])
	xl=[a[1] for a in tmp]
	yl=[a[2] for a in tmp]
    !omit_lines ? lines!(xl,yl, color = :black, linewidth = 0.5) : nothing

    tmp=circshift(-179.5:1.0:179.5,-lon0)
    ii=[i for i in tmp, j in -75:15:75];
    jj=[j for i in tmp, j in -75:15:75];
    xl=vcat([[ii[:,i]; NaN] for i in 1:size(ii,2)]...)
    yl=vcat([[jj[:,i]; NaN] for i in 1:size(ii,2)]...)
    tmp=trans.(xl[:],yl[:])
	xl=[a[1] for a in tmp]
	yl=[a[2] for a in tmp]
    !omit_lines ? lines!(xl,yl, color = :black, linewidth = 0.5) : nothing

    hidespines!(ax)
    hidedecorations!.(ax)

    po=data.polygons #LineSplitting.LineSplit(data.polygons,lon0)
    po=[[Point2(trans(p[1],p[2])) for p in k] for k in po]
    [lines!(ax,l,color=:black,linewidth=0.5) for l in po]

	#add colorbar
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

##

"""
	read_polygons(file="countries.geojson")

Call `json_to_Makie` or `shp_to_Makie` depending on file extension.

"""
function read_polygons(file::String)
	if !isfile(file)
		error("file not found ($file)")
	elseif occursin(".geojson",file)&&file[end-7:end]==".geojson"
		json_to_Makie(file)
	elseif occursin(".shp",file)&&file[end-3:end]==".shp"
		shp_to_Makie(file)
	else
		error("unknown file extension ($file)")
	end
end

## read data from file

"""
	json_to_Makie(file="countries.geojson")

Call `GeoJSON.read` then `geo2basic`. Return a vector of `LineString`.
"""
function json_to_Makie(file="countries.geojson")
	tmp2=MeshArrays.read_json(file)

	tmp22=Vector{Point2{Float64}}[]
	for l1 in tmp2
		if isa(l1[1][1],Vector{Float64})
			push!(tmp22,geo2basic(l1[1]))
		else
			for l2 in l1
				push!(tmp22,geo2basic(l2[1]))
			end
		end
	end
	
	LineString.(tmp22)
end

"""
	shp_to_Makie(file="countries.geojson")

Done via `Shapefile.shapes` then `geo2basic`. Return a vector of `LineString`.
"""
function shp_to_Makie(file="ne_110m_admin_0_countries.shp")
	tmp2=MeshArrays.read_shp(file)

	tmp22=Vector{Point2{Float64}}[]
	[[[push!(tmp22,geo2basic(l3)) for l3 in l2] for l2 in l1] for l1 in tmp2]
	
	LineString.(tmp22)
end

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


############################################################
#                                                          #
#     Splitting Of LineString at chosen longitude          #
#                                                          #
############################################################

#This is needed to fix e.g. coast line displays when lon_0 is not 0 but cutting polygons at lon_0+-180.

import Base: split

function regroup(tmp::Vector)
	coastlines_custom=LineString[]
	for ii in 1:length(tmp)
		push!(coastlines_custom,tmp[ii][:]...)
	end
	coastlines_custom
end

"""
	split(tmp::Vector{<:LineString},lon0::Real)

Split each `LineString` at `lon0`.
"""
function split(tmp::Vector{<:LineString},lon0::Real)
	regroup([split(a,lon0) for a in tmp])
end

"""
	split(tmp::LineString,lon0::Real)

Split `LineString` at `lon0`.
"""
function split(tmp::LineString,lon0::Real)
	lon0<0.0 ? lon1=lon0+180 : lon1=lon0-180 
	np=length(tmp)
	tmp2=fill(0,np)
	for p in 1:np
		tmp1=tmp[p]
		tmp2[p]=maximum( [(tmp1[1][1]<=lon1)+2*(tmp1[2][1]>=lon1) , (tmp1[2][1]<=lon1)+2*(tmp1[1][1]>=lon1)] )
	end
	if sum(tmp2.==3)==0
		[tmp]
	else
		#println("splitting here")
		jj=[0;findall(tmp2.==3)...;np+1]
		[LineString(tmp[jj[ii]+1:jj[ii+1]-1]) for ii in 1:length(jj)-1]
	end
end

split(tmp::Vector{<:LineString},dest::Observable) = tmp

function split(tmp::Vector{<:LineString},dest::String)
	if occursin("+lon_0",dest)
		tmp1=split(dest)
		tmp2=findall(occursin.(Ref("+lon_0"),tmp1))[1]
		lon_0=parse(Float64,split(tmp1[tmp2],"=")[2])
		regroup(split(tmp,lon_0))
	else
		tmp
	end
end

end # module

