
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

## AbstractPrAxis, PrAxis, and ProjAxis

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
	#contourf!(pr_ax.ax,x,y,csfield, kwargs...)
end

function contour!(pr_ax::PrAxis,lon,lat,field; kwargs...)
	x,y,csfield=circshift_etc(pr_ax::PrAxis,lon,lat,field)
	surface!(pr_ax.ax,x,y,0*csfield; color=csfield, kwargs...)
	#contour!(pr_ax.ax,x,y,csfield, kwargs...)
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
