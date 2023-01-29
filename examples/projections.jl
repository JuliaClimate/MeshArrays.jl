module ProjMakie

using CairoMakie, Proj

"""
    projmap(dat,proj=1)

See `ClimateModels.jl/examples/IPCC.jl` for example/use case.
"""
function projmap(data,proj=1)
	
    lon0=data.meta.lon0
	dx=-calc_shift(data.lon[:,1],data.meta.lon0)

    lons = circshift(data.lon[:,1],dx)
	lats = data.lat[1,:]
	field = circshift(data.var,(dx,0))

	source="+proj=longlat +datum=WGS84"
	if proj==2 
		dest="+proj=eqearth +lon_0=$(lon0) +lat_1=0.0 +x_0=0.0 +y_0=0.0 +ellps=GRS80"
	elseif proj==1
		dest="+proj=wintri +lon_0=$(lon0) "
	elseif proj==3
		dest="+proj=longlat +datum=WGS84 +lon_0=$(lon0)"
	end
	trans = Proj.Transformation(source,dest, always_xy=true) 

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
    proj<3 ? lines!(xl,yl, color = :black, linewidth = 0.5) : nothing

    tmp=circshift(-179.5:1.0:179.5,-lon0)
    ii=[i for i in tmp, j in -75:15:75];
    jj=[j for i in tmp, j in -75:15:75];
    xl=vcat([[ii[:,i]; NaN] for i in 1:size(ii,2)]...)
    yl=vcat([[jj[:,i]; NaN] for i in 1:size(ii,2)]...)
    tmp=trans.(xl[:],yl[:])
	xl=[a[1] for a in tmp]
	yl=[a[2] for a in tmp]
    proj<3 ? lines!(xl,yl, color = :black, linewidth = 0.5) : nothing

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

end #module GeoProjMakie
