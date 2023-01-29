module ProjMakie

using CairoMakie, Proj

"""
    projmap(dat,proj=1)

See `ClimateModels.jl/examples/IPCC.jl` for example/use case.
"""
function projmap(dat,proj=1)
	
	dx=dat.meta.shift
	lons = circshift(dat.lon[:,1],dx)
	lats = dat.lat[1,:]
	field = circshift(dat.var,(dx,0))

	source="+proj=longlat +datum=WGS84"
	if proj==2 
		dest="+proj=eqearth +lon_0=200.0 +lat_1=0.0 +x_0=0.0 +y_0=0.0 +ellps=GRS80"
		lon0=-160.0 #or 200.0 ?
	elseif proj==1
		dest="+proj=wintri"
		lon0=0.0
	elseif proj==3
		dest="+proj=longlat +datum=WGS84 +lon_0=-160.0"
		lon0=-160.0
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
    ax = f[1, 1] = Axis(f, aspect = DataAspect(), title = dat.meta.ttl)
    
    surf = surface!(ax,x,y,0*x; color=field, 
	colorrange=dat.meta.colorrange, colormap=dat.meta.cmap,
        shading = false)

	ii=[i for i in -180:45:180, j in -78.5:1.0:78.5]';
    jj=[j for i in -180:45:180, j in -78.5:1.0:78.5]';
    xl=vcat([[ii[:,i]; NaN] for i in 1:size(ii,2)]...)
    yl=vcat([[jj[:,i]; NaN] for i in 1:size(ii,2)]...)
    tmp=trans.(xl[:],yl[:])
	xl=[a[1] for a in tmp]
	yl=[a[2] for a in tmp]
    proj<3 ? lines!(xl,yl, color = :black, linewidth = 0.5) : nothing

	if proj==2 
	    tmp=circshift(-179.5:1.0:179.5,(-200))
	elseif proj==1
	    tmp=(-179.5:1.0:179.5)
	elseif proj==3
	    tmp=circshift(-179.5:1.0:179.5,(-200))
	end
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

#	all_lines=demo.LineSplit(GeoMakie.coastlines(),lon0)
#	[lines!(ax, l,color=:black,linewidth=1.0) for l in all_lines]
	
	#add colorbar
	Colorbar(f[1,2], surf, height = Relative(0.5))

	f
end

end #module GeoProjMakie
