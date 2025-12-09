
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
	elseif ID==:meridional_overturning
		meridional_overturning(stuff...)
	elseif ID==:gradient_EN
		gradient_EN(stuff...)
	elseif ID==:gradient_xy
		gradient_xy(stuff...)
	elseif ID==:basemap
		basemap(stuff...)
	elseif ID==:baseproj
		baseproj(stuff...;kwargs...)
	elseif ID==:polygons_plot_dev1
		polygons_plot_dev1(stuff...;kwargs...)
	elseif ID==:polygons_plot
		polygons_plot(stuff...;kwargs...)
	else
		println("unknown plot ID")
	end
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

function meridional_overturning(Γ,ov)
	x=vec(-89.0:89.0); y=reverse(vec(Γ.RF[1:end-1])); #coordinate variables
	z=reverse(ov,dims=2); #z[z.==0.0].=NaN

	fig1 = Figure(size = (900,400),markersize=0.1)
	ax1 = Axis(fig1[1,1], title="Meridional Overturning Streamfunction (in Sv)")
	hm1=contourf!(ax1,x,y,1e-6*z,levels=(-40.0:5.0:40.0))
	contour!(ax1,x,y,1e-6*z,color=:black)
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

##

function polygons_plot_dev1_basis(sphere_view)
	f=Figure(size=(600,400))
	if sphere_view
		scene = LScene(f[1, 1])
		cam = Makie.Camera3D(scene.scene, projectiontype = Makie.Perspective)
		sphere = Sphere(Point3f(0), 0.99f0)
		mesh!(scene,sphere, color = :black)
		zoom!(scene.scene,cam,2.0)
		a=scene
	else
		a=Axis(f[1,1])
	end	
	f,a
end

cols=[:blue :green :orange :black :red :violet]

"""
    polygons_plot_dev1(pols,pols3D, facets=1:6; 
		colors=cols, sphere_view=true)

```
polygons_plot_dev1(pols,pols3D,1:2)
```
"""
function polygons_plot_dev1(pols,pols3D; facets=1:6, colors=cols, sphere_view=true)
	f,a=polygons_plot_dev1_basis(sphere_view)
	if sphere_view
		xyz=vcat([p[:] for p in pols3D[facets]]...)
		col=vcat([fill(c,32*32) for c in facets]...)
#		lin=vcat([fill(0.5,32*32) for c in facets]...)
		plot!(a,xyz,color=col)
	else
		xyz=vcat([p[:] for p in pols[facets]]...)
#		col=vcat([fill(c,32*32) for c in facets]...)
#		lin=vcat([fill(0.5,32*32) for c in facets]...)
		lines!(a,xyz)
	end
	f
end
