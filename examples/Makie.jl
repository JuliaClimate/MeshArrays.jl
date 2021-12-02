
#possible issues
#1. intermediate yy needed; MeshArrays errors otherwise
#2. lack of exch_z, which would allow using corner points 
#4. Float32 / Float64 dispatch may be an issue in e.g.
#   `pi/2 .-Γ.YC*pi/180`

#    XC=γ.read(γ.write(Γ.XC),MeshArray(γ,Float64))
#    YC=γ.read(γ.write(Γ.YC),MeshArray(γ,Float64))

using MeshArrays
import Meshes, MeshViz
import GLMakie as Mkie
#import CairoMakie as Mkie

col=[:blue,:greenyellow,:magenta,:yellow2,:tomato,:black]

"""
    plot_as_sphere(Γ)

Plot mesh (Γ.XC,Γ.YC) and depth (Γ.Depth) on the surface of the sphere in 3D.
"""
function plot_as_sphere(Γ; title="", az=-0.2π, el=0.2π)
    yy=pi/2 .-Γ.YC*pi/180
    x=sin.(yy)*cos.(Γ.XC*pi/180)
    y=sin.(yy)*sin.(Γ.XC*pi/180)
    z=cos.(yy)
    d=Γ.Depth
    crng=(minimum(d),maximum(d))
    minimum(d)==maximum(d) ? crng=(0.8*crng[1],1.1*crng[1]) : nothing

    fig = Mkie.Figure(resolution = (900,900), backgroundcolor = :grey95)
    ax = Mkie.Axis3(fig, aspect = :data, viewmode = :fitzoom, #perspectiveness = 0.5,
        azimuth = az, elevation = el)
    for i=1:length(z.fSize)
        Mkie.surface!(ax, x[i], y[i], z[i], color = d[i], colorrange = crng, colormap = :grays)
        Mkie.wireframe!(x[i],y[i],z[i], overdraw = false, linewidth = 0.25,color=col[i])
    end
    
    Mkie.hidedecorations!(ax)
    Mkie.hidespines!(ax)
    ax.title=title
    fig[1,1] = ax
    fig
end

"""
    plot_as_plane(Γ)

Plot mesh (Γ.XC,Γ.YC) and depth (Γ.Depth) on a 2D plane.
"""
function plot_as_plane(Γ; title="")
    x=Float64.(Γ.XC)
    y=Float64.(Γ.YC)
	z=0*y
    d=Float64.(Γ.Depth)
    crng=(minimum(d),maximum(d))
    minimum(d)==maximum(d) ? crng=(0.8*crng[1],1.1*crng[1]) : nothing

    fig = Mkie.Figure(resolution = (900,900), backgroundcolor = :grey95)

    ax = Mkie.Axis(fig[1,1])
    for i=1:length(x.fSize)
        Mkie.surface!(ax, x[i], y[i], z[i], color = d[i], colorrange = crng, colormap = :grays)
		#contourf!(x[i], y[i], d[i], colorrange = (0.0, 5e3), colormap = :grays)
		Mkie.wireframe!(x[i],y[i], z[i], overdraw = false, linewidth = 0.25,color=col[i])
    end
    #hidedecorations!(ax)
    #hidespines!(ax)
    ax.title=title
    fig[1,1] = ax
    fig
end


"""
    plot_as_scatter(Γ; title="")

Plot mesh (Γ.XC,Γ.YC) as scatter of points.
"""
function plot_as_scatter(Γ; title="")
    x=Float64.(Γ.XC)
    y=Float64.(Γ.YC)

    fig = Figure(resolution = (900,900), backgroundcolor = :grey95)

    ax = Axis(fig[1,1])
    for i=1:length(x.fSize)
        lines!(ax, x[i][:], y[i][:], color = col[i])
        scatter!(ax, x[i][:], y[i][:], color = col[i], marker = 'o')
    end
    #hidedecorations!(ax)
    #hidespines!(ax)
    ax.title=title
    fig[1,1] = ax
    fig
end

#Γ=GridLoad(GridSpec("LatLonCap",MeshArrays.GRID_LLC90))
function plot_as_plane_b()
    Mkie.set_theme!(Mkie.theme_light())
    fig = Mkie.Figure(resolution = (900,600),markersize=0.1)
    ax = Mkie.Axis(fig[1,1],xlabel="longitude",ylabel="latitude")

    col=[:red,:green,:blue,:magenta,:cyan,:black]

    for f in 1:length(Γ.XC)
        tmp=[Meshes.Point2(Γ.XC[f][i],Γ.YC[f][i]) for i in eachindex(Γ.XC[f])]
        MeshViz.viz!(ax,tmp,elementcolor=col[f],markersize=2.0)
    end

    fig,ax
end

#Γ=GridLoad(GridSpec("LatLonCap",MeshArrays.GRID_LLC90))
function plot_as_sphere_b()
    Mkie.set_theme!(Mkie.theme_light())
    fig = Mkie.Figure(resolution = (900,600),markersize=0.1)
    ax = Mkie.Axis3(fig[1,1])

    col=[:red,:green,:blue,:magenta,:cyan,:black]

    yy=pi/2 .-Γ.YC*pi/180
    x=sin.(yy)*cos.(Γ.XC*pi/180)
    y=sin.(yy)*sin.(Γ.XC*pi/180)
    z=cos.(yy)
    #d=Γ.Depth

    #sphr=sample(Meshes.Sphere((0.0,0.0,0.0), 0.99), Meshes.RegularSampling(360, 180))
    #sphr=Meshes.Sphere((0.0,0.0,0.0), 0.99)
    #MeshViz.viz!(ax,sphr)

    sphere = Mkie.Sphere(Mkie.Point3f(0.0), 0.99)
    Mkie.mesh!(ax,sphere, color=:white)

    for f in 1:length(Γ.XC)
        tmp=[Meshes.Point3(x[f][i],y[f][i],z[f][i]) for i in eachindex(Γ.XC[f])]
        MeshViz.viz!(ax,tmp,color=col[f],pointsize=1000.0)
    end

    Mkie.hidedecorations!(ax)
    Mkie.hidespines!(ax)
    Mkie.ylims!(-0.7,0.7)
    Mkie.xlims!(-0.7,0.7)
    Mkie.zlims!(-0.7,0.7)
    
    fig,ax
end
