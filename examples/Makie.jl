
#possible issues
#1. intermediate yy needed; MeshArrays errors otherwise
#2. lack of exch_z, which would allow using corner points 
#4. Float32 / Float64 dispatch may be an issue in e.g.
#   `pi/2 .-Γ.YC*pi/180`

#    XC=γ.read(γ.write(Γ.XC),MeshArray(γ,Float64))
#    YC=γ.read(γ.write(Γ.YC),MeshArray(γ,Float64))

using GLMakie

function plot_as_sphere(Γ)
    yy=pi/2 .-Γ.YC*pi/180
    x=sin.(yy)*cos.(Γ.XC*pi/180)
    y=sin.(yy)*sin.(Γ.XC*pi/180)
    z=cos.(yy)
    d=Γ.Depth
    crng=(minimum(d),maximum(d))
    minimum(d)==maximum(d) ? crng=(0.9*crng[1],1.1*crng[1]) : nothing

    c=[:gold,:magenta,:DeepSkyBlue,:cyan,:red,:black]
    az=Node(0.4π)
    el=Node(0.2π)
    fig = Figure(resolution = (900,900), backgroundcolor = :grey90)
    ax = Axis3(fig, aspect = :data, viewmode = :fitzoom, #perspectiveness = 0.5,
        azimuth = az, elevation = el,)
    for i=1:length(z.fSize)
        surface!(ax, x[i], y[i], z[i], color = d[i], colorrange = crng, colormap = :grays)
        wireframe!(x[i],y[i],z[i], overdraw = false, linewidth = 1.0,color=c[i])
    end
    #hidedecorations!(ax)
    #hidespines!(ax)
    fig[1,1] = ax
    fig
end

function plot_as_plane(Γ)
    x=Float64.(Γ.XC)
    y=Float64.(Γ.YC)
	z=0*y
    d=Float64.(Γ.Depth)
    crng=(minimum(d),maximum(d))
    minimum(d)==maximum(d) ? crng=(0.9*crng[1],1.1*crng[1]) : nothing

    c=[:gold,:magenta,:blue,:cyan,:red,:black]
    fig = Figure(resolution = (900,900), backgroundcolor = :grey90)
    ax = Axis3(fig[1,1])
    for i=1:length(x.fSize)
        surface!(ax, x[i], y[i], z[i], color = d[i], colorrange = crng, colormap = :grays)
		#contourf!(x[i], y[i], d[i], colorrange = (0.0, 5e3), colormap = :grays)
		wireframe!(x[i],y[i], z[i], overdraw = false, linewidth = 1.0,color=c[i])
    end
    #hidedecorations!(ax)
    #hidespines!(ax)
    fig[1,1] = ax
    fig
end