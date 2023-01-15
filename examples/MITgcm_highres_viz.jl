
using MeshArrays, GLMakie, Colors

"""
Main workflow.

```
#grid and SSH field
pth="MITgcm_highres_sample/"
γ,Γ=grid_highres_load(pth)

#interpolation coefficients for plotting
dx=0.1; 
lat=[j for i=-179.95:dx:179.95, j=-89.95:dx:89.95]; 
lon=[i for i=-179.95:dx:179.95, j=-89.95:dx:89.95];
(f,i,j,c)=knn(Γ.XC,Γ.YC,vec(lon),vec(lat));

#get data and color range
log_grad_rng=calc_log_grad("THETA")

#plotting
fig=earth_view(log_grad_rng...)

#saving to file
save(joinpath(tempdir(),"tmp.png"),fig)
```
"""
function grid_highres(pth="./")
    grDir=pth
    nFaces=5
    grTopo="LatLonCap"
    ioPrec=Float32
    ioSize=[90 1170]*24
    facesSize=[24 .*(90, 270), 24 .*(90, 270), 24 .*(90, 90), 24 .*(270, 90), 24 .*(270, 90)]
    gcmgrid(grDir,grTopo,nFaces,facesSize, ioSize, ioPrec, read, write)
end

function grid_highres_load(pth="./")
    γ=grid_highres(joinpath(pth,"grid"))
    Γ=GridLoad(γ) #,option="full"
    DXC=γ.read(joinpath(pth,"grid","DXC.data"),MeshArray(γ))
    DYC=γ.read(joinpath(pth,"grid","DYC.data"),MeshArray(γ))
    AngleCS=γ.read(joinpath(pth,"grid","AngleCS.data"),MeshArray(γ))
    AngleSN=γ.read(joinpath(pth,"grid","AngleSN.data"),MeshArray(γ))
    XG=γ.read(joinpath(pth,"grid","XG.data"),MeshArray(γ))
    YG=γ.read(joinpath(pth,"grid","YG.data"),MeshArray(γ))
    Γ=merge(Γ,(DXC=DXC,DYC=DYC,AngleCS=AngleCS,AngleSN=AngleSN,XG=XG,YG=YG))
    γ,Γ
end

## Helper functions

landmsk!(SSH) = SSH[findall(SSH.==0.0)].=NaN

function to_logvel(u,v)
    vel=sqrt.(u.^2 + v.^2)
    γ.read(log10.(γ.write(vel)),vel)
end

function vel_ocn()
    u=γ.read(joinpath(pth,"mit_output","U","U.0000700720.data"),MeshArray(γ))
    landmsk!(u)
    v=γ.read(joinpath(pth,"mit_output","V","V.0000700720.data"),MeshArray(γ))
    landmsk!(v)

    (u, v)=UVtoUEVN(u, v,Γ)
    logvel=to_logvel(u,v)
    logvel_lonlat=reshape(Interpolate_knn(logvel,f,i,j),size(lon))
    logvel_lonlat,(-2.0,0.0)
end

function vel_atm()
    u=γ.read(joinpath(pth,"mit_output","geo5_u10m",
      "geo5_u10m.0000700720.data"),MeshArray(γ))
    v=γ.read(joinpath(pth,"mit_output","geo5_v10m",
      "geo5_v10m.0000700720.data"),MeshArray(γ))

    logvel=to_logvel(u,v)
    logvel_lonlat=reshape(Interpolate_knn(logvel,f,i,j),size(lon))
    logvel_lonlat,(0.0,1.2)
end

function sst_ocn()
    Θ=γ.read(joinpath(pth,"mit_output","THETA","THETA.0000700720.data"),MeshArray(γ))
    landmsk!(Θ)
    Θ_lonlat=reshape(Interpolate_knn(Θ,f,i,j),size(lon))
    Θ_lonlat,(2.0,30.0)
end

function sss_ocn()
    Θ=γ.read(joinpath(pth,"mit_output","SALT","SALT.0000700720.data"),MeshArray(γ))
    landmsk!(Θ)
    Θ_lonlat=reshape(Interpolate_knn(Θ,f,i,j),size(lon))
    Θ_lonlat,(32.0,36.0)
end

function oceQnet()
    Θ=γ.read(joinpath(pth,"mit_output","oceQnet",
      "oceQnet.0000700720.data"),MeshArray(γ))
    landmsk!(Θ)
    Θ_lonlat=reshape(Interpolate_knn(Θ,f,i,j),size(lon))
    Θ_lonlat,(-1.0,1.0)
end

function SSH()
    Θ=γ.read(joinpath(pth,"mit_output","SSH",
      "SSH.0000700720.data"),MeshArray(γ))
    landmsk!(Θ)
    Θ_lonlat=reshape(Interpolate_knn(Θ,f,i,j),size(lon))
    Θ_lonlat,(-1.0,1.0)
end

function calc_log_grad(v="SSH")
    #read variable
    SSH=γ.read(joinpath(pth,"mit_output",v,v*".0000700720.data"),MeshArray(γ))

    #land masking
    SSH[findall(SSH.==0.0)].=NaN;

    #compute gradient magnitude
    (dDdx, dDdy)=gradient(SSH,Γ)
    (dDdx, dDdy)=UVtoUEVN(dDdx, dDdy,Γ)
    dD=sqrt.(dDdx.^2 + dDdy.^2)

    if v=="SSH"
      rng=(-7.0,-5.0)
    elseif v=="THETA"
      rng=(-6.0,-4.0)
    elseif v=="SALT"
      rng=(-7.0,-5.0)
    else
      error("unknown color range")
    end

    #return log10
    γ.read(log10.(γ.write(dD)),dD),rng
end

## interpolate / knn

"""
    Interpolate_knn(A,f,i,j)

```
dx=0.1; 
lat=[j for i=-179.95:dx:179.95, j=-89.95:dx:89.95]; 
lon=[i for i=-179.95:dx:179.95, j=-89.95:dx:89.95];
(f,i,j,c)=knn(Γ.XC,Γ.YC,vec(lon),vec(lat));
logdD_lonlat=reshape(Interpolate_knn(logdD,f,i,j),size(lon));
```
"""
Interpolate_knn(A::MeshArray,f,i,j) = [A[f[jj]][i[jj],j[jj]] for jj in 1:length(i)]

function interp_coeffs(Γ)
    lon=[i for i=-179.75:0.5:179.75, j=-89.75:0.5:89.75]
    lat=[j for i=-179.75:0.5:179.75, j=-89.75:0.5:89.75]
    
    (f,i,j,w)=InterpolationFactors(Γ,vec(lon),vec(lat))

    return (lon=lon, lat=lat, f=f, i=i, j=j, w=w)
end

## plot interpolated result

if false
    ii=1:size(lon,1)
    jj=1:size(lon,2)
    heatmap(lon[ii,1],lat[1,jj],tmp[ii,jj],colorrange=(-7.0,-5.0))
end

##

using GLMakie, FileIO
using Downloads: download

function earth_view(tmp,rng)

    isa(tmp,MeshArray) ? tmp=reshape(Interpolate_knn(tmp,f,i,j),size(lon)) : nothing

    earth_tmp=reverse(permutedims(tmp),dims=1)
    #earth_img = load(download("https://upload.wikimedia.org/wikipedia/commons/5/56/Blue_Marble_Next_Generation_%2B_topography_%2B_bathymetry.jpg"))
    earth_img=load(joinpath(pth,"images",
      "Blue_Marble_Next_Generation_topography_bathymetry.jpg"))
    n = 1024 ÷ 4 # 2048
    θ = LinRange(0, π, n)
    φ = LinRange(0, 2π, 2 * n)
    x = [cos(φ) * sin(θ) for θ in θ, φ in φ]
    y = [sin(φ) * sin(θ) for θ in θ, φ in φ]
    z = [cos(θ) for θ in θ, φ in φ]
    fig = Figure(resolution = (1200, 800), backgroundcolor = :grey80)
    ax = LScene(fig[1, 1], show_axis = false)
    surface!(ax, 0.99*x, 0.99*y, 0.99*z; color = Gray.(earth_img))
    surface!(ax, x, y, z; color = earth_tmp, shading = false,
        lightposition = Vec3f(-2, -3, -3), ambient = Vec3f(0.8, 0.8, 0.8),
        backlight = 1.5f0, colorrange=rng, colormap=:thermal)
    zoom!(ax.scene, cameracontrols(ax.scene), 0.6)
    rotate!(ax.scene, Vec3f(0.4, -0.8, 1), 0.6)
    fig
end

function latlon_view(tmp,rng)
    isa(tmp,MeshArray) ? tmp=reshape(Interpolate_knn(tmp,f,i,j),size(lon)) : nothing

    earth_img=load(joinpath(pth,"images",
      "Blue_Marble_Next_Generation_topography_bathymetry.jpg"))
    earth_img=reverse(permutedims(earth_img),dims=2)

    fig = Figure(resolution = (1200, 800), backgroundcolor = :grey80)
    ax = Axis(fig[1, 1])
    image!(ax,lon[:,1],lat[1,:],Gray.(earth_tmp))
    heatmap!(ax,lon[:,1], lat[1,:], tmp; colorrange=rng, colormap=:thermal)

    #hidedecorations!(current_axis())
    hidedecorations!(ax)

    fig
end


## display one face data

function plot_face_grad(SSH)
    tmp2=SSH[5];
    tmp2[tmp2.==0].=NaN;
    tmp3=sqrt.(diff(tmp2,dims=1)[:,1:end-1].^2+diff(tmp2,dims=2)[1:end-1,:].^2);
    ii=findall((!isnan).(log10.(tmp3)));
    scatter(x[ii],y[ii],color=log10.(tmp3)[ii],colorrange=(-3.0,-1.0),markersize=0.1)
end

function plot_face_scatter(Γ,dD)
    x=Float64.(Γ.XC[5])
    y=Float64.(Γ.YC[5])
    ii=findall((!isnan).(dD[5]))
    ii=shuffle(ii)[1:1000000]
    scatter(x[ii],y[ii],color=log10.(dD[5][ii]),colorrange=(-7.0,-5.0),markersize=2.0)
end
