
include("MITgcm_highres_viz.jl")

o(n=1) = println("o"^n)

pth="MITgcm_highres_sample/"
γ,Γ=grid_highres_load(pth)

o(1)

## read velocity field

function logvel(pth)
    u=γ.read(joinpath(pth,"mit_output","U","U.0000700720.data"),MeshArray(γ))
    landmsk!(u)
    v=γ.read(joinpath(pth,"mit_output","V","V.0000700720.data"),MeshArray(γ))
    landmsk!(v)

    (u, v)=UVtoUEVN(u, v,Γ)
    vel=sqrt.(u.^2 + v.^2)
    γ.read(log10.(γ.write(vel)),vel)
end

o(2)

## interpolation coefficients (to lat-lon grid)

dx=0.02
lat=[j for i=-180.0:dx:180.0, j=-90.0:dx:90.0] 
lon=[i for i=-180.0:dx:180.0, j=-90.0:dx:90.0]
(f,i,j,c)=knn(Γ.XC,Γ.YC,vec(lon),vec(lat))

save(joinpath(tempdir(),"knn_tmp.jld2"),f,i,j,c)

o(3)

## rescale to 0-1 range mapped to z0-z1

z=logvel(pth)
(z0,z1)=(-2.0,0.0)

z_latlon=reshape(Interpolate_knn(z,f,i,j),size(lon))
z_latlon=[min(i,z1) for i in z_latlon]
z_latlon=[max(i,z0) for i in z_latlon]
z_latlon=(z_latlon .-z0)./(z1-z0)
z_latlon[isnan.(z_latlon)].=0.5

o(4)

## prepare for interpolation

import Interpolations
nodes=(-180.0:dx:180.0, -90.0:dx:90.0)
itp = Interpolations.interpolate(nodes, z_latlon, 
    Interpolations.Gridded(Interpolations.Linear()))

save(joinpath(tempdir(),"itp_tmp.jld2"),"itp",itp)

o(5)

## test

if false
    using Tyler, MapTiles
    tile=Tile(5,5,5)
    (lon,lat)=Tyler.tile2positions(tile)
    Depth_tile=itp.(lon,lat)
end

