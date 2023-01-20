if false
    include("MITgcm_highres_viz.jl")

    pth="MITgcm_highres_sample/"
    γ,Γ=grid_highres_load(pth)

    dx=0.02
    lat=[j for i=-180.0:dx:180.0, j=-90.0:dx:90.0] 
    lon=[i for i=-180.0:dx:180.0, j=-90.0:dx:90.0]

    i=load("knn_coeffs.jld2","i")
    j=load("knn_coeffs.jld2","j");
    f=load("knn_coeffs.jld2","f")
    c=load("knn_coeffs.jld2","c");
end

if true
    (z,(z0,z1))=SSH();
    filename="itp_SSH.jld2";
else
    (z,(z0,z1))=sst_ocn();
    filename="itp_sst_ocn.jld2";
end    

z_latlon=z
z_latlon=[min(i,z1) for i in z_latlon]
z_latlon=[max(i,z0) for i in z_latlon]
z_latlon=(z_latlon .-z0)./(z1-z0)
z_latlon[isnan.(z_latlon)].=0.5

import Interpolations

nodes=(-180.0:dx:180.0, -90.0:dx:90.0)

itp = Interpolations.interpolate(nodes, z_latlon,
    Interpolations.Gridded(Interpolations.Linear()))
           
save(filename,"itp",itp)


