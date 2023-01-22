
pth="MITgcm_highres_sample/"
γ,Γ=grid_highres_load(pth)

knn_file=joinpath(tempdir(),"knn_tmp.jld2")
knn_api =JLD2.jldopen(knn_file, "r");

(i,j,f,c,lon,lat)=(knn_api["i"],knn_api["j"],knn_api["f"],knn_api["c"],
            knn_api["lon"],knn_api["lat"])

include("MITgcm_highres_viz.jl")
import Interpolations, JLD2

"""
```
to_itp(log_grad_rng,"itp_log_grad_rng.jld2")
```
"""
function to_itp(in_func,out_file)

    filename=joinpath(tempdir(),out_file)
    (z,(z0,z1))=in_func();

    z_latlon=z
    z_latlon=[min(i,z1) for i in z_latlon]
    z_latlon=[max(i,z0) for i in z_latlon]
    z_latlon=(z_latlon .-z0)./(z1-z0)
    z_latlon[isnan.(z_latlon)].=0.5

    dx=knn_api["dx"]
    nodes=(Float32.(-180.0:dx:180.0), Float32.(-90.0:dx:90.0))
    itp = Interpolations.interpolate(nodes, Float32.(reverse(z_latlon,dims=2)),
            Interpolations.Gridded(Interpolations.Linear()))
           
    save(filename,"itp",itp)
end


