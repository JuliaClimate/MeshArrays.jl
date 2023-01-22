
using Revise
using Tyler, GLMakie, JLD2

#file="MITgcm_highres_sample/mit_output_Interpolations/itp_sss_ocn.jld2"
#file=joinpath(tempdir(),"itp_sst_ocn.jld2")

pth0="MITgcm_highres_sample/mit_output_Interpolations/"
list0=CSV.read("MITgcm_highres_interp_list.csv",DataFrame)
file=joinpath(pth0,list0.filename[ii])

f = jldopen(file, "r");

#Projector(lon,lat) = f["itp"](lon,-lat)

itp=Tyler.ForInterpolations.Interpolator(
     f["itp"], 
     Dict(:min_zoom => 1,:max_zoom => 19)
    ); #pb with show method

m = Tyler.Map( Rect2f(-30.0, 38.0, 5.0, 5.0), provider=itp)

