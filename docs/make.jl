using Documenter, MeshArrays
import PlutoSliderServer, DataDeps, CairoMakie 
import JLD2, GeometryOps, Shapefile, GeoJSON
#import Climatology, MITgcm

ENV["DATADEPS_ALWAYS_ACCEPT"]=true
#MITgcm.getdata("mitgcmsmall")
#Climatology.get_ecco_velocity_if_needed()
MeshArrays.Dataset("GRID_LL360")
MeshArrays.Dataset("GRID_LLC90")
MeshArrays.Dataset("GRID_CS32")
MeshArrays.Dataset("countries_shp1")

MeshArrays.interpolation_setup()

makedocs(;
    sitename = "MeshArrays",
    format   = Documenter.HTML(),
    modules  = [MeshArrays, Base.get_extension(MeshArrays, :MeshArraysDataDepsExt)],
    warnonly = [:cross_references,:missing_docs],
    pages = [
    "Home" => "index.md",
    "Get Started" => "start.md",
    "Main Features" => "main.md",
    "Notebook Tutorials" => "tutorials.md",
    "API documentation" => "API.md",
    "Miscellaneous" => ["detail.md", "videos.md", "dev.md"],
    ],
    authors="gaelforget <gforget@mit.edu>",
)

lst=("basics.jl","geography.jl","vectors.jl",
	"dev/MeshArrays_to_Polygons.jl","dev/GeometryOps_exploration.jl")
pth=("tutorials","tutorials","tutorials",
	"dev","dev")

for ii in 1:length(lst)
    i=lst[ii]
    fil_in=joinpath(@__DIR__,"..", "examples",i)
    fil_out=joinpath(@__DIR__,"build", pth[ii] ,basename(i)[1:end-2]*"html")
    PlutoSliderServer.export_notebook(fil_in)
    mv(fil_in[1:end-2]*"html",fil_out)
    cp(fil_in,fil_out[1:end-4]*"jl")
end

deploydocs(
    repo = "github.com/JuliaClimate/MeshArrays.jl",
)
