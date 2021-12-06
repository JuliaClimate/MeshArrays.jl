using Documenter, PlutoSliderServer, MeshArrays

import OceanStateEstimation
OceanStateEstimation.get_ecco_velocity_if_needed()
MeshArrays.GRID_LL360_download()
MeshArrays.GRID_LLC90_download()
MeshArrays.GRID_CS32_download()

makedocs(
    sitename = "MeshArrays",
    format   = Documenter.HTML(),
    modules  = [MeshArrays],
    pages = [
    "Home" => "index.md",
    "Get Started" => "start.md",
    "Tutorials" => "tutorials.md",
    "Main Features" => "main.md",
    "API documentation" => "API.md",
    "Miscellaneous" => "detail.md",
    "Video Examples" => "videos.md",
    ]
)

lst=("basics.jl","geography.jl","vectors.jl")
pth=("tutorials","tutorials","tutorials")

for ii in 1:length(lst)
    i=lst[ii]
    fil_in=joinpath(@__DIR__,"..", "examples",i)
    fil_out=joinpath(@__DIR__,"build", pth[ii] ,i[1:end-2]*"html")
    PlutoSliderServer.export_notebook(fil_in)
    mv(fil_in[1:end-2]*"html",fil_out)
    cp(fil_in,fil_out[1:end-4]*"jl")
end

deploydocs(
    repo = "github.com/JuliaClimate/MeshArrays.jl",
)
