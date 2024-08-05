
include("ERA5_module.jl")
import Main.ERA5_to_ECCO4: loop_over_years

v=1 #parse(Int64, ARGS[1])
S=loop_over_years(v)

if false
    using MeshArrays, CairoMakie, DataDeps, JLD2
    γ=GridSpec("LatLonCap",MeshArrays.GRID_LLC90)
    λ=interpolation_setup()
    heatmap(read(S.y-S.yi,γ),interpolation=λ)
end