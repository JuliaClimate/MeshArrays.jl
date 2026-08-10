using Test, Documenter, Suppressor, MeshArrays, CairoMakie
import DataDeps, JLD2, Shapefile, GeoJSON, Proj, GeometryOps
import MITgcm, Climatology
import MeshArrays: GI, NEMO_GRID, Grids_simple

MeshArrays.Dataset("GRID_LL360")
MeshArrays.Dataset("GRID_LLC90")
MeshArrays.Dataset("GRID_LLC270")
MeshArrays.Dataset("GRID_CS32")

pol_shp=MeshArrays.Dataset("countries_shp1")
pol_json=MeshArrays.Dataset("oceans_geojson1")

p=dirname(pathof(MeshArrays))
include(joinpath(p,"../examples/Demos.jl"))
