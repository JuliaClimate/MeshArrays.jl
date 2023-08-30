## Module associated with the MeshArrays package

module MeshArrays

import Pkg
thistoml=joinpath(dirname(pathof(MeshArrays)), "..", "Project.toml")
thisversion=Pkg.TOML.parsefile(thistoml)["version"]

using LazyArtifacts
p=dirname(pathof(MeshArrays))
artifact_toml = joinpath(p, "../Artifacts.toml")

GRID_LLC90_hash = artifact_hash("GRID_LLC90", artifact_toml)
GRID_LLC90 = joinpath(artifact_path(GRID_LLC90_hash)*"/","GRID_LLC90-1.1/")
GRID_LLC90_download() = artifact"GRID_LLC90"
GRID_LL360_hash = artifact_hash("GRID_LL360", artifact_toml)
GRID_LL360 = joinpath(artifact_path(GRID_LL360_hash)*"/","GRID_LL360-1.0/")
GRID_LL360_download() = artifact"GRID_LL360"
GRID_CS32_hash = artifact_hash("GRID_CS32", artifact_toml)
GRID_CS32 = joinpath(artifact_path(GRID_CS32_hash)*"/","GRID_CS32-1.1/")
GRID_CS32_download() = artifact"GRID_CS32"

include("Types.jl")
include("Grids.jl")
include("GridPaths.jl")
include("Operations.jl")
include("Exchanges.jl")
include("ReadWrite.jl")
include("Solvers.jl")
include("ReIndexing.jl")
include("Interpolation.jl")
include("VerticalDimension.jl")

export AbstractMeshArray, MeshArray, gcmgrid, varmeta
export GridSpec, GridLoad, GridLoadVar, UnitGrid, simple_periodic_domain
export exchange, Tiles, Tiles!, Interpolate, InterpolationFactors, knn
#The following exch_UV differs from normal exchange; incl. exch_UV_N
export exch_UV
export nansum, nanmean, nanmax, nanmin

#export InnerArray, OuterArray
#export smooth, mask
#export GridOfOnes 
#export GridAddWS!

export ScalarPotential, VectorPotential
export UVtoUEVN, UVtoTransport, UVtoTransport!
export gradient, curl, convergence
export LatitudeCircles, Transect, ThroughFlow
export StereographicProjection, isosurface

#export location_is_out, NeighborTileIndices_dpdo, NeighborTileIndices_cs, RelocationFunctions_cs
#export update_location_cs!, update_location_llc!, update_location_dpdo!

function unzip end
function read_polygons end
function read_shp end
function read_json end
function read_JLD2 end
function write_JLD2 end
function download_polygons end
function ocean_sections end
function one_section end
function ocean_basins end
function interpolation_setup end
function land_mask end

function examples_plot end

end # module
