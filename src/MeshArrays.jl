## Module associated with the MeshArrays package

module MeshArrays

using Pkg, Pkg.Artifacts
thistoml=joinpath(dirname(pathof(MeshArrays)), "..", "Project.toml")
thisversion=Pkg.TOML.parsefile(thistoml)["version"]

p=dirname(pathof(MeshArrays))
artifact_toml = joinpath(p, "../Artifacts.toml")

GRID_LLC90_hash = artifact_hash("GRID_LLC90", artifact_toml)
GRID_LLC90 = joinpath(artifact_path(GRID_LLC90_hash)*"/","GRID_LLC90-1.1/")
GRID_LL360_hash = artifact_hash("GRID_LL360", artifact_toml)
GRID_LL360 = joinpath(artifact_path(GRID_LL360_hash)*"/","GRID_LL360-1.0/")
GRID_CS32_hash = artifact_hash("GRID_CS32", artifact_toml)
GRID_CS32 = joinpath(artifact_path(GRID_CS32_hash)*"/","GRID_CS32-1.1/")

include("Types.jl");
include("Grids.jl");
include("Operations.jl");
include("Exchanges.jl");
include("ReadWrite.jl");
include("Solvers.jl");
include("Interpolation.jl");

export AbstractMeshArray, MeshArray, InnerArray, OuterArray, varmeta
export gcmgrid, exchange, gradient, convergence, smooth, mask
export simple_periodic_domain, GridSpec, GridLoad, GridOfOnes, GridAddWS!
export Tiles, Interpolate, InterpolationFactors, knn
export ScalarPotential, VectorPotential, ThroughFlow
export StereographicProjection, LatitudeCircles

#The following exch_UV differs from normal exchange; incl. exch_UV_N
export exch_UV

end # module
