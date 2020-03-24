## Module associated with the MeshArrays package

module MeshArrays

using Pkg
thisversion=Pkg.TOML.parsefile(joinpath(dirname(pathof(MeshArrays)), "..", "Project.toml"))["version"]

include("Types.jl");
include("Grids.jl");
include("Operations.jl");
include("Exchanges.jl");
include("ReadWrite.jl");
include("Solvers.jl");

export AbstractMeshArray, MeshArray, gcmgrid
export exchange, gradient, convergence, smooth, mask
export GridSpec, GridLoad, GridOfOnes, GridAddWS!
export TileMap, LatitudeCircles, ThroughFlow
export ScalarPotential, VectorPotential

#The following exch_UV differs from normal exchange; incl. exch_UV_N
export exch_UV

end # module
