## Module associated with the MeshArrays package

module MeshArrays

using Pkg
thistoml=joinpath(dirname(pathof(MeshArrays)), "..", "Project.toml")
thisversion=Pkg.TOML.parsefile(thistoml)["version"]

include("Types.jl");
include("Grids.jl");
include("Operations.jl");
include("Exchanges.jl");
include("ReadWrite.jl");
include("Solvers.jl");
include("Interpolation.jl");

export AbstractMeshArray, MeshArray, gcmgrid
export exchange, gradient, convergence, smooth, mask
export GridSpec, GridLoad, GridOfOnes, GridAddWS!
export Tiles, InterpolationFactors, knn
export ScalarPotential, VectorPotential, ThroughFlow
export StereographicProjection, LatitudeCircles
export PolygonAngle, QuadArrays, QuadCoeffs, ParaCoeffs

#The following exch_UV differs from normal exchange; incl. exch_UV_N
export exch_UV

end # module
