## Module associated with the MeshArrays package

module MeshArrays

using LazyArtifacts
p=dirname(pathof(MeshArrays))
artifact_toml = joinpath(p, "../Artifacts.toml")

function read_polygons end
function read_shp end
function read_json end
function read_JLD2 end
function write_JLD2 end
function ProjAxis end
function grid_lines! end

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
include("demo.jl")
include("Integration.jl")

export AbstractMeshArray, MeshArray, gcmgrid, varmeta, gridpath, gridmask
export GridSpec, GridLoad, GridLoadVar, UnitGrid, simple_periodic_domain
export exchange, Tiles, Tiles!, Interpolate, InterpolationFactors, knn, interpolation_setup
#The following exch_UV differs from normal exchange; incl. exch_UV_N
export exch_UV
export nansum, nanmean, nanmax, nanmin
export demo, land_mask
export Integration

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

function MA_datadep end
export MA_datadep

function plot_examples end; export plot_examples
examples_plot=plot_examples

end # module
