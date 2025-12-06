## Module associated with the MeshArrays package

module MeshArrays

using LazyArtifacts
p=dirname(pathof(MeshArrays))
artifact_toml = joinpath(p, "../Artifacts.toml")

import GeoInterface as GI

function read_shp end
function read_json end
function read_jld2 end
function write_jld2 end
function ProjAxis end
function grid_lines! end
function within_pol end

include("types/main.jl")
include("grids/main.jl")
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
include("Datasets.jl")
include("Polygons.jl")

export AbstractMeshArray, MeshArray, MeshArray_wh
export gcmgrid, varmeta, gridpath, gridmask
export GridSpec, GridLoad, GridLoadVar, Grids_simple, NEMO_GRID
export exchange, Tiles, Tiles!, Interpolate, InterpolationFactors, knn, interpolation_setup
#The following exch_UV differs from normal exchange; incl. exch_UV_N
export exch_UV
export nansum, nanmean, nanmax, nanmin
export demo, land_mask, edge_path
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

function mydatadep end

MA_datadep=mydatadep
export MA_datadep

function plot_examples end; export plot_examples
examples_plot=plot_examples

end # module
