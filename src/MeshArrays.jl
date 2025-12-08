## Module associated with the MeshArrays package

module MeshArrays

using LazyArtifacts
p=dirname(pathof(MeshArrays))
artifact_toml = joinpath(p, "../Artifacts.toml")

import GeoInterface as GI

include("extensions.jl")
include("types/main.jl")
include("ReadWrite.jl")
include("grids/main.jl")
include("GridPaths.jl")
include("Operations.jl")
include("Exchanges.jl")
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
export GridSpec, GridLoad, GridLoadVar
export exchange, Tiles, Tiles!
export Interpolate, InterpolationFactors, knn, interpolation_setup
#The following exch_UV differs from normal exchange; incl. exch_UV_N
export exch_UV
export nansum, nanmean, nanmax, nanmin
export demo, land_mask, edge_path
export Integration

#export InnerArray, OuterArray
#export smooth, mask

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
