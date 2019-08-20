## Module associated with the MeshArrays package

module MeshArrays

#using Printf

include("Types.jl");
include("gcmfaces_grids.jl");
include("gcmfaces_calc.jl");
include("gcmfaces_exch.jl");
include("gcmfaces_convert.jl");
include("gcmfaces_IO.jl");
include("gcmfaces_demo.jl");

export AbstractMeshArray, MeshArray, gcmgrid
export exchange, gradient, convergence, smooth, mask
export GCMGridSpec, GCMGridLoad, GCMGridOnes
export nFacesEtc, fijind, findtiles, LatitudeCircles, ThroughFlow
#The following exch_UV differs from normal exchange; incl. exch_UV_N.
export exch_UV
#The following codes add dependencies to Plots & NetCDF.
#include("gcmfaces_plot.jl");
#include("gcmfaces_nctiles.jl");

end # module
