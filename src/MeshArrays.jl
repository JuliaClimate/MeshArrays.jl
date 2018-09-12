## Module associated with the MeshArrays package

module MeshArrays

using Printf

include("gcmfaces_type.jl");
include("gcmfaces_grids.jl");
include("gcmfaces_calc.jl");
include("gcmfaces_exch.jl");
include("gcmfaces_convert.jl");
include("gcmfaces_IO.jl");
include("gcmfaces_demo.jl");

export gcmfaces, exchange, gradient, smooth, mask, fsize
export GCMGridSpec, GCMGridLoad, GCMGridOnes
export demo1, demo2
#The following functions rely on grid specs; currently via global vars.
export read_bin, convert2array, convert2gcmfaces
#The following exch_UV differs from normal exchange; incl. exch_UV_N.
export exch_UV
#The following codes add dependencies to Plots & NetCDF.
#include("gcmfaces_plot.jl");
#include("gcmfaces_nctiles.jl");

end # module
