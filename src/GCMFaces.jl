# GCMFaces.jl
#
#	First Draft Implementation
#
# gaelforget (https://github.com/gaelforget/gcmfaces_jl)
# Julia 0.6.2
# Created: 02.01.18
# Last Edit: 02.01.18

#Activate type definitions and methods

__precompile__()

module GCMFaces

# the following is needed for embedded calls to e.g. heatmap in qwckplot
using Plots
# the following is needed for embedded calls to e.g. ncread in read_nctiles
using NetCDF

include("gcmfaces_type.jl");
include("gcmfaces_exch.jl");
include("gcmfaces_convert.jl");
include("gcmfaces_plot.jl");
include("gcmfaces_IO.jl");
include("gcmfaces_grids.jl");

export gcmfaces, exchange, exch_T_N, convert2array, qwckplot
export GCMGridSpec, GCMGridLoad, read_nctiles, read_bin

end # module
