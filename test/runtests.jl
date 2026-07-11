include("setup.jl")

# Select testsets via env var, e.g.: JULIA_TESTSETS=transport,gcmfaces julia runtests.jl
# Omit the env var (or set to "all") to run everything.
const _RUN = split(get(ENV, "JULIA_TESTSETS", "all"), ",")
include_maybe(name) = ("all" in _RUN || name in _RUN) && include("testsets/$name.jl")

include_maybe("mesharray_basic")
include_maybe("vertical_dim")
include_maybe("regional_integration")
include_maybe("transport")
include_maybe("gcmfaces")
include_maybe("unitgrid")
include_maybe("nemo_grid")
include_maybe("gridspec")
include_maybe("plotting_makie")
include_maybe("nanmath")
include_maybe("plotting_basemap")
include_maybe("polygon_ops")
include_maybe("doctests")
