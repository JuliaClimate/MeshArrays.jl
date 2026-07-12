# Internals

Functions like `GridSpec("LLC90")` return a `gcmgrid` struct that contains the basic specification of a global grid. This is not the grid itself -- just a few parameters, ranges, and possibly a path to grid files. A `gcmgrid` is embeded in each `MeshArray` instance for which it provides a blueprint. It specifies how an array collection forms a global mesh and allows e.g. the `exchange` function to dispatch to the appropriate method. 

Various configurations that are commonly used in `Earth System Models` are readily implemented using the concrete type called `MeshArray`. This type is in fact an alias for more specific types that can be used interchangeably via `MeshArray` (initially: `gcmfaces` or `gcmarray`).

Within a `MeshArray`, a whole Earth System Model grid is represented as an array of elementary arrays. Each one of these represents a subdomain. For example, a `gcmarray` instance for one Earth map `x` has a column array `x.f` of elementary 2D arrays of various sizes. 

The [basics tutorial](@ref id_Basics) illustrates how standard operations apply to `MeshArray` like as they do to common `Array`. More specialized functions and distinctive features, such as domain decomposition or plotting maps, are demo'ed in the [geography tutorial](@ref id_Geography) and [vector tutorial](@ref id_Vectors).

## Exchange Methods

Two families of methods connect grid cells across subdomain boundaries:

**Finding neighbors** — given a position `(i, j, face)`, locate the corresponding cell in an adjacent subdomain. The top-level entry point is `update_location!`, which dispatches to grid-specific implementations:
- `update_location_cs!` — cubed-sphere and LLC grids, via `RelocationFunctions_cs`
- `update_location_PeriodicDomain!` — periodic tiled grids, via `NeighborTileIndices_PeriodicDomain`

**Adding halo rows/columns** — `exchange` (public API) calls `exchange_main`, which wraps each face with extra rows/columns copied from neighboring subdomains, returning a `MeshArray_wh` (with halo). The per-topology implementations are:
- `exch_T_N_cs` — cubed-sphere / LLC, using `exch_cs_target` and `exch_cs_sources` to compute source/target index ranges for each face
- `exch_T_N_PeriodicChannel`, `exch_T_N_PeriodicDomain` — periodic topologies

The [Finding and Adding Neighbors](dev/exchange_methods.html) developer notebook demonstrates these methods interactively across several grid configurations.

# Background

The origin of `MeshArrays.jl` is rooted in a [Matlab / Octave package](https://gcmfaces.readthedocs.io/en/latest/) called `gcmfaces`, which was introduced in [Forget et al., 2015](http://www.geosci-model-dev.net/8/3071/2015/) (`doi:10.5194/gmd-8-3071-2015`). `GCM` is an acronym for [General Circulation Model](https://en.wikipedia.org/wiki/General_circulation_model), or Global Climate Model, and `faces` can be equivalent to meshes, arrays, facets, or subdomains (these are the elements of `x.f` in a `MeshArray ` instance `x`).

# Testing

The test suite lives in `test/` and is split into per-topic files under `test/testsets/`. The shared setup (package imports, dataset downloads, `Demos.jl`) is in `test/setup.jl`. The top-level `test/runtests.jl` runs all testsets and is used by `Pkg.test`.

The `JULIA_TESTSETS` environment variable selects a subset when running `runtests.jl` directly:

```
JULIA_TESTSETS=mesharray_basic,transport julia --project=test test/runtests.jl
```

Currently available testset names (roughly in order of importance): `mesharray_basic`, `transport`, `regional_integration`, `interpolation`, `vertical_dim`, `polygon_ops`, `plotting_makie`, `plotting_basemap`, `unitgrid`, `gridspec`, `nemo_grid`, `nanmath`, `datasets`, `doctests`.

A convenience script (not part of the package) can also be used from the repository root. It activates the test environment and forwards command-line arguments as the testset filter:

```julia
import Pkg
Pkg.activate(joinpath(@__DIR__, "MeshArrays.jl", "test"))

!isempty(ARGS) && (ENV["JULIA_TESTSETS"] = join(ARGS, ","))

include(joinpath(@__DIR__, "MeshArrays.jl", "test", "runtests.jl"))
```

```
julia run_tests.jl                              # all testsets
julia run_tests.jl transport                    # one testset
julia run_tests.jl mesharray_basic transport    # subset
```
