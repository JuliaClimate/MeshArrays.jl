# Internals

Functions like `GridSpec("LLC90")` return a `gcmgrid` struct that contains the basic specification of a global grid. This is not the grid itself -- just a few parameters, ranges, and possibly a path to grid files. A `gcmgrid` is embeded in each `MeshArray` instance for which it provides a blueprint. It specifies how an array collection forms a global mesh and allows e.g. the `exchange` function to dispatch to the appropriate method. 

Various configurations that are commonly used in `Earth System Models` are readily implemented using the concrete type called `MeshArray`. This type is in fact an alias for more specific types that can be used interchangeably via `MeshArray` (initially: `gcmfaces` or `gcmarray`).

Within a `MeshArray`, a whole Earth System Model grid is represented as an array of elementary arrays. Each one of these represents a subdomain. For example, a `gcmarray` instance for one Earth map `x` has a column array `x.f` of elementary 2D arrays of various sizes. 

The [basics tutorial](@ref id_Basics) illustrates how standard operations apply to `MeshArray` like as they do to common `Array`. More specialized functions and distinctive features, such as domain decomposition or plotting maps, are demo'ed in the [geography tutorial](@ref id_Geography) and [vector tutorial](@ref id_Vectors).

# Background

The origin of `MeshArrays.jl` is rooted in a [Matlab / Octave package](https://gcmfaces.readthedocs.io/en/latest/) called `gcmfaces`, which was introduced in [Forget et al., 2015](http://www.geosci-model-dev.net/8/3071/2015/) (`doi:10.5194/gmd-8-3071-2015`). `GCM` is an acronym for [General Circulation Model](https://en.wikipedia.org/wiki/General_circulation_model), or Global Climate Model, and `faces` can be equivalent to meshes, arrays, facets, or subdomains (these are the elements of `x.f` in a `MeshArray ` instance `x`).
