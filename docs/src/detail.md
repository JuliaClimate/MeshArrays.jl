# Miscellaneous

### Details

Functions like `GridSpec("LLC90")` return a `gcmgrid` struct that contains the basic specification of a global grid. This is not the grid itself -- just a few parameters, ranges, and possibly a path to grid files. A `gcmgrid` is embeded in each `MeshArray` instance for which it provides a blueprint. It specifies how an array collection forms a global mesh and allows e.g. the `exchange` function to dispatch to the appropriate method. 

Various configurations that are commonly used in `Earth System Models` are readily implemented using the concrete type called `MeshArray`. This type is in fact an alias for more specific types that can be used interchangeably via `MeshArray` (initially: `gcmfaces` or `gcmarray`).

Within a `MeshArray`, a whole Earth System Model grid is represented as an array of elementary arrays. Each one of these represents a subdomain. For example, a `gcmarray` instance for one Earth map `x` has a column array `x.f` of elementary 2D arrays of various sizes. `demo1` illustrates how one easily operates `MeshArray` structs via standard and specialized functions. In brief, a `MeshArray` should be used just like a common `Array`.

### Background

`MeshArrays.jl` is rooted in a [Matlab / Octave package](https://gcmfaces.readthedocs.io/en/latest/) called `gcmfaces`, which was introduced in [Forget et al., 2015](http://www.geosci-model-dev.net/8/3071/2015/) (`doi:10.5194/gmd-8-3071-2015`). `GCM` is an acronym for [General Circulation Model](https://en.wikipedia.org/wiki/General_circulation_model), or Global Climate Model, and `faces` can mean meshes, arrays, facets, or subdomains (these are the elements of `x.f` in a `MeshArray ` instance `x`).
