## Additional Detail

Functions like `GCMGridSpec("LLC90")` return a `gcmgrid` struct that just contains the basic specifications of a global grid -- not the grid itself; just a few parameters, ranges, and possibly a file location. A `gcmgrid` is embeded in each `MeshArray` objects for which they serve as a blueprint. The `gcmgrid` determines how an array collection forms a global mesh and allows `exchange` to dispatch to the appropriate method. 

Various configurations that are commonly used in `Earth System Models` are implemented using the concrete type called `MeshArray`. This type is used as an alias for other types that can thus used interchangeably. Initially, `gcmfaces` or `gcmarray` are the two available options for `MeshArray`. 

Within a `MeshArray`, a whole Earth System Model grid is represented as an array of elemental arrays. Each one of these typically represents a subdomain mesh. For example, `gcmarray` represents a climate state variable `x` as an array of 2D ragged arrays `x.f`. `MeshArrays.demo1` illustrates the `MeshArray` data structures and basic methods. In brief, a `MeshArray` is used simply like any other Array type.

_Background:_ `MeshArrays.jl` is rooted in a [previous Matlab / Octave package](https://gcmfaces.readthedocs.io/en/latest/), `gcmfaces` introduced in [Forget et al., 2015](http://www.geosci-model-dev.net/8/3071/2015/) (`doi:10.5194/gmd-8-3071-2015`). `GCM` is an acronym for General Circulation Model, or Global Climate Model (see [this wikipedia entry](https://en.wikipedia.org/wiki/General_circulation_model)), and `faces` can mean meshes, arrays, facets, or subdomains (`x.f` in a `MeshArray ` instance `x`).

