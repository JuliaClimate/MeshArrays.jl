## Additional Detail

Functions like `GCMGridSpec("LLC90")` return a `gcmgrid` struct that just contains the basic specifications of a global grid -- not the grid itself; just a few parameters, ranges, and possibly a file location. A `gcmgrid` is embeded in each `MeshArray` instance for which it provides a blueprint. It specifies how an array collection forms a global mesh and allows `exchange` to dispatch to the appropriate method. 

Various configurations that are commonly used in `Earth System Models` are implemented using the concrete type called `MeshArray`. This type is in fact an alias for other types that can thus be used interchangeably (initially: `gcmfaces` or `gcmarray`).

Within a `MeshArray`, a whole Earth System Model grid is represented as an array of elementary arrays. Each one typically represents a subdomain mesh. For example, `gcmarray` represents a global 2D state variable `x` as a column array of 2D ragged arrays `x.f`. `MeshArrays.demo1` illustrates how the `MeshArray` data structures operate via standard functions. In brief, a `MeshArray` is used like any other Array type.

_Background:_ `MeshArrays.jl` is rooted in a [previous Matlab / Octave package](https://gcmfaces.readthedocs.io/en/latest/), `gcmfaces` introduced in [Forget et al., 2015](http://www.geosci-model-dev.net/8/3071/2015/) (`doi:10.5194/gmd-8-3071-2015`). `GCM` is an acronym for General Circulation Model, or Global Climate Model (see [this wikipedia entry](https://en.wikipedia.org/wiki/General_circulation_model)), and `faces` can mean meshes, arrays, facets, or subdomains (`x.f` in a `MeshArray ` instance `x`).

