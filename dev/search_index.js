var documenterSearchIndex = {"docs":
[{"location":"videos/#Video-Examples","page":"Video Examples","title":"Video Examples","text":"","category":"section"},{"location":"videos/","page":"Video Examples","title":"Video Examples","text":"The videos below introduce MeshArrays.jl and present two use cases :","category":"page"},{"location":"videos/","page":"Video Examples","title":"Video Examples","text":"introduction of MeshArrays.jl at JuliaCon 2018, for context.\ninteractive analysis of gridded ocean model output for ECCO.\nsimulation of plastic particles following ocean currents.","category":"page"},{"location":"videos/#JuliaCon-2018","page":"Video Examples","title":"JuliaCon 2018","text":"","category":"section"},{"location":"videos/","page":"Video Examples","title":"Video Examples","text":"MeshArrays.jl was first introduced at JuliaCon in 2018. This video provides some motivation and context for the project.","category":"page"},{"location":"videos/","page":"Video Examples","title":"Video Examples","text":"(Image: JuliaCon-2018 presentation)","category":"page"},{"location":"videos/#Interactive-Model-Analysis","page":"Video Examples","title":"Interactive Model Analysis","text":"","category":"section"},{"location":"videos/","page":"Video Examples","title":"Video Examples","text":"In this demo of interactive visualization of ocean variables, the notebook provides various options for choosing variables and viewing them. Plots etc react to these user selections as illustrated in the video. Code is based on  MeshArrays.jl (grids, arrays), OceanStateEstimation.jl (velocity fields, etc), Pluto.jl (notebook), Makie.jl (plotting), and other open source packages.","category":"page"},{"location":"videos/","page":"Video Examples","title":"Video Examples","text":"(Image: ocean state estimate analysis)","category":"page"},{"location":"videos/","page":"Video Examples","title":"Video Examples","text":"Another method for interacting with data sets is to use GLMakie.jl or WGLMakie.jl. For example, Tyler.jl can be used to explore high-resolution model output.","category":"page"},{"location":"videos/","page":"Video Examples","title":"Video Examples","text":"(Image: high resolution model analysis)","category":"page"},{"location":"videos/#Particle-Tracking-and-Modeling","page":"Video Examples","title":"Particle Tracking and Modeling","text":"","category":"section"},{"location":"videos/","page":"Video Examples","title":"Video Examples","text":"Here we visualize a simulation of particles moving at a fixed depth in the Ocean (300m depth). This uses IndividualDisplacements.jl, and MeshArrays.jl underneath, to simulate particle trajectories. ","category":"page"},{"location":"videos/","page":"Video Examples","title":"Video Examples","text":"(Image: simulated particle movie (300m))","category":"page"},{"location":"videos/","page":"Video Examples","title":"Video Examples","text":"More examples like this, using Julia to run models, and related work in JuliaOcean are available in the longer video below. ","category":"page"},{"location":"videos/","page":"Video Examples","title":"Video Examples","text":"(Image: Modeling Marine Ecosystems)","category":"page"},{"location":"tutorials/#Tutorials","page":"Notebook Tutorials","title":"Tutorials","text":"","category":"section"},{"location":"tutorials/#id_Basics","page":"Notebook Tutorials","title":"Basics","text":"","category":"section"},{"location":"tutorials/","page":"Notebook Tutorials","title":"Notebook Tutorials","text":"The basics tutorial (code link) illustrates how the MeshArrays.jl data structures let us write generic code readily applicable to whole families of grids. ","category":"page"},{"location":"tutorials/","page":"Notebook Tutorials","title":"Notebook Tutorials","text":"It focuses on a global workflow (smoothing) that requires communication across the entire gridded domain – a key feature provided by MeshArrays.jl.  The same workflow is repeated three times, for different grid configurations commonly used in numerical models.","category":"page"},{"location":"tutorials/","page":"Notebook Tutorials","title":"Notebook Tutorials","text":"Grid scale noise Smoothed noise\n(Image: raw) (Image: smooth)","category":"page"},{"location":"tutorials/#id_Geography","page":"Notebook Tutorials","title":"Geography","text":"","category":"section"},{"location":"tutorials/","page":"Notebook Tutorials","title":"Notebook Tutorials","text":"The geography tutorial (code link) deals with interpolation, projection, and vizualization of gridded fields in geographic coordinates. ","category":"page"},{"location":"tutorials/","page":"Notebook Tutorials","title":"Notebook Tutorials","text":"Ocean Depth Map Great Circle Path\n(Image: map projection) (Image: great circle)","category":"page"},{"location":"tutorials/#id_Vectors","page":"Notebook Tutorials","title":"Vector Fields","text":"","category":"section"},{"location":"tutorials/","page":"Notebook Tutorials","title":"Notebook Tutorials","text":"The vectors tutorial (code link) illustrates how MeshArrays.jl represents gridded, vector fields. This enables  analyses of various quantities, like heat, flow with oceanic currents and atmospheric winds within the climate system. ","category":"page"},{"location":"tutorials/","page":"Notebook Tutorials","title":"Notebook Tutorials","text":"Horizontal Streamfunction Vertical Streamfunction\n(Image: raw) (Image: smooth)","category":"page"},{"location":"detail/#Internal-Details","page":"Miscellaneous","title":"Internal Details","text":"","category":"section"},{"location":"detail/","page":"Miscellaneous","title":"Miscellaneous","text":"Functions like GridSpec(\"LLC90\") return a gcmgrid struct that contains the basic specification of a global grid. This is not the grid itself – just a few parameters, ranges, and possibly a path to grid files. A gcmgrid is embeded in each MeshArray instance for which it provides a blueprint. It specifies how an array collection forms a global mesh and allows e.g. the exchange function to dispatch to the appropriate method. ","category":"page"},{"location":"detail/","page":"Miscellaneous","title":"Miscellaneous","text":"Various configurations that are commonly used in Earth System Models are readily implemented using the concrete type called MeshArray. This type is in fact an alias for more specific types that can be used interchangeably via MeshArray (initially: gcmfaces or gcmarray).","category":"page"},{"location":"detail/","page":"Miscellaneous","title":"Miscellaneous","text":"Within a MeshArray, a whole Earth System Model grid is represented as an array of elementary arrays. Each one of these represents a subdomain. For example, a gcmarray instance for one Earth map x has a column array x.f of elementary 2D arrays of various sizes. ","category":"page"},{"location":"detail/","page":"Miscellaneous","title":"Miscellaneous","text":"The basics tutorial illustrates how standard operations apply to MeshArray like as they do to common Array. More specialized functions and distinctive features, such as domain decomposition or plotting maps, are demo'ed in the geography tutorial and vector tutorial.","category":"page"},{"location":"detail/#Background","page":"Miscellaneous","title":"Background","text":"","category":"section"},{"location":"detail/","page":"Miscellaneous","title":"Miscellaneous","text":"The origin of MeshArrays.jl is rooted in a Matlab / Octave package called gcmfaces, which was introduced in Forget et al., 2015 (doi:10.5194/gmd-8-3071-2015). GCM is an acronym for General Circulation Model, or Global Climate Model, and faces can be equivalent to meshes, arrays, facets, or subdomains (these are the elements of x.f in a MeshArray instance x).","category":"page"},{"location":"main/#Main-Features","page":"Main Features","title":"Main Features","text":"","category":"section"},{"location":"main/","page":"Main Features","title":"Main Features","text":"MeshArray, gcmgrid, varmeta\nfull Earth grid examples (C-grids)\nvector fields, transports, budgets\ninterpolation, distances, collocation\nvisualization (via Makie extension)\nparticle tracking via IndividualDisplacements.jl","category":"page"},{"location":"main/#Summary","page":"Main Features","title":"Summary","text":"","category":"section"},{"location":"main/","page":"Main Features","title":"Main Features","text":"The MeshArray type is a sub-type of AbstractArray with an outer array where each element is itself a 2D inner array. By default, outer and inner arrays are of all of the standard Array type. However, this setup potentially allows different choices for the outer and inner arrays – for example DistributedArrays and AxisArrays, respectively, could be an option. MeshArrays.jl thus provides a simple but general solution to analyze or e.g. simulate climate system variables. By default the MeshArray type is an alias to the MeshArrays.gcmarray type.","category":"page"},{"location":"main/","page":"Main Features","title":"Main Features","text":"The internals of a MeshArray are regulated by its gcmgrid – a struct containing just a few index ranges, array size specifications, and inter-connection rules. A second  lightweight struct, varmeta, contains the MeshArray variable name, its unit, time, and position in grid space. A general approach like this is useful because climate models often involve advanced domain decompositions (see, e.g., Grids), and many variables, which can put a burden on users. ","category":"page"},{"location":"main/","page":"Main Features","title":"Main Features","text":"Encoding the grid specification inside the MeshArray data type allows user to manipulate MeshArrays just like they would manipulate Arrays without having to invoke model grid details explicitely. In addition, the provided exchange methods readily transfer data between connected subdomains to extend them at the sides. This makes it easy to compute e.g. partial derivatives and related operators like gradients, curl, or divergences over subdomain edges as often needed for precise computation of transports, budgets, etc using climate model output (see, e.g., Tutorials).","category":"page"},{"location":"main/#Data-Structures","page":"Main Features","title":"Data Structures","text":"","category":"section"},{"location":"main/","page":"Main Features","title":"Main Features","text":"The elements of a MeshArray /  are arrays. These elementary arrays typically represent subdomains inter-connected at their edges. The organization and connections between subdomains is determined by a user-specified gcmgrid which is embeded inside each MeshArrays.gcmarray instance. ","category":"page"},{"location":"main/","page":"Main Features","title":"Main Features","text":"Interpolate can be used to interpolate a MeshArray to any location (i.e. arbitrary longitude, latitude pair). Exchange methods transfer data between neighboring arrays to extend computational subdomains – this is often needed in analyses of climate or ocean model output. ","category":"page"},{"location":"main/","page":"Main Features","title":"Main Features","text":"The current default for MeshArray is the gcmarray type, with various examples provided in the Tutorials.","category":"page"},{"location":"main/","page":"Main Features","title":"Main Features","text":"One of the examples is based on a grid known as LatLonCap where each global map is associated with 5 subdomains of different sizes. The grid has 50 depth levels. Such a MeshArray has a size of (5, 50) (see Tutorials).","category":"page"},{"location":"main/","page":"Main Features","title":"Main Features","text":"The underlying, MeshArray, data structure is:","category":"page"},{"location":"main/","page":"Main Features","title":"Main Features","text":"struct gcmarray{T, N} <: AbstractMeshArray{T, N}\n   grid::gcmgrid\n   meta::varmeta\n   f::Array{Array{T,2},N}\n   fSize::Array{NTuple{2, Int}}\n   fIndex::Array{Int,1}\n   version::String\nend","category":"page"},{"location":"main/","page":"Main Features","title":"Main Features","text":"A MeshArray generally behaves just like an Array and the broadcasting of operations has notably been customized so that it reaches elements of each elementary array (i.e. within f[i] for each index of f).","category":"page"},{"location":"main/","page":"Main Features","title":"Main Features","text":"In addition, Mesharray specific functions like exchange can alter the internal structure of a MeshArray by adding rows and columns at the periphery of subdomains. ","category":"page"},{"location":"main/#Embedded-Metadata","page":"Main Features","title":"Embedded Metadata","text":"","category":"section"},{"location":"main/","page":"Main Features","title":"Main Features","text":"A MeshArrays.gcmarray includes a gcmgrid specification which can be constructed as outlined below.","category":"page"},{"location":"main/","page":"Main Features","title":"Main Features","text":"struct gcmgrid\n  path::String\n  class::String\n  nFaces::Int\n  fSize::Array{NTuple{2, Int},1}\n  ioSize::Union{NTuple{2, Int},Array{Int64,2}}\n  ioPrec::Type\n  read::Function\n  write::Function\nend","category":"page"},{"location":"main/","page":"Main Features","title":"Main Features","text":"The grid class can be set to LatLonCap, CubeSphere, PeriodicChannel, or PeriodicDomain. For example, A PeriodicChannel (periodic in the x direction) of size 360 by 160, can be defined as follows.","category":"page"},{"location":"main/","page":"Main Features","title":"Main Features","text":"pth=MeshArrays.GRID_LL360\nclass=\"PeriodicChannel\"\nioSize=(360, 160)\nioPrec=Float32\n\nγ=gcmgrid(pth, class, 1, [ioSize], ioSize, ioPrec, read, write)","category":"page"},{"location":"main/","page":"Main Features","title":"Main Features","text":"Importantly, a gcmgrid does not contain any actual grid data – hence its memory footprint is minimal. Grid variables are instead read to memory only when needed e.g. as shown below. To make this easy, each gcmgrid includes a pair of read / write methods to allow for basic I/O at any time. These methods are typically specified by the user although defaults are provided. ","category":"page"},{"location":"main/","page":"Main Features","title":"Main Features","text":"using MeshArrays, Unitful\nγ=GridSpec(\"LatLonCap\",MeshArrays.GRID_LLC90)\nm=MeshArrays.varmeta(u\"m\",fill(0.5,2),missing,\"Depth\",\"Depth\")\nD=γ.read(γ.path*\"Depth.data\",MeshArray(γ,Float64;meta=m))","category":"page"},{"location":"main/","page":"Main Features","title":"Main Features","text":"The above commands define a MeshArray called D which is the one displayed at the top of this section. A definition of the varmeta structure is reported below. The position of a D point within its grid cell is given as x ∈ [0. 1.] in each direction.","category":"page"},{"location":"main/","page":"Main Features","title":"Main Features","text":"struct varmeta\n  unit::Union{Unitful.Units,Number,Missing}\n  position::Array{Float64,1}\n  time::Union{DateTime,Missing,Array{DateTime,1}}\n  name::String\n  long_name::String\nend","category":"page"},{"location":"main/","page":"Main Features","title":"Main Features","text":"(Image: OceanDepthMap)","category":"page"},{"location":"main/#Visualization,-Particles,-Transports","page":"Main Features","title":"Visualization, Particles, Transports","text":"","category":"section"},{"location":"main/","page":"Main Features","title":"Main Features","text":"A simple way to plot a MeshArray consists in using the Makie extension. ","category":"page"},{"location":"main/","page":"Main Features","title":"Main Features","text":"By default, for a MeshArray the heatmap command plots each elementary array separately. This is illustrated in Grids and in the Tutorials. If an interpolation scheme is provided then heatmap produces a global map instead. See the geography tutorial for examples. The JuliaClimate Notebooks provide more examples","category":"page"},{"location":"main/","page":"Main Features","title":"Main Features","text":"The vectors tutorial illustrates a common Earth System use case – using gridded flow fields to integrate transports, streamfunctions, budgets, etc. Particle trajectories are readily computed with IndividualDisplacements.jl when velocity fields are provided as MeshArrays. ","category":"page"},{"location":"main/","page":"Main Features","title":"Main Features","text":"The MITgcm.jl examples shows how MeshArrays.jl can ingest any standard grid from the MIT general circulation model.","category":"page"},{"location":"main/","page":"Main Features","title":"Main Features","text":"(Image: OceanMOC)","category":"page"},{"location":"start/#Get-Started","page":"Get Started","title":"Get Started","text":"","category":"section"},{"location":"start/","page":"Get Started","title":"Get Started","text":"To create your first MeshArray, open the Julia REPL and type:","category":"page"},{"location":"start/","page":"Get Started","title":"Get Started","text":"using MeshArrays\nC=MeshArray(randn(21,10))","category":"page"},{"location":"start/","page":"Get Started","title":"Get Started","text":"using CairoMakie\nheatmap(C,x=-10:10,y=1:10)","category":"page"},{"location":"start/#Grids","page":"Get Started","title":"Grids","text":"","category":"section"},{"location":"start/","page":"Get Started","title":"Get Started","text":"Three grids are available directly via this package the examples (GRID_LL360, GRID_CS32, and GRID_LLC90).","category":"page"},{"location":"start/","page":"Get Started","title":"Get Started","text":"Lat-lon Cube Sphere Lat-Lon-Cap\n(Image: Lat-Lon) (Image: Cube Sphere) (Image: Lat-Lon-Cap)","category":"page"},{"location":"start/","page":"Get Started","title":"Get Started","text":"In the example below we read and display GRID_CS32. This is a cubed-sphere grid where the EArth surface is projected on the six faces of a cube.","category":"page"},{"location":"start/","page":"Get Started","title":"Get Started","text":"using MeshArrays,CairoMakie\nγ=GridSpec(\"CubeSphere\",MeshArrays.GRID_CS32)\nΓ=GridLoad(γ)\nheatmap(Γ.YC,title=\"grid point latitudes\")","category":"page"},{"location":"start/#Install","page":"Get Started","title":"Install","text":"","category":"section"},{"location":"start/","page":"Get Started","title":"Get Started","text":"To install MeshArrays.jl and verify it works as expected, open the Julia REPL and type:","category":"page"},{"location":"start/","page":"Get Started","title":"Get Started","text":"using Pkg\nPkg.add(\"MeshArrays\")\nPkg.test(\"MeshArrays\")","category":"page"},{"location":"#MeshArrays.jl","page":"Home","title":"MeshArrays.jl","text":"","category":"section"},{"location":"","page":"Home","title":"Home","text":"MeshArrays.jl defines an array type that can contain / organize / distribute collections of inter-connected arrays as generally done in climate models (see Grids). Data structures defined in MeshArrays.jl can be used to simulate and analyze key variables of the climate system such as particles and transports.","category":"page"},{"location":"","page":"Home","title":"Home","text":"note: Note\nThe basics tutorial provides a walk through of core features. The geography tutorial deals with interpolation and mapping. The vectors tutorial demonstrates functionalities related to vector fields and transport processes. The JuliaClimate Notebooks link to a collection of use cases.","category":"page"},{"location":"API/#API-documentation","page":"API documentation","title":"API documentation","text":"","category":"section"},{"location":"API/#1.-Data-Structures","page":"API documentation","title":"1. Data Structures","text":"","category":"section"},{"location":"API/","page":"API documentation","title":"API documentation","text":"By default, the MeshArray type is an alias to MeshArrays.gcmarray.","category":"page"},{"location":"API/","page":"API documentation","title":"API documentation","text":"AbstractMeshArray\nMeshArrays.gcmarray\ngcmgrid\nvarmeta","category":"page"},{"location":"API/#MeshArrays.AbstractMeshArray","page":"API documentation","title":"MeshArrays.AbstractMeshArray","text":"AbstractMeshArray{T, N}\n\nSubtype of AbstractArray{T, N}\n\n\n\n\n\n","category":"type"},{"location":"API/#MeshArrays.gcmarray","page":"API documentation","title":"MeshArrays.gcmarray","text":"gcmarray{T, N, AT}\n\ngcmarray data structure. Available constructors:\n\ngcmarray{T,N,AT}(grid::gcmgrid,meta::varmeta,f::Array{AT,N},\n         fSize::Array{NTuple{N, Int}},fIndex::Array{Int,1},v::String)\n\ngcmarray(grid::gcmgrid,f::Array{Array{T,2},N}) where {T,N}\ngcmarray(grid::gcmgrid,f::Array{Array{T,N},1}) where {T,N}\n\ngcmarray(grid::gcmgrid,fSize::Array{NTuple{N, Int}},fIndex::Array{Int,1})\ngcmarray(<same as above>,n3::Int)\ngcmarray(<same as above>,n3::Int,n4::Int)\n\ngcmarray(grid::gcmgrid)\ngcmarray(grid::gcmgrid,::Type{T})\ngcmarray(grid::gcmgrid,::Type{T},n3::Int)\ngcmarray(grid::gcmgrid,::Type{T},n3::Int,n4::Int)\n\ngcmarray(A::Array{T,N};meta::varmeta=defaultmeta)\n\n\n\n\n\n","category":"type"},{"location":"API/#MeshArrays.gcmgrid","page":"API documentation","title":"MeshArrays.gcmgrid","text":"gcmgrid\n\ngcmgrid data structure. Available constructors:\n\ngcmgrid(path::String, class::String, nFaces::Int, \n        fSize::Array{NTuple{2, Int},1},\n        ioSize::Union{NTuple{2, Int},Array{Int64,2}},\n        ioPrec::Type, read::Function, write::Function)\n\nThe class can be set to \"LatLonCap\", \"CubeSphere\", \"PeriodicChannel\", \"PeriodicDomain\".\n\nFor example, A periodic channel (periodic in the x direction) of size 360 by 160, can be defined as follows.\n\npth=MeshArrays.GRID_LL360\nclass=\"PeriodicChannel\"\nioSize=(360, 160)\nioPrec=Float32\n\nγ=gcmgrid(pth, class, 1, [ioSize], ioSize, ioPrec, read, write)\n\nΓ=GridLoad(γ)    \n\nPlease refer to GridSpec and UnitGrid for more info related to gcmgrid options.\n\n\n\n\n\n","category":"type"},{"location":"API/#MeshArrays.varmeta","page":"API documentation","title":"MeshArrays.varmeta","text":"varmeta\n\nvarmeta data structure. By default, unit is missing (non-dimensional), position is fill(0.5,3) (cell center), time is missing, and name / long_name is unknown.\n\nAvailable constructors:\n\nvarmeta(unit::Union{Unitful.Units,Number,Missing},position::Array{Float64,1},\n        time::Union{DateTime,Missing,Array{DateTime,1}},name::String,long_name::String)\n\nAnd:\n\ndefaultmeta = varmeta(missing,fill(0.5,3),missing,\"unknown\",\"unknown\")\n\n\n\n\n\n","category":"type"},{"location":"API/","page":"API documentation","title":"API documentation","text":"More : ","category":"page"},{"location":"API/","page":"API documentation","title":"API documentation","text":"MeshArrays.gcmvector\nMeshArrays.gcmfaces","category":"page"},{"location":"API/#MeshArrays.gcmvector","page":"API documentation","title":"MeshArrays.gcmvector","text":"gcmvector{T, N}\n\ngcmvector data structure that can be used for   subsetting and indexing into a gcmarray.\n\ngcmvector{T,N}(grid::gcmgrid,f::Array{Array{T,1},N},\n         fSize::Array{NTuple{N, Int}},fIndex::Array{Int,1})\n\n\n\n\n\n\n","category":"type"},{"location":"API/#MeshArrays.gcmfaces","page":"API documentation","title":"MeshArrays.gcmfaces","text":"gcmfaces{T, N}\n\ngcmfaces data structure. Available constructors:\n\ngcmfaces{T,N}(grid::gcmgrid,f::Array{Array{T,N},1},\n         fSize::Array{NTuple{N, Int}}, aSize::NTuple{N,Int})\n\ngcmfaces(grid::gcmgrid,v1::Array{Array{T,N},1}) where {T,N}\ngcmfaces(grid::gcmgrid,::Type{T},\n         fSize::Array{NTuple{N, Int}}, aSize::NTuple{N,Int}) where {T,N}\n\ngcmfaces(grid::gcmgrid)\ngcmfaces(grid::gcmgrid,::Type{T})\ngcmfaces(grid::gcmgrid,::Type{T},n3::Int)\n\n\n\n\n\n","category":"type"},{"location":"API/#2.-Grids-And-I/O","page":"API documentation","title":"2. Grids And I/O","text":"","category":"section"},{"location":"API/","page":"API documentation","title":"API documentation","text":"UnitGrid\nsimple_periodic_domain\nGridSpec\nGridLoad\nGridLoadVar","category":"page"},{"location":"API/#MeshArrays.UnitGrid","page":"API documentation","title":"MeshArrays.UnitGrid","text":"UnitGrid(γ::gcmgrid;option=\"minimal\")\n\nGenerate a unit grid, where every grid spacing and area is 1, according to γ. \n\n\n\n\n\nUnitGrid(ioSize, tileSize; option=\"minimal\")\n\nGenerate a unit grid, where every grid spacing and area is 1, according to ioSize, tileSize. \n\nSince ioSize, tileSize defines a one to one mapping from global to tiled array, here we use  read_tiles, write_tiles instead of the default read, write. And we overwrite local XC,YC etc accordingly with global XC,YC etc . \n\n\n\n\n\n","category":"function"},{"location":"API/#MeshArrays.simple_periodic_domain","page":"API documentation","title":"MeshArrays.simple_periodic_domain","text":"simple_periodic_domain(np::Integer,nq=missing)\n\nSet up a simple periodic domain of size np x nq\n\nusing MeshArrays\nnp=16 #domain size is np x np\nΓ=simple_periodic_domain(np)\nisa(Γ.XC,MeshArray)\n\n\n\n\n\n","category":"function"},{"location":"API/#MeshArrays.GridSpec","page":"API documentation","title":"MeshArrays.GridSpec","text":"GridSpec(category=\"PeriodicDomain\",path=tempname(); ID=:unknown)\n\nSelect one of the pre-defined grids either by ID (keyword) or by category.\nReturn the corresponding gcmgrid specification, including the path where grid files can be accessed (path).\n\nselection by ID\n\n:LLC90\n:CS32\n:onedegree\n:default\n\nExample:\n\nusing MeshArrays\ng = GridSpec(ID=:LLC90)\n\nnote : the path to these fully supported grids are handled internally in MeshArrays.jl.\n\nby category and path\n\n\"PeriodicDomain\"\n\"PeriodicChannel\"\n\"CubeSphere\"\n\"LatLonCap\"`\n\nExamples:\n\nusing MeshArrays\ng = GridSpec()\ng = GridSpec(\"PeriodicChannel\",MeshArrays.GRID_LL360)\ng = GridSpec(\"CubeSphere\",MeshArrays.GRID_CS32)\ng = GridSpec(\"LatLonCap\",MeshArrays.GRID_LLC90)\nisa(g,gcmgrid)\n\n\n\n\n\n","category":"function"},{"location":"API/#MeshArrays.GridLoad","page":"API documentation","title":"MeshArrays.GridLoad","text":"GridLoad(γ=GridSpec(); ID=:default, option=:minimal)\n\nReturn a NamedTuple of grid variables read from files located in γ.path (see ?GridSpec).\noption : \n(default) option=:minimal means that only grid cell center positions (XC, YC) are loaded. \noption=:full provides a complete set of 2D grid variables. \noption=:full provides a complete set of 2D & 3d grid variables. \n\nBased on the MITgcm naming convention, grid variables are:\n\nXC, XG, YC, YG, AngleCS, AngleSN, hFacC, hFacS, hFacW, Depth.\nRAC, RAW, RAS, RAZ, DXC, DXG, DYC, DYG.\nDRC, DRF, RC, RF (one-dimensional)\n\nMITgcm documentation : \n\nhttps://mitgcm.readthedocs.io/en/latest/algorithm/algorithm.html#spatial-discretization-of-the-dynamical-equations\n\nusing MeshArrays\nγ = GridSpec(ID=:LLC90)\nΓ = GridLoad(γ;option=\"full\")\n\nisa(Γ.XC,MeshArray)\n\n\n\n\n\n","category":"function"},{"location":"API/#MeshArrays.GridLoadVar","page":"API documentation","title":"MeshArrays.GridLoadVar","text":"GridLoadVar(nam::String,γ::gcmgrid)\n\nReturn a grid variable read from files located in γ.path (see ?GridSpec, ?GridLoad).\n\nBased on the MITgcm naming convention, grid variables are:\n\nXC, XG, YC, YG, AngleCS, AngleSN, hFacC, hFacS, hFacW, Depth.\nRAC, RAW, RAS, RAZ, DXC, DXG, DYC, DYG.\nDRC, DRF, RC, RF (one-dimensional)\n\nMITgcm documentation : \n\nhttps://mitgcm.readthedocs.io/en/latest/algorithm/algorithm.html#spatial-discretization-of-the-dynamical-equations\n\nusing MeshArrays\n\nγ = GridSpec(\"CubeSphere\",MeshArrays.GRID_CS32)\nXC = GridLoadVar(\"XC\",γ)\n\nisa(XC,MeshArray)\n\n\n\n\n\n","category":"function"},{"location":"API/","page":"API documentation","title":"API documentation","text":"More : ","category":"page"},{"location":"API/","page":"API documentation","title":"API documentation","text":"Tiles\nTiles!\nexchange\nMeshArrays.read\nMeshArrays.read!\nMeshArrays.write","category":"page"},{"location":"API/#MeshArrays.Tiles","page":"API documentation","title":"MeshArrays.Tiles","text":"Tiles(γ::gcmgrid,ni::Int,nj::Int)\n\nDefine sudomain tiles of size ni,nj. Each tile is defined by a Dict where tile,face,i,j correspond to tile ID, face ID, index ranges.\n\nusing MeshArrays\nγ=GridSpec(\"LatLonCap\",MeshArrays.GRID_LLC90)\nτ=Tiles(γ,30,30)\n\nisa(τ[1],NamedTuple)\n\n\n\n\n\nTiles(τ::Array{Dict},x::MeshArray)\n\nReturn an Array of tiles which cover x according to tile partition τ.\n\nusing MeshArrays\nγ=GridSpec(\"LatLonCap\",MeshArrays.GRID_LLC90)\nd=γ.read(γ.path*\"Depth.data\",MeshArray(γ,γ.ioPrec))\nτ=Tiles(γ,30,30)\ntd=Tiles(τ,d)\n\nD=similar(d)\nTiles!(τ,td,D)\n\nisa(td[1],Array)\n\n\n\n\n\n","category":"function"},{"location":"API/#MeshArrays.Tiles!","page":"API documentation","title":"MeshArrays.Tiles!","text":"Tiles!(τ::Array,tx::Array,x::MeshArrays)\n\nMap tiles in tx according to tile partition τ into x.\n\n\n\n\n\n","category":"function"},{"location":"API/#MeshArrays.exchange","page":"API documentation","title":"MeshArrays.exchange","text":"exchange(fld::MeshArray)\n\nExchange / transfer data between neighboring arrays. Other methods are\n\nexchange(fld::MeshArray,N::Integer)\nexchange(u::MeshArray,v::MeshArray)\nexchange(u::MeshArray,v::MeshArray,N::Integer)\n\n\n\n\n\n","category":"function"},{"location":"API/#Base.read","page":"API documentation","title":"Base.read","text":"read(fil::String,x::MeshArray)\n\nRead array from file and return as a MeshArray. \n\nThe second argument (MeshArray or gcmgrid) provides the grid specifications (x.grid.ioSize). ```\n\n\n\n\n\nread(xx::Array,γ::gcmgrid)\n\nReformat Array data into a MeshArray shaped after γ.\n\n\n\n\n\nread(xx::Array,x::MeshArray)\n\nReformat Array data into a MeshArray similar to x. \n\n\n\n\n\n","category":"function"},{"location":"API/#Base.read!","page":"API documentation","title":"Base.read!","text":"read!(xx::Array,x::MeshArray)\n\nReformat array into MeshArray and write into x.\n\n\n\n\n\n","category":"function"},{"location":"API/#Base.write","page":"API documentation","title":"Base.write","text":"write(fil::String,x::MeshArray)\n\nWrite MeshArray to binary file. Other methods:\n\nwrite(xx::Array,x::MeshArray) #to Array\n\n\n\n\n\n","category":"function"},{"location":"API/#3.-Interpolation","page":"API documentation","title":"3. Interpolation","text":"","category":"section"},{"location":"API/","page":"API documentation","title":"API documentation","text":"interpolation_setup\nInterpolate\nInterpolationFactors\nStereographicProjection\nknn","category":"page"},{"location":"API/#MeshArrays.interpolation_setup","page":"API documentation","title":"MeshArrays.interpolation_setup","text":"interpolation_setup(fil::String)\n\nRead e.g. interp_coeffs_halfdeg.jld2\n\nfil=joinpath(tempdir(),\"interp_coeffs_halfdeg.jld2\")\nλ=MeshArrays.interpolation_setup(fil)\n\n\n\n\n\ninterpolation_setup(;Γ,lon,lat,path,url)\n\nDownload or recompute interpolation coefficients.\n\nλ=interpolation_setup() to download \"interpcoeffshalfdeg.jld2\" \nλ=interpolation_setup(Γ=Γ) to recompute interpolation to lon,lat\n\n\n\n\n\n","category":"function"},{"location":"API/#MeshArrays.Interpolate","page":"API documentation","title":"MeshArrays.Interpolate","text":"Interpolate(z_in::MeshArray,f,i,j,w)\n\nusing MeshArrays\n\nγ=GridSpec(\"LatLonCap\",MeshArrays.GRID_LLC90)\nΓ=GridLoad(γ; option=\"full\")\n\nlon=[i for i=-179.5:1.0:179.5, j=-89.5:1.0:89.5]\nlat=[j for i=-179.5:1.0:179.5, j=-89.5:1.0:89.5]\n(f,i,j,w,j_f,j_x,j_y)=InterpolationFactors(Γ,vec(lon),vec(lat))\nDD=Interpolate(Γ.Depth,f,i,j,w)\n\nusing CairoMakie\nheatmap(vec(lon[:,1]),vec(lat[1,:]),DD,colorrange=(0.,6000.))\n\n\n\n\n\nInterpolate(z_in::MeshArray,f,i,j,w)\n\nusing MeshArrays\n\nγ=GridSpec(\"LatLonCap\",MeshArrays.GRID_LLC90)\nΓ=GridLoad(γ; option=\"full\")\n\nlon=[i for i=-179.5:1.0:179.5, j=-89.5:1.0:89.5]\nlat=[j for i=-179.5:1.0:179.5, j=-89.5:1.0:89.5]\n(f,i,j,w,j_f,j_x,j_y)=InterpolationFactors(Γ,vec(lon),vec(lat))\nL=(lon=lon, lat=lat, f=f, i=i, j=j, w=w)\n\nusing CairoMakie\nheatmap(Interpolate(Γ.Depth,L)...,colorrange=(0.,6000.))\n\nor \n\nheatmap(Γ.Depth,interpolation=L,colorrange=(0.,6000.))\n\n\n\n\n\n","category":"function"},{"location":"API/#MeshArrays.InterpolationFactors","page":"API documentation","title":"MeshArrays.InterpolationFactors","text":"InterpolationFactors(Γ,lon::Array{T,1},lat::Array{T,1})\n\nCompute interpolation coefficients etc from grid Γ to lon,lat\n\nusing MeshArrays\nγ=GridSpec(\"CubeSphere\",MeshArrays.GRID_CS32)\nΓ=GridLoad(γ; option=\"full\")\nlon=collect(45.:0.1:46.); lat=collect(60.:0.1:61.)\n(f,i,j,w,j_f,j_x,j_y)=InterpolationFactors(Γ,lon,lat)\nYC=Interpolate(Γ.YC,f,i,j,w)\nextrema(i)==(9,10)\n\n\n\n\n\n","category":"function"},{"location":"API/#MeshArrays.StereographicProjection","page":"API documentation","title":"MeshArrays.StereographicProjection","text":"StereographicProjection(XC0,YC0,XC,YC)\n\nApply stereographic projection that puts XC0,YC0 at 0.0,0.0 to target point(s) XC,YC\n\nlon=collect(45.:0.1:46.); lat=collect(60.:0.1:61.)\nx,y=StereographicProjection(45.,60.,lon,lat)\n\n\n\n\n\n","category":"function"},{"location":"API/#NearestNeighbors.knn","page":"API documentation","title":"NearestNeighbors.knn","text":"knn(xgrid,ygrid::MeshArray,x,y::Array{T,1},k::Int)\n\nFind k nearest neighbors to each point in x,y on xgrid,ygrid\n\nlon=collect(0.1:0.5:2.1); lat=collect(0.1:0.5:2.1);\n(f,i,j,c)=knn(Γ.XC,Γ.YC,lon,lat)\n\n\n\n\n\n","category":"function"},{"location":"API/#4.-Vector-Fields","page":"API documentation","title":"4. Vector Fields","text":"","category":"section"},{"location":"API/","page":"API documentation","title":"API documentation","text":"curl\nconvergence\ngradient\nScalarPotential\nVectorPotential\nThroughFlow\nUVtoTransport\nUVtoUEVN","category":"page"},{"location":"API/#MeshArrays.curl","page":"API documentation","title":"MeshArrays.curl","text":"curl(u::MeshArray,v::MeshArray,Γ::NamedTuple)\n\nCompute curl of a velocity field.\n\n\n\n\n\n","category":"function"},{"location":"API/#MeshArrays.convergence","page":"API documentation","title":"MeshArrays.convergence","text":"convergence(uFLD::MeshArray,vFLD::MeshArray)\n\nCompute convergence of a vector field\n\n\n\n\n\n","category":"function"},{"location":"API/#MeshArrays.gradient","page":"API documentation","title":"MeshArrays.gradient","text":"gradient(inFLD::MeshArray,Γ::NamedTuple)\n\nCompute spatial derivatives. Other methods:\n\ngradient(inFLD::MeshArray,Γ::NamedTuple,doDIV::Bool)\ngradient(inFLD::MeshArray,iDXC::MeshArray,iDYC::MeshArray)\n\n\n\n\n\n","category":"function"},{"location":"API/#MeshArrays.ScalarPotential","page":"API documentation","title":"MeshArrays.ScalarPotential","text":"ScalarPotential(TrspCon)\n\nScalar potential inversion.\n\nTrspPot=ScalarPotential(TrspCon)\n\n\n\n\n\n","category":"function"},{"location":"API/#MeshArrays.VectorPotential","page":"API documentation","title":"MeshArrays.VectorPotential","text":"VectorPotential(TrspX,TrspY,Γ,method::Int=1)\n\nVector potential inversion.\n\nTrspPot=ScalarPotential(TrspCon)\n\n\n\n\n\n","category":"function"},{"location":"API/#MeshArrays.ThroughFlow","page":"API documentation","title":"MeshArrays.ThroughFlow","text":"ThroughFlow(VectorField,IntegralPath,Γ::NamedTuple)\n\nCompute transport through an integration path\n\n\n\n\n\n","category":"function"},{"location":"API/#MeshArrays.UVtoTransport","page":"API documentation","title":"MeshArrays.UVtoTransport","text":"UVtoTransport(U,V,G::NamedTuple)\n\nConvert e.g. velocity (m/s) to transport (m^3/s) by multiplying by DRF and DXG,DYG.\n\n\n\n\n\n","category":"function"},{"location":"API/#MeshArrays.UVtoUEVN","page":"API documentation","title":"MeshArrays.UVtoUEVN","text":"UVtoUEVN(u,v,G::NamedTuple)\n\nInterpolate to grid cell centers (uC,vC)\nConvert to Eastward/Northward components (uE,vN)\n\nNote: land masking u,v with NaNs preemptively can be adequate.\n\n\n\n\n\n","category":"function"},{"location":"API/#5.-Other","page":"API documentation","title":"5. Other","text":"","category":"section"},{"location":"API/","page":"API documentation","title":"API documentation","text":"LatitudeCircles\nTransect\nisosurface\nMA_datadep","category":"page"},{"location":"API/#MeshArrays.LatitudeCircles","page":"API documentation","title":"MeshArrays.LatitudeCircles","text":"LatitudeCircles(LatValues,Γ::NamedTuple; format=:gridpath, range=(0.0,360.0))\n\nCompute integration paths that follow latitude circles, within the specified longitude range.\n\n\n\n\n\n","category":"function"},{"location":"API/#MeshArrays.Transect","page":"API documentation","title":"MeshArrays.Transect","text":"Transect(name,lons,lats,Γ; segment=:short, format=:gridpath)\n\nCompute integration paths that follow a great circle between two geolocations give by lons, lats.\n\n\n\n\n\n","category":"function"},{"location":"API/#MeshArrays.isosurface","page":"API documentation","title":"MeshArrays.isosurface","text":"isosurface(θ,T,Γ)\n\nDepth of isosurface θ=T\n\n\n\n\n\n","category":"function"},{"location":"API/#MeshArrays.MA_datadep","page":"API documentation","title":"MeshArrays.MA_datadep","text":"MA_datadep(nam=\"countries_shp1\")\n\nDownload data dependency with predefined name; currently :\n\n\"countries_shp1\"\n\"countries_geojson1\"\n\"basemap_jpg1\"\n\"interp_halfdeg\"\n\n\n\n\n\n","category":"function"}]
}
