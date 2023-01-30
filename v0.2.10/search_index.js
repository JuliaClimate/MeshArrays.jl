var documenterSearchIndex = {"docs":
[{"location":"detail/#Additional-Detail-1","page":"Detail","title":"Additional Detail","text":"","category":"section"},{"location":"detail/#","page":"Detail","title":"Detail","text":"Functions like GridSpec(\"LLC90\") return a gcmgrid struct that contains the basic specification of a global grid. This is not the grid itself – just a few parameters, ranges, and possibly a path to grid files. A gcmgrid is embeded in each MeshArray instance for which it provides a blueprint. It specifies how an array collection forms a global mesh and allows e.g. the exchange function to dispatch to the appropriate method. ","category":"page"},{"location":"detail/#","page":"Detail","title":"Detail","text":"Various configurations that are commonly used in Earth System Models are readily implemented using the concrete type called MeshArray. This type is in fact an alias for more specific types that can be used interchangeably via MeshArray (initially: gcmfaces or gcmarray).","category":"page"},{"location":"detail/#","page":"Detail","title":"Detail","text":"Within a MeshArray, a whole Earth System Model grid is represented as an array of elementary arrays. Each one of these represents a subdomain. For example, a gcmarray instance for one Earth map x has a column array x.f of elementary 2D arrays of various sizes. demo1 illustrates how one easily operates MeshArray structs via standard and specialized functions. In brief, a MeshArray should be used just like a common Array.","category":"page"},{"location":"detail/#","page":"Detail","title":"Detail","text":"Background: MeshArrays.jl is rooted in a Matlab / Octave package called gcmfaces, which was introduced in Forget et al., 2015 (doi:10.5194/gmd-8-3071-2015). GCM is an acronym for General Circulation Model, or Global Climate Model, and faces can mean meshes, arrays, facets, or subdomains (these are the elements of x.f in a MeshArray instance x).","category":"page"},{"location":"main/#Main-Features-1","page":"Main","title":"Main Features","text":"","category":"section"},{"location":"main/#","page":"Main","title":"Main","text":"The elements of a MeshArray are arrays. These elementary arrays typically represent subdomains inter-connected at their edges. The organization and connections between subdomains is determined by a user-specified gcmgrid which is embeded inside each MeshArray instance. ","category":"page"},{"location":"main/#","page":"Main","title":"Main","text":"Interpolate can be used to interpolate a MeshArray to any location (i.e. arbitrary longitude, latitude pair). Exchange methods transfer data between neighboring arrays to extend computational subdomains – this is often needed in analyses of climate or ocean model output. ","category":"page"},{"location":"main/#","page":"Main","title":"Main","text":"The current default for MeshArray is the gcmarray type and an instance H is shown below. This example is based on a grid known as LatLonCap where each global map is associated with 5 subdomains. Hence, H.f is a (5, 50) array when H represents a gridded variable on 50 depth levels, and elements of  H.f are arrays of size (90, 270), (90, 90), or (270, 90). ","category":"page"},{"location":"main/#","page":"Main","title":"Main","text":"julia> show(D)\n  name        = Depth\n  unit        = m\n  data type   = Float64\n  cell pos.   = [0.5, 0.5]\n\n  tile array  = (5,)\n  tile sizes  = (90, 270)\n                (90, 270)\n                (90, 90)\n                (270, 90)\n                (270, 90)\n\n  grid class  = LatLonCap\n  MeshArray   = gcmarray \n  version     = 0.2.7 ","category":"page"},{"location":"main/#","page":"Main","title":"Main","text":"The underlying, MeshArray, data structure is:","category":"page"},{"location":"main/#","page":"Main","title":"Main","text":"struct gcmarray{T, N} <: AbstractMeshArray{T, N}\n   grid::gcmgrid\n   meta::varmeta\n   f::Array{Array{T,2},N}\n   fSize::Array{NTuple{2, Int}}\n   fIndex::Array{Int,1}\n   version::String\nend","category":"page"},{"location":"main/#","page":"Main","title":"Main","text":"A MeshArray generally behaves just like an Array including for operations listed below. The broadcasting function has been customized so that it reaches elements of each elementary array.","category":"page"},{"location":"main/#","page":"Main","title":"Main","text":"size(D)\neltype(D)\nview(D,:)\n\nD .* 1.0\nD .* D\n1000*D\nD*1000\n\nD[findall(D .> 300.)] .= NaN\nD[findall(D .< 1.)] .= NaN\n\nD[1]=0.0 .+ D[1]\ntmp=cos.(D)","category":"page"},{"location":"main/#","page":"Main","title":"Main","text":"In addition, Mesharray specific functions like exchange can alter the internal structure of a MeshArray. Elementary array sizes are thus larger in show(exchange(D)) than show(D).","category":"page"},{"location":"main/#","page":"Main","title":"Main","text":"julia> show(exchange(D))\n  tile sizes  = (92, 272)\n                (92, 272)\n                (92, 92)\n                (272, 92)\n                (272, 92)","category":"page"},{"location":"main/#Embedded-Meta-Data-1","page":"Main","title":"Embedded Meta Data","text":"","category":"section"},{"location":"main/#","page":"Main","title":"Main","text":"A MeshArray includes a gcmgrid specification which can be constructed as outlined below.","category":"page"},{"location":"main/#","page":"Main","title":"Main","text":"gcmgrid(path::String, class::String, \n        nFaces::Int, fSize::Array{NTuple{2, Int},1}, \n        ioSize::Array{Int64,2}, ioPrec::Type, \n        read::Function, write::Function)","category":"page"},{"location":"main/#","page":"Main","title":"Main","text":"Importantly, a gcmgrid does not contain any actual grid data – hence its memory footprint is minimal. Grid variables are instead read to memory only when needed e.g. as shown below. To make this easy, each gcmgrid includes a pair of read / write methods to allow for basic I/O at any time. These methods are typically specified by the user although defaults are provided. ","category":"page"},{"location":"main/#","page":"Main","title":"Main","text":"using MeshArrays, Unitful\nγ=GridSpec(\"LatLonCap\",\"GRID_LLC90/\")\nm=MeshArrays.varmeta(u\"m\",fill(0.5,2),\"Depth\",\"Depth\")\nD=γ.read(γ.path*\"Depth.data\",MeshArray(γ,Float64;meta=m))","category":"page"},{"location":"main/#","page":"Main","title":"Main","text":"The above commands define a MeshArray called D which is the one displayed at the top of this section. A definition of the varmeta structure is reported below. The position of a D point within its grid cell is given as x ∈ [0. 1.] in each direction.","category":"page"},{"location":"main/#","page":"Main","title":"Main","text":"varmeta(unit::Union{Unitful.AbstractQuantity,Number},\n        position::Array{Float64,1},\n        name::String,long_name::String)","category":"page"},{"location":"main/#Examples,-Interpolation,-and-Plots-1","page":"Main","title":"Examples, Interpolation, & Plots","text":"","category":"section"},{"location":"main/#","page":"Main","title":"Main","text":"The JuliaCon-2018 presentation relied on two Jupyter notebooks available in GlobalOceanNotebooks/DataStructures which are also included here in examples/Demos.jl. GlobalOceanNotebooks/OceanTransports provides use case examples related to Earth System transports.","category":"page"},{"location":"main/#","page":"Main","title":"Main","text":"A simple way to plot a MeshArray consists in plotting each elementary array separately. Other methods that e.g. produce global maps and projections are illustrated in the notebooks. A simple one is shown below that demonstrates the included interpolation scheme.","category":"page"},{"location":"main/#","page":"Main","title":"Main","text":"p=dirname(pathof(MeshArrays));\nusing Plots; include(joinpath(p,\"../examples/Plots.jl\"));\nheatmap(D,title=\"Ocean Depth\",clims=(0.,6000.))","category":"page"},{"location":"main/#","page":"Main","title":"Main","text":"(Image: OceanDepthMap)","category":"page"},{"location":"main/#","page":"Main","title":"Main","text":"lon=[i for i=-179.5:1.0:179.5, j=-89.5:1.0:89.5]\nlat=[j for i=-179.5:1.0:179.5, j=-89.5:1.0:89.5]\n\nΓ=GridLoad(GridSpec(\"LatLonCap\",\"GRID_LLC90/\"))\n(f,i,j,w)=InterpolationFactors(Γ,vec(lon),vec(lat))\nDD=Interpolate(Γ[\"Depth\"],f,i,j,w)\n\ncontourf(vec(lon[:,1]),vec(lat[1,:]),DD,clims=(0.,6000.))","category":"page"},{"location":"main/#","page":"Main","title":"Main","text":"(Image: OceanDepthMap)","category":"page"},{"location":"#MeshArrays.jl-documentation-1","page":"Home","title":"MeshArrays.jl documentation","text":"","category":"section"},{"location":"#","page":"Home","title":"Home","text":"MeshArrays.jl defines an array type that can contain / organize / distribute collections of inter-connected arrays as done in climate models (see Earth Model Grids below). The MeshArray type is a sub-type of AbstractArray with an outer array where each element is itself a 2D inner array. This setup potentially allows different choices for the outer and inner arrays – for example DistributedArrays and AxisArrays, respectively, could be an option.","category":"page"},{"location":"#","page":"Home","title":"Home","text":"MeshArrays.jl thus provides a simple but general solution to analyze or e.g. simulate climate system variables. The internals of a MeshArray are regulated by a few index ranges, array size specifications, and inter-connection rules that are encoded in a gcmgrid. A second, also lightweight, structure called varmeta contains the MeshArray variable name, unit, and position on the grid. A general approach like this is useful because climate models often involve advanced domain decompositions (see Earth Model Grids), and many variables, which puts a burden on users. ","category":"page"},{"location":"#","page":"Home","title":"Home","text":"Encoding the grid specification inside the MeshArray data type allows user to manipulate MeshArrays just like they would manipulate Arrays without having to keep track of model grid details. In addition, the provided exchange methods readily transfer data between connected subdomains to extend them at the sides. This makes it easy to compute e.g. partial derivatives and related operators like gradients, curl, or divergences over subdomain edges as often needed for precise computation of transports, budgets, etc using climate model output.","category":"page"},{"location":"#","page":"Home","title":"Home","text":"MeshArrays.jl was first introduced as as gcmfaces.jl in a JuliaCon-2018 presentation. This notebook folder demonstrates how its data structures can be used to accurately analyze the General Ocean Circulation. Examples include computations of ocean heat transport and streamfunctions that are important and widely studied aspects of the climate system.","category":"page"},{"location":"#","page":"Home","title":"Home","text":"Contents:","category":"page"},{"location":"#","page":"Home","title":"Home","text":"Pages = [\"index.md\",\"main.md\",\"detail.md\",\"API.md\"]\nDepth = 3","category":"page"},{"location":"#Install-and-Test-1","page":"Home","title":"Install & Test","text":"","category":"section"},{"location":"#","page":"Home","title":"Home","text":"using Pkg\nPkg.add(\"MeshArrays\")\nPkg.test(\"MeshArrays\")","category":"page"},{"location":"#","page":"Home","title":"Home","text":"Julia's package manager, Pkg.jl, is documented the main Julia doc and here in details.","category":"page"},{"location":"#Basic-Examples-1","page":"Home","title":"Basic Examples","text":"","category":"section"},{"location":"#","page":"Home","title":"Home","text":"Examples below (1) generate a grid configuration, (2) seed a 2D field of random noise, (3) smooth out this field, and (4) plot subdomain arrays. Smoothing is done via a lateral diffusion equation through time to illustrate how MeshArray computes partial derivatives & transfers data between neighboring subdomains. Examples 2 & 3 illustrate grid configurations commonly used in global models.","category":"page"},{"location":"#","page":"Home","title":"Home","text":"[A] 16 subdomains, with 40x40 grid points each, covering a doubly periodic domain","category":"page"},{"location":"#","page":"Home","title":"Home","text":"using MeshArrays; p=dirname(pathof(MeshArrays))\nγ,Γ=GridOfOnes(\"PeriodicDomain\",16,20)\n\ninclude(joinpath(p,\"../examples/Demos.jl\"))\n(xi,xo,_,_)=demo2(Γ);\nshow(xo)\n\nusing Plots; plotlyjs()\ninclude(joinpath(p,\"../examples/Plots.jl\"))\nheatmap(xo,clims=(-0.25,0.25))","category":"page"},{"location":"#","page":"Home","title":"Home","text":"Grid scale noise Smoothed noise\n(Image: raw) (Image: smooth)","category":"page"},{"location":"#","page":"Home","title":"Home","text":"[B] 6 subdomains, with 100x100 points each, covering the six faces of a cube","category":"page"},{"location":"#","page":"Home","title":"Home","text":"γ,Γ=GridOfOnes(\"CubeSphere\",6,100)\nD=demo2(Γ)","category":"page"},{"location":"#","page":"Home","title":"Home","text":"[C] Global Model Grid with 5 uneven subdomains, variable spacing, & continents","category":"page"},{"location":"#","page":"Home","title":"Home","text":"This requires downloading a pre-defined global ocean grid from the MITgcm community.","category":"page"},{"location":"#","page":"Home","title":"Home","text":"#run(`git clone https://github.com/gaelforget/GRID_LLC90`)\nΓ=GridLoad(GridSpec(\"LatLonCap\",\"GRID_LLC90/\"))\nD=demo2(Γ)\nheatmap(D[2],clims=(-0.25,0.25))","category":"page"},{"location":"#Earth-Model-Grids-1","page":"Home","title":"Earth Model Grids","text":"","category":"section"},{"location":"#","page":"Home","title":"Home","text":"-(Image: EarthGrids)","category":"page"},{"location":"API/#API-Guide-1","page":"API","title":"API Guide","text":"","category":"section"},{"location":"API/#","page":"API","title":"API","text":"","category":"page"},{"location":"API/#","page":"API","title":"API","text":"Modules = [MeshArrays]\nOrder   = [:type,:function]","category":"page"},{"location":"API/#MeshArrays.AbstractMeshArray","page":"API","title":"MeshArrays.AbstractMeshArray","text":"AbstractMeshArray{T, N}\n\nSubtype of AbstractArray{T, N}\n\n\n\n\n\n","category":"type"},{"location":"API/#MeshArrays.gcmgrid","page":"API","title":"MeshArrays.gcmgrid","text":"gcmgrid\n\ngcmgrid data structure. Available constructors:\n\ngcmgrid(path::String, class::String,\n        nFaces::Int, fSize::Array{NTuple{2, Int},1},\n        ioSize::Array{Int64,2}, ioPrec::Type,\n        read::Function, write::Function)\n\n\n\n\n\n","category":"type"},{"location":"API/#MeshArrays.GridAddWS!-Tuple{Dict}","page":"API","title":"MeshArrays.GridAddWS!","text":"GridAddWS!(Γ::Dict)\n\nCompute XW, YW, XS, and YS (vector field locations) from XC, YC (tracer field locations) and add them to Γ.\n\nΓ=GridLoad(GridSpec(\"LatLonCap\",\"GRID_LLC90/\"))\nGridAddWS!(Γ)\n\n\n\n\n\n","category":"method"},{"location":"API/#MeshArrays.GridLoad-Tuple{gcmgrid}","page":"API","title":"MeshArrays.GridLoad","text":"GridLoad(γ::gcmgrid)\n\nReturn a Dict of grid variables read from files located in γ.path (see ?GridSpec).\n\nBased on the MITgcm naming convention, grid variables are:\n\nXC, XG, YC, YG, AngleCS, AngleSN, hFacC, hFacS, hFacW, Depth.\nRAC, RAW, RAS, RAZ, DXC, DXG, DYC, DYG.\nDRC, DRF, RC, RF (one-dimensional)\n\n\n\n\n\n","category":"method"},{"location":"API/#MeshArrays.GridOfOnes-Tuple{Any,Any,Any}","page":"API","title":"MeshArrays.GridOfOnes","text":"GridOfOnes(grTp,nF,nP)\n\nDefine all-ones grid variables instead of using GridSpec & GridLoad. E.g.\n\nγ,Γ=GridOfOnes(\"CubeSphere\",6,20);\n\n\n\n\n\n","category":"method"},{"location":"API/#MeshArrays.GridSpec","page":"API","title":"MeshArrays.GridSpec","text":"GridSpec(GridName,GridParentDir=\"./\")\n\nReturn a gmcgrid specification that provides grid files path, class, nFaces, ioSize, facesSize, ioPrec, & a read function (not yet) using hard-coded values for \"PeriodicDomain\", \"PeriodicChannel\", \"CubeSphere\", and `\"LatLonCap\" for now.\n\n\n\n\n\n","category":"function"},{"location":"API/#MeshArrays.Interpolate-Tuple{MeshArrays.gcmarray,Any,Any,Any,Any}","page":"API","title":"MeshArrays.Interpolate","text":"Interpolate(z_in::MeshArray,f,i,j,w)\n\nlon=[i for i=-179.5:1.0:179.5, j=-89.5:1.0:89.5]\nlat=[j for i=-179.5:1.0:179.5, j=-89.5:1.0:89.5]\n\nΓ=GridLoad(GridSpec(\"LatLonCap\",\"GRID_LLC90/\"))\n(f,i,j,w,j_f,j_x,j_y)=InterpolationFactors(Γ,vec(lon),vec(lat))\nDD=Interpolate(Γ[\"Depth\"],f,i,j,w)\n\nusing Plots\ncontourf(vec(lon[:,1]),vec(lat[1,:]),DD,clims=(0.,6000.))\n\n\n\n\n\n","category":"method"},{"location":"API/#MeshArrays.InterpolationFactors-Union{Tuple{T}, Tuple{Any,Array{T,1},Array{T,1}}} where T","page":"API","title":"MeshArrays.InterpolationFactors","text":"InterpolationFactors(Γ,lon::Array{T,1},lat::Array{T,1})\n\nCompute interpolation coefficients etc from grid Γ to lon,lat\n\nlon=collect(45.:0.1:46.); lat=collect(60.:0.1:61.)\n(f,i,j,w,j_f,j_x,j_y)=InterpolationFactors(Γ,lon,lat)\n\n\n\n\n\n","category":"method"},{"location":"API/#MeshArrays.LatitudeCircles-Tuple{Any,Dict}","page":"API","title":"MeshArrays.LatitudeCircles","text":"LatitudeCircles(LatValues,Γ::Dict)\n\nCompute integration paths that follow latitude circles\n\n\n\n\n\n","category":"method"},{"location":"API/#MeshArrays.ScalarPotential-Tuple{Any}","page":"API","title":"MeshArrays.ScalarPotential","text":"ScalarPotential(TrspCon)\n\nScalar potential inversion.\n\nTrspPot=ScalarPotential(TrspCon)\n\n\n\n\n\n","category":"method"},{"location":"API/#MeshArrays.StereographicProjection-NTuple{4,Any}","page":"API","title":"MeshArrays.StereographicProjection","text":"StereographicProjection(XC0,YC0,XC,YC)\n\nApply stereographic projection that puts XC0,YC0 at 0.0,0.0 to target point(s) XC,YC\n\nlon=collect(45.:0.1:46.); lat=collect(60.:0.1:61.)\nx,y=StereographicProjection(45.,60.,lon,lat)\n\n\n\n\n\n","category":"method"},{"location":"API/#MeshArrays.ThroughFlow-Tuple{Any,Any,Dict}","page":"API","title":"MeshArrays.ThroughFlow","text":"ThroughFlow(VectorField,IntegralPath,Γ::Dict)\n\nCompute transport through an integration path\n\n\n\n\n\n","category":"method"},{"location":"API/#MeshArrays.Tiles-Tuple{Array{Dict,N} where N,MeshArrays.gcmarray}","page":"API","title":"MeshArrays.Tiles","text":"Tiles(τ::Array{Dict},x::MeshArray)\n\nReturn an Array of tiles which cover x according to tile partition τ.\n\nγ=GridSpec(\"LatLonCap\",\"GRID_LLC90/\")\nd=γ.read(γ.path*\"Depth.data\",MeshArray(γ,γ.ioPrec))\nτ=Tiles(γ,30,30)\ntd=Tiles(τ,d)\n\n\n\n\n\n","category":"method"},{"location":"API/#MeshArrays.Tiles-Tuple{gcmgrid,Int64,Int64}","page":"API","title":"MeshArrays.Tiles","text":"Tiles(γ::gcmgrid,ni::Int,nj::Int)\n\nDefine sudomain tiles of size ni,nj. Each tile is defined by a Dict where tile,face,i,j correspond to tile ID, face ID, index ranges.\n\nγ=GridSpec(\"LatLonCap\",\"GRID_LLC90/\")\nτ=Tiles(γ,30,30)\n\n\n\n\n\n","category":"method"},{"location":"API/#MeshArrays.VectorPotential","page":"API","title":"MeshArrays.VectorPotential","text":"VectorPotential(TrspX,TrspY,Γ,method::Int=1)\n\nVector potential inversion.\n\nTrspPot=ScalarPotential(TrspCon)\n\n\n\n\n\n","category":"function"},{"location":"API/#MeshArrays.convergence-Tuple{MeshArrays.gcmarray,MeshArrays.gcmarray}","page":"API","title":"MeshArrays.convergence","text":"convergence(uFLD::MeshArray,vFLD::MeshArray)\n\nCompute convergence of a vector field\n\n\n\n\n\n","category":"method"},{"location":"API/#MeshArrays.exchange-Tuple{MeshArrays.gcmarray}","page":"API","title":"MeshArrays.exchange","text":"exchange(fld::MeshArray)\n\nExchange / transfer data between neighboring arrays. Other methods are\n\nexchange(fld::MeshArray,N::Integer)\nexchange(u::MeshArray,v::MeshArray)\nexchange(u::MeshArray,v::MeshArray,N::Integer)\n\n\n\n\n\n","category":"method"},{"location":"API/#MeshArrays.gradient-Tuple{MeshArrays.gcmarray,Dict}","page":"API","title":"MeshArrays.gradient","text":"gradient(inFLD::MeshArray,Γ::Dict)\n\nCompute spatial derivatives. Other methods:\n\ngradient(inFLD::MeshArray,Γ::Dict,doDIV::Bool)\ngradient(inFLD::MeshArray,iDXC::MeshArray,iDYC::MeshArray)\n\n\n\n\n\n","category":"method"},{"location":"API/#MeshArrays.mask-Tuple{MeshArrays.gcmarray,Number}","page":"API","title":"MeshArrays.mask","text":"mask(fld::MeshArray, val::Number)\n\nReplace non finite values with val. Other methods:\n\nmask(fld::MeshArray)\nmask(fld::MeshArray, val::Number, noval::Number)\n\n\n\n\n\n","category":"method"},{"location":"API/#MeshArrays.smooth-Tuple{MeshArrays.gcmarray,MeshArrays.gcmarray,MeshArrays.gcmarray,Dict}","page":"API","title":"MeshArrays.smooth","text":"smooth(FLD::MeshArray,DXCsm::MeshArray,DYCsm::MeshArray,Γ::Dict)\n\nSmooth out scales below DXCsm / DYCsm via diffusion\n\n\n\n\n\n","category":"method"},{"location":"API/#NearestNeighbors.knn-Union{Tuple{T}, Tuple{MeshArrays.gcmarray,MeshArrays.gcmarray,Array{T,1},Array{T,1}}, Tuple{MeshArrays.gcmarray,MeshArrays.gcmarray,Array{T,1},Array{T,1},Any}} where T","page":"API","title":"NearestNeighbors.knn","text":"knn(xgrid,ygrid::MeshArray,x,y::Array{T,1},k::Int)\n\nFind k nearest neighbors to each point in x,y on xgrid,ygrid\n\nlon=collect(0.1:0.5:2.1); lat=collect(0.1:0.5:2.1);\n(f,i,j,c)=knn(Γ[\"XC\"],Γ[\"YC\"],lon,lat)\n\n\n\n\n\n","category":"method"},{"location":"API/#MeshArrays.gcmarray","page":"API","title":"MeshArrays.gcmarray","text":"gcmarray{T, N, AT}\n\ngcmarray data structure. Available constructors:\n\ngcmarray{T,N,AT}(grid::gcmgrid,meta::varmeta,f::Array{AT,N},\n         fSize::Array{NTuple{N, Int}},fIndex::Array{Int,1},v::String)\n\ngcmarray(grid::gcmgrid,f::Array{Array{T,2},N}) where {T,N}\ngcmarray(grid::gcmgrid,f::Array{Array{T,N},1}) where {T,N}\n\ngcmarray(grid::gcmgrid,fSize::Array{NTuple{N, Int}},fIndex::Array{Int,1})\ngcmarray(<same as above>,n3::Int)\ngcmarray(<same as above>,n3::Int,n4::Int)\n\ngcmarray(grid::gcmgrid)\ngcmarray(grid::gcmgrid,::Type{T})\ngcmarray(grid::gcmgrid,::Type{T},n3::Int)\ngcmarray(grid::gcmgrid,::Type{T},n3::Int,n4::Int)\n\n\n\n\n\n","category":"type"},{"location":"API/#MeshArrays.gcmfaces","page":"API","title":"MeshArrays.gcmfaces","text":"gcmfaces{T, N}\n\ngcmfaces data structure. Available constructors:\n\ngcmfaces{T,N}(grid::gcmgrid,f::Array{Array{T,N},1},\n         fSize::Array{NTuple{N, Int}}, aSize::NTuple{N,Int})\n\ngcmfaces(grid::gcmgrid,v1::Array{Array{T,N},1}) where {T,N}\ngcmfaces(grid::gcmgrid,::Type{T},\n         fSize::Array{NTuple{N, Int}}, aSize::NTuple{N,Int}) where {T,N}\n\ngcmfaces(grid::gcmgrid)\ngcmfaces(grid::gcmgrid,::Type{T})\ngcmfaces(grid::gcmgrid,::Type{T},n3::Int)\n\n\n\n\n\n","category":"type"},{"location":"API/#MeshArrays.gcmsubset","page":"API","title":"MeshArrays.gcmsubset","text":"gcmsubset{T, N}\n\ngcmsubset data structure for subsets of gcmfaces. Available constructors:\n\ngcmsubset{T,N}(grid::gcmgrid,f::Array{Array{T,N},1},\n               fSize::Array{NTuple{N, Int}},aSize::NTuple{N, Int},\n               i::Array{Array{T,N},1},iSize::Array{NTuple{N, Int}})\ngcmsubset(grid::gcmgrid,::Type{T},fSize::Array{NTuple{N, Int}},\n          aSize::NTuple{N,Int},dims::NTuple{N,Int}) where {T,N}\n\n\n\n\n\n","category":"type"},{"location":"API/#MeshArrays.gcmvector","page":"API","title":"MeshArrays.gcmvector","text":"gcmvector{T, N}\n\ngcmvector data structure that can be used for   subsetting and indexing into a gcmarray.\n\ngcmvector{T,N}(grid::gcmgrid,f::Array{Array{T,1},N},\n         fSize::Array{NTuple{N, Int}},fIndex::Array{Int,1})\n\n\n\n\n\n\n","category":"type"},{"location":"API/#MeshArrays.varmeta","page":"API","title":"MeshArrays.varmeta","text":"varmeta\n\nvarmeta data structure. By default, unit is 1.0 (non-dimensional), position is fill(0.5,3) (cell center), and name / long_name is unknown.\n\nAvailable constructors:\n\nvarmeta(unit::Union{Unitful.AbstractQuantity,Number},position::Array{Float64,1},\n        name::String,long_name::String)\n\nAnd:\n\ndefaultmeta = varmeta(1.0,fill(0.5,3),\"unknown\",\"unknown\")\n\n\n\n\n\n","category":"type"},{"location":"API/#Base.findall-Union{Tuple{MeshArrays.gcmarray{Bool,N,AT} where AT}, Tuple{N}} where N","page":"API","title":"Base.findall","text":"findall(A::gcmarray{Bool})\n\nReturn a gcmvector of the true indices in A. This allows:\n\nfindall(A.<0) #gcmvector of CartesianIndex{2}\nA[findall(A.<0)] #gcmvector of eltype(A)\nview(A,findall(A.<0)) #CatView of eltype(A)\n\nA[findall(A.<0)]=B[findall(A.<0)]\nA[findall(A.<0)].=view(B,findall(A.<0))\nA[findall(A.<0)].=NaN\n\n\n\n\n\n","category":"method"},{"location":"API/#Base.read-Tuple{String,MeshArrays.gcmarray}","page":"API","title":"Base.read","text":"read(fil::String,x::MeshArray)\n\nRead binary file to MeshArray. Other methods:\n\nread(xx::Array,x::MeshArray) #from Array\n\n\n\n\n\n","category":"method"},{"location":"API/#Base.write-Tuple{String,MeshArrays.gcmarray}","page":"API","title":"Base.write","text":"write(fil::String,x::MeshArray)\n\nWrite MeshArray to binary file. Other methods:\n\nwrite(xx::Array,x::MeshArray) #to Array\n\n\n\n\n\n","category":"method"},{"location":"API/#MeshArrays.MapWetPoints-Tuple{Any}","page":"API","title":"MeshArrays.MapWetPoints","text":"MapWetPoints(mskWet)\n\nMapping from global array to global ocean vector.\n\n(Kvec,Lvec,Kmap,Lmap)=MapWetPoints(mskWet)\n\n\n\n\n\n","category":"method"},{"location":"API/#MeshArrays.MaskWetPoints-Tuple{Any}","page":"API","title":"MeshArrays.MaskWetPoints","text":"MaskWetPoints(TrspCon)\n\nMask land points with NaN.\n\n(TrspCon, mskWet, mskDry)=MaskWetPoints(TrspCon)\n\n\n\n\n\n","category":"method"},{"location":"API/#MeshArrays.MatrixForPoisson-NTuple{7,Any}","page":"API","title":"MeshArrays.MatrixForPoisson","text":"MatrixForPoisson(TrspCon,mskWet,mskDry,Kvec,Lvec,Kmap,Lmap)\n\nAssemble sparse matrix using mskWet, Kvec, Lvec directly and Kmap, Lmap via SeedWetPoints\n\nA=MatrixForPoisson(TrspCon,mskWet,mskDry,Kvec,Lvec,Kmap,Lmap)\n\n\n\n\n\n","category":"method"},{"location":"API/#MeshArrays.ParaCoeffs","page":"API","title":"MeshArrays.ParaCoeffs","text":"ParaCoeffs(px,py,ox=[],oy=[])\n\nCompute bilinear interpolation coefficients for ox,oy within px,py by remapping these parallelograms to the unit square.\n\npx,py are Mx4 matrices where each line specifies one quadrilateral.\nox,oy are MxP position matrices\npw (output) are the MxPx4 bilinear interpolation weights\n\nx=1.0; y=1.0 #Try send the corners to unit square corners?\nprintln(vec(ParaCoeffs([0., 2., 3., 1.]',[0., 0., 1., 1.]',x,y)))\nprintln(vec(QuadCoeffs([0., 2.01, 3., 1.]',[0., 0., 1., 1.]',x,y)))\n\n\n\n\n\n","category":"function"},{"location":"API/#MeshArrays.PolygonAngle","page":"API","title":"MeshArrays.PolygonAngle","text":"PolygonAngle(px,py,x=[],y=[])\n\nCompute sum of interior angles for polygons or points-to-polygons (when px,py,x,y is provided as input). px,py are MxN matrices where each line specifies one polygon. (optional) x,y are position vectors.\n\npx=[0. 0. 1. 1.]; py=[0. 1. 1. 0.];\nx=collect(-1.0:0.25:2.0); y=x;\nPolygonAngle(px,py,x,y)\n\n\n\n\n\n","category":"function"},{"location":"API/#MeshArrays.QuadArrays-Union{Tuple{T}, Tuple{Array{T,2},Array{T,2}}} where T","page":"API","title":"MeshArrays.QuadArrays","text":"QuadArrays(x_grid,y_grid)\n\nTransform xgrid,ygrid (size ni+2,nj+2) into xquad,yquad,iquad,jquad quadrilaterals (size ni+1*nj+1,4) where iquad,jquad are point indices\n\n\n\n\n\n","category":"method"},{"location":"API/#MeshArrays.QuadCoeffs","page":"API","title":"MeshArrays.QuadCoeffs","text":"QuadCoeffs(px,py,ox=[],oy=[])\n\nCompute bilinear interpolation coefficients for ox,oy within px,py by remapping these quadrilaterals to the unit square.\n\npx,py are Mx4 matrices where each line specifies one quadrilateral.\nox,oy are MxP position matrices\now (output) are the MxPx4 bilinear interpolation weights\n\nQuadCoeffs([-1., 8., 13., -4.]',[-1., 3., 11., 8.]',0.,6.)\nQuadCoeffs([0., 2., 3., 1.]',[0., 0., 1., 1.]',0.1,0.1)\n\n\n\n\n\n","category":"function"},{"location":"API/#MeshArrays.SeedWetPoints-Tuple{MeshArrays.gcmarray,MeshArrays.gcmarray,MeshArrays.gcmarray,Vararg{Any,N} where N}","page":"API","title":"MeshArrays.SeedWetPoints","text":"SeedWetPoints(tmp::MeshArray,Kmap::MeshArray,Lmap::MeshArray,I...)\n\nSeed a subset of grid points.\n\n(FLDones,FLDkkFROM)=SeedWetPoints(tmp::MeshArray,Kmap::MeshArray,Lmap::MeshArray,I...)\n\n\n\n\n\n","category":"method"},{"location":"API/#MeshArrays.fijind-Tuple{MeshArrays.gcmfaces,Int64}","page":"API","title":"MeshArrays.fijind","text":"fijind(A::gcmfaces,ij::Int)\n\nCompute face and local indices (f,j,k) from global index (ij).\n\n(needed in other types?)\n\n\n\n\n\n","category":"method"},{"location":"API/#MeshArrays.getindexetc-Union{Tuple{N}, Tuple{T}, Tuple{MeshArrays.gcmarray{T,N,Array{T,2}},Vararg{Union{Colon, Int64, AbstractUnitRange, Array{Int64,N} where N},N}}} where N where T","page":"API","title":"MeshArrays.getindexetc","text":"getindexetc(A::gcmarray, I::Vararg{_}) where {T,N}\n\nSame as getindex but also returns the face size and index\n\n\n\n\n\n","category":"method"},{"location":"API/#MeshArrays.nFacesEtc-Tuple{MeshArrays.gcmarray}","page":"API","title":"MeshArrays.nFacesEtc","text":"nFacesEtc(a::gcmarray)\n\nReturn nFaces, n3 (1 in 2D case; >1 otherwise)\n\n\n\n\n\n","category":"method"}]
}