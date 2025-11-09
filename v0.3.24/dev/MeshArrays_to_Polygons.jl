### A Pluto.jl notebook ###
# v0.20.20

using Markdown
using InteractiveUtils

# This Pluto notebook uses @bind for interactivity. When running this notebook outside of Pluto, the following 'mock version' of @bind gives bound variables a default value (instead of an error).
macro bind(def, element)
    #! format: off
    return quote
        local iv = try Base.loaded_modules[Base.PkgId(Base.UUID("6e696c72-6542-2067-7265-42206c756150"), "AbstractPlutoDingetjes")].Bonds.initial_value catch; b -> missing; end
        local el = $(esc(element))
        global $(esc(def)) = Core.applicable(Base.get, el) ? Base.get(el) : iv(el)
        el
    end
    #! format: on
end

# ╔═╡ 540641e2-14e7-47da-9f8e-ac37cd499886
begin
	using MeshArrays, MITgcm, Glob
	using CairoMakie
	using PlutoUI
	import GeoInterface as GI
	import GeometryOps as GO
	import GeoInterfaceMakie
	md"""### Julia packages"""
end

# ╔═╡ 6da3809b-4bf4-45d6-bbd2-382312e947b7
TableOfContents()

# ╔═╡ 7aca186c-efd7-4312-afd2-37fdc1a2bbb6
md"""# Polygons, 3D LineStrings, and Visualization

In this example we focus on the cube sphere grid at coarse resolution (`cs32`).
"""

# ╔═╡ 3e5aac32-c82a-45e2-8d7b-05cf96a4ce3c
begin
	ui_do_sphere = @bind do_sphere CheckBox(default=true)
	md"""Show on the sphere? (or 2D plane view)

	$(ui_do_sphere)
	"""
end

# ╔═╡ 7bace9d8-4e9e-451d-aef4-5093c62b08e7
md"""## Main Calculations"""

# ╔═╡ 53b62449-7006-4f2a-a337-0b0e45d82f5e
function to_sphere(XC,YC)
	x=sin.(pi/2 .-YC*pi/180).*cos.(XC*pi/180);
	y=sin.(pi/2 .-YC*pi/180).*sin.(XC*pi/180);
	z=cos.(pi/2 .-YC*pi/180);
	(x,y,z)
end

# ╔═╡ 65dcb624-7f27-4247-bc73-446bf069a67a
function treat_180lon!(x::Matrix{Float64})
	if (maximum(x)-minimum(x))>180
		x[x.>0].=x[x.>0].-360
	end
end

# ╔═╡ 7d0afd5e-fa88-4f40-a1f0-1a2a944341cd
function to_LineStrings3D(pol)
	arr2=Array{Any}(undef,32,32)
	b0=Array{Any}(undef,5)
	for ij in eachindex(pol)
		geom=GI.getgeom(pol[ij])
		a=GI.getpoint.(geom)
		b=GI.coordinates.(geom)
		b1=GI.coordinates.(a[1])
		for p in 1:5
#			b0[p,:].=to_sphere(GI.coordinates(a[1][p])...)
			b0[p]=GI.Point(to_sphere(GI.coordinates(a[1][p])...))
		end
#		arr2[ij]=deepcopy(b0)
		arr2[ij]=deepcopy(GI.LineString(b0))
	end
	arr2
end

# ╔═╡ 915a5f42-9a31-4837-a891-d22f8825930b
md"""## Read Grid Definition (MITgcm input)

- should have a data structure for these
- for now just use an array of arrays
- in the plots : 
   - `YG` is latitude
   - `XG` is longitude
"""

# ╔═╡ 3c2f45b3-b3c3-4704-ba0e-c1ada7c48c13
begin
	path_MITgcm=MITgcm.getdata("mitgcmsmallverif")
	path_grid=joinpath(path_MITgcm,"MITgcm","verification","tutorial_held_suarez_cs","input")
	files_grid=glob("grid_cs32.face00?.bin",path_grid)
end

# ╔═╡ 8f3ffcf6-9aa8-4882-b055-583ff2bd82f4
begin
	list_fields=["XC","YC","DXF","DYF","RAC","XG","YG","DXV","DYU","RAZ",
		"DXC","DYC","RAW","RAS","DXG","DYG","AngleCS","AngleSN"]

#	ni=32; nj=32
#	list_ni=[ni,ni,ni,ni,ni,ni+1,ni+1,ni+1,ni+1,ni+1,ni+1,ni,ni+1,ni,ni,ni+1]
#	list_nj=[nj,nj,nj,nj,nj,nj+1,nj+1,nj+1,nj+1,nj+1,nj,nj+1,nj,nj+1,nj+1,nj]
#	sum(list_ni.*list_nj)
	
	function read_native_grid(ff=1,va="XG")
	  fil=files_grid[ff]
	  fid = open(fil)
      (n1,n2,n3)=(33,33,18)
	  xx = Array{Float64,1}(undef,n1*n2*n3)
	  read!(fid,xx)
	  xx = reshape(hton.(xx),(n1,n2,n3))
	  close(fid)
      vv=findall(list_fields.==va)[1]
      xx[:,:,vv]
	end

	#156816/(33*33)/8
	#length(list_fields)

	native_grid = begin
		XG=Array{Any}(undef,6)
		YG=Array{Any}(undef,6)
		for ff=1:6
			XG[ff]=read_native_grid(ff,"XG")
			YG[ff]=read_native_grid(ff,"YG")
		end
		(XG=XG,YG=YG)
	end
end

# ╔═╡ a248a2b1-ed63-48a0-922c-b00a3dbfc294
function to_Polygons(ff=1)
	arr2=Array{Any}(undef,32,32)
	IJ=1:32
	for i in IJ
		for j in IJ
			x=[XG[ff][i,j] XG[ff][i+1,j] XG[ff][i+1,j+1] XG[ff][i,j+1] XG[ff][i,j]]
			treat_180lon!(x)
			y=[YG[ff][i,j] YG[ff][i+1,j] YG[ff][i+1,j+1] YG[ff][i,j+1] YG[ff][i,j]]  
			arr2[i,j]=GI.Polygon([ GI.LinearRing(GI.Point.(zip(vec(x),vec(y)))) ])
#			arr2[i,j]=GI.LineString(GI.Point.(zip(vec(x),vec(y))))
		end
	end
	arr2
end

# ╔═╡ 17dd67f4-0c98-408a-9cf9-b02e9a00b5a3
pols=[to_Polygons(ff) for ff in 1:6];

# ╔═╡ 2fd20020-8b29-4262-92f3-aabf9bc609fa
pols3D=[to_LineStrings3D(p) for p in pols];

# ╔═╡ f112ecfa-8a66-4ae1-aac3-7b354bb5863f
let
	fi=Figure(); ax=Axis(fi[1,1])
	for ff in 1:6
		[plot!(p) for p in pols[ff]]
	end
	fi
end

# ╔═╡ 87bb371f-c08e-4a31-b7bf-56b58fe9cbda
function to_LineStrings2D(ff=1;do_sphere=true)
	arr2=Array{Any}(undef,32,32)
	IJ=1:32
	for i in IJ
		for j in IJ
			x=[XG[ff][i,j] XG[ff][i+1,j] XG[ff][i+1,j+1] XG[ff][i,j+1] XG[ff][i,j]] 
			treat_180lon!(x)
			y=[YG[ff][i,j] YG[ff][i+1,j] YG[ff][i+1,j+1] YG[ff][i,j+1] YG[ff][i,j]]  
			arr2[i,j]=GI.LineString(GI.Point.(zip(vec(x),vec(y))))
		end
	end
	arr2
end

# ╔═╡ 523311f3-8d0a-4b57-b376-da945c4f20a7
md"""## Read MITgcm Output for this Grid"""

# ╔═╡ 4208e083-2971-4b1c-bc72-5df1caaf3d36
begin
	γ=GridSpec(ID=:CS32)
	Γ=GridLoad(γ,option=:light)
end

# ╔═╡ 199f9dfc-c393-47f9-95d8-c0a783ae47e6
md"""# Appendix"""

# ╔═╡ 3549d2f2-ab0f-44a0-a0dd-b9aa64354f31
md"""### Plotting Functions"""

# ╔═╡ 6a354d1f-eba1-40f6-bdfc-1ecd8d96a0ef
function fig2_basis(do_sphere)
	f=Figure()
	if do_sphere
		scene = LScene(f[1, 1])
		cam = Makie.Camera3D(scene.scene, projectiontype = Makie.Perspective)
		sphere = Sphere(Point3f(0), 0.99f0)
		mesh!(scene,sphere, color = :white)
		zoom!(scene.scene,cam,2.0)
		a=scene
	else
		a=Axis(f[1,1])
	end	
	f,a
end

# ╔═╡ 798f41ed-cd77-4890-9148-be278df52ba4
function fig2(facets=1:6,cc=[:blue :green :orange :black :red :violet];
	do_sphere=true)
	f,a=fig2_basis(do_sphere)
	for ff in facets
		#here in 3D can we do LineRings or generalized Polygons?
		if do_sphere
#			xyz=to_sphere.(vec(Γ.XC[ff]),vec(Γ.YC[ff]))
#			scatter!(a,xyz,markersize=2,color=cc[ff])
			[plot!(a,pols3D[ff][i,j],color=cc[ff],linewidth=0.5) for i in 1:32, j in 1:32]
		else
			#arr=to_LineStrings2D(ff,do_sphere=false)
			[lines!(a,pols[ff][i,j],color=cc[ff],linewidth=0.5) for i in 1:32, j in 1:32]
		end
	end
	f
end

# ╔═╡ 1b5b7cdc-9528-4afc-b282-2828da34ae2f
f2=fig2(1:3,do_sphere=do_sphere)

# ╔═╡ 1effaab8-8237-4236-ad75-755795f898fc
function fig1(grid::NamedTuple,va=:XG,na="latitude",cr=(-90,90),cm=:thermal)
	f=Figure()
	ax=[[[1,ff] for ff in 1:3] [[2,ff] for ff in 1:3]]
	for ff in 1:length(ax)
		axi=Axis(f[ax[ff]...],title="$(va) ($(na)) face $(ff)")
		heatmap!(axi,grid[va][ff],colorrange=cr,colormap=cm)
	end
	Colorbar(f[1:2,4],colorrange = cr,colormap=cm)
	f
end

# ╔═╡ a5fcb7cb-a86c-4f19-a3b9-1e280ce6159b
fig1(native_grid,:YG)

# ╔═╡ ac64b4cf-a6b5-4dad-a90c-bf1f483690e5
fig1(native_grid,:XG,"longitude",(-180,180))

# ╔═╡ 98496427-b774-4d9b-b386-ffe5ea775d37
function fig0()
	f=Figure()
	a=Axis3(f[1,1])
	sphere = Sphere(Point3f(0), 0.99f0)
	mesh!(a,sphere, color = :white)
	cc=[:blue :green :orange :black :red :violet]
	for f in 1:6
		XC=Γ.XC[f][:]
		YC=Γ.YC[f][:]
		scatter!(to_sphere(XC,YC)...,markersize=5.0,color=cc[f])
	end
	f
end

# ╔═╡ fc2e3116-8a7e-405d-ad6d-bdf5d4ae45c9
md"""### Incomplete Solution

- not really useful 
- using grid output does not really work for this
- incl due missing corner data points in CS grid output
- even if vorticity point was available (as in gcmfaces)
- need to use the input to MITgcm instead (33x33 v 32x32)
"""

# ╔═╡ f4d0e938-e544-4add-be50-7942df72a3b8
incomplete_solution="""
	arr=get_LineString()
	arr=Array{Any}(undef,32,32)
	for i in 2:33
	for j in 2:33
	if j<33
		x=[XGe.MA[1][i,j] XGe.MA[1][i+1,j] XGe.MA[1][i+1,j+1] XGe.MA[1][i,j+1] XGe.MA[1][i,j]] 
		y=[YGe.MA[1][i,j] YGe.MA[1][i+1,j] YGe.MA[1][i+1,j+1] YGe.MA[1][i,j+1] YGe.MA[1][i,j]]  
	else #should be exch_Z
		x=[XGe.MA[1][i,j] XGe.MA[1][i+1,j] XGe.MA[1][i+1-1,j+1] XGe.MA[1][i-1,j+1] XGe.MA[1][i,j]]
		y=[YGe.MA[1][i,j] YGe.MA[1][i+1,j] YGe.MA[1][i+1-1,j+1] YGe.MA[1][i-1,j+1] YGe.MA[1][i,j]]
	end
		#arr[i,j]=zip.(x,y)
		arr[i-1,j-1]=GI.LineString(GI.Point.(zip(vec(x),vec(y))))
	end
	end
"""

# ╔═╡ 5f5e6be3-a7ad-4772-a52c-acb765106593
begin 
	XCe=exchange(Γ.XC)
	YCe=exchange(Γ.YC)
	XGe=exchange(Γ.XG)
	YGe=exchange(Γ.YG)
	all_points=(XCe,YCe,XGe,YGe)
	md"Exchange function calls for incomplete solution"
end

# ╔═╡ a4418a2e-831b-416f-8348-e282b08cfa2d
md"""### Basic Tests

- `LinearRing` does not seem to work in 3D (`do_two_D=false`)
- `LineString` works in 3D. 
- `(x,y,z)` is wrapped in a couple levels of types. 
"""

# ╔═╡ 5365dc71-0ddc-4d40-8f38-7c494b3305d9
begin
	ui_do_two_D = @bind do_two_D CheckBox(default=true)
	md""" `do_two_D` = $(ui_do_two_D) (false => 3D => error message)"""
end

# ╔═╡ b494d021-3727-4736-8031-9587d0d542a2
let #Let's try constructing from scratch again
	a=( (0.5773502691896258,-0.5773502691896257,-0.5773502691896257),
		(0.5837494212853372,-0.5643343213894291,-0.5837494212853374),
		(0.5977118906557506,-0.566895270649142,-0.566895270649142),
		(0.5837494212853375,-0.5837494212853372,-0.5643343213894292),
		(0.5773502691896258,-0.5773502691896257,-0.5773502691896257),
		)
	#do_two_D=false
	if do_two_D
		aa=[b[1:2] for b in a]
		bb=GI.LinearRing(aa)
		#	GI.polygon([bb])
		#	GI.LinearRing(a)
		#	aa,GI.getpoint(bb)
		cc=GI.Polygon([bb])
		plot(cc,color=:cyan)
		lines!(aa,color=:red)
		current_figure()
	else
		bb=GI.LinearRing(a)
	end
end	

# ╔═╡ e4d61551-c0e7-4a61-b38e-774897862fed
typeof(to_LineStrings2D(1)[1,1])

# ╔═╡ 50272674-6935-4239-89bf-29997341148b
b=let #Let's dig into a LineString
	arr=to_LineStrings2D(1)
#	polygon1 = GI.Polygon([arr]);
	typeof(arr[1,1])
	GI.getgeom(arr[1,1])
	a=GI.getpoint(arr[1,1])
	###GI.LinearRing(a)
	b=a[1]
	c=b.geom[1:end]
	(a,b,c)
end


# ╔═╡ 63c58499-c353-4124-9270-6c22b0e6bcd7
let
	#[GO.area(pol0) GO.area(pol1)]
	#GO.intersection(pol0, pol1)
	#typeof(pol1)
	line1 = GI.Line([(124.584961,-12.768946), (126.738281,-17.224758)])
	line2 = GI.Line([(123.354492,-15.961329), (127.22168,-14.008696)])
	inter_points = GO.intersection(line1, line2; target = GI.PointTrait())
	GI.coordinates.(inter_points)
end

# ╔═╡ 9405eb97-1130-427d-a9f4-e45efea361b2
let
	#GI.getgeom.
	geom=GI.getgeom(pols[1][1,1])
#	ext = exterior(geom)
	a=GI.getpoint.(geom)
	b=GI.coordinates.(geom)
	b1=GI.coordinates.(a[1])
	b0=zeros(5,2)
	for p in 1:5
		b0[p,:].=GI.coordinates(a[1][p])
	end
	b0
end

# ╔═╡ 00000000-0000-0000-0000-000000000001
PLUTO_PROJECT_TOML_CONTENTS = """
[deps]
GeometryOps = "3251bfac-6a57-4b6d-aa61-ac1fef2975ab"
Glob = "c27321d9-0574-5035-807b-f59d2c89b15c"
MITgcm = "dce5fa8e-68ce-4431-a242-9469c69627a0"
MeshArrays = "cb8c808f-1acf-59a3-9d2b-6e38d009f683"
PlutoUI = "7f904dfe-b85e-4ff6-b463-dae2292396a8"

[compat]
GeometryOps = "~0.1.13"
Glob = "~1.3.1"
MITgcm = "~0.5.0"
MeshArrays = "~0.3.17"
PlutoUI = "~0.7.60"
"""

# ╔═╡ 00000000-0000-0000-0000-000000000002
PLUTO_MANIFEST_TOML_CONTENTS = """
# This file is machine-generated - editing it directly is not advised

julia_version = "1.12.1"
manifest_format = "2.0"
project_hash = "3652a10fe9714d25bb73f7bab17a3139a7b6dfd7"

[[deps.AbstractPlutoDingetjes]]
deps = ["Pkg"]
git-tree-sha1 = "6e1d2a35f2f90a4bc7c2ed98079b2ba09c35b83a"
uuid = "6e696c72-6542-2067-7265-42206c756150"
version = "1.3.2"

[[deps.AbstractTrees]]
git-tree-sha1 = "2d9c9a55f9c93e8887ad391fbae72f8ef55e1177"
uuid = "1520ce14-60c1-5f80-bbc7-55ef81b5835c"
version = "0.4.5"

[[deps.AdaptivePredicates]]
git-tree-sha1 = "7e651ea8d262d2d74ce75fdf47c4d63c07dba7a6"
uuid = "35492f91-a3bd-45ad-95db-fcad7dcfedb7"
version = "1.2.0"

[[deps.ArgTools]]
uuid = "0dad84c5-d112-42e6-8d28-ef12dabb789f"
version = "1.1.2"

[[deps.Artifacts]]
uuid = "56f22d72-fd6d-98f1-02f0-08ddc0907c33"
version = "1.11.0"

[[deps.Base64]]
uuid = "2a0f44e3-6c83-55bd-87e4-b1978d98bd5f"
version = "1.11.0"

[[deps.BitFlags]]
git-tree-sha1 = "0691e34b3bb8be9307330f88d1a3c3f25466c24d"
uuid = "d1d4a3ce-64b1-5f1a-9ba4-7e7e69966f35"
version = "0.1.9"

[[deps.CFTime]]
deps = ["Dates", "Printf"]
git-tree-sha1 = "2ed76cf4ac70526e3df565435d65e7c7b5c7a77a"
uuid = "179af706-886a-5703-950a-314cd64e0468"
version = "0.2.4"

[[deps.CRlibm]]
deps = ["CRlibm_jll"]
git-tree-sha1 = "66188d9d103b92b6cd705214242e27f5737a1e5e"
uuid = "96374032-68de-5a5b-8d9e-752f78720389"
version = "1.0.2"

[[deps.CRlibm_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "e329286945d0cfc04456972ea732551869af1cfc"
uuid = "4e9b3aee-d8a1-5a3d-ad8b-7d824db253f0"
version = "1.0.1+0"

[[deps.CSV]]
deps = ["CodecZlib", "Dates", "FilePathsBase", "InlineStrings", "Mmap", "Parsers", "PooledArrays", "PrecompileTools", "SentinelArrays", "Tables", "Unicode", "WeakRefStrings", "WorkerUtilities"]
git-tree-sha1 = "deddd8725e5e1cc49ee205a1964256043720a6c3"
uuid = "336ed68f-0bac-5ca0-87d4-7b16caf5d00b"
version = "0.10.15"

[[deps.CatViews]]
deps = ["Random", "Test"]
git-tree-sha1 = "23d1f1e10d4e24374112fcf800ac981d14a54b24"
uuid = "81a5f4ea-a946-549a-aa7e-2a7f63a27d31"
version = "1.0.0"

[[deps.ChunkCodecCore]]
git-tree-sha1 = "51f4c10ee01bda57371e977931de39ee0f0cdb3e"
uuid = "0b6fb165-00bc-4d37-ab8b-79f91016dbe1"
version = "1.0.0"

[[deps.ChunkCodecLibZlib]]
deps = ["ChunkCodecCore", "Zlib_jll"]
git-tree-sha1 = "cee8104904c53d39eb94fd06cbe60cb5acde7177"
uuid = "4c0bbee4-addc-4d73-81a0-b6caacae83c8"
version = "1.0.0"

[[deps.ChunkCodecLibZstd]]
deps = ["ChunkCodecCore", "Zstd_jll"]
git-tree-sha1 = "34d9873079e4cb3d0c62926a225136824677073f"
uuid = "55437552-ac27-4d47-9aa3-63184e8fd398"
version = "1.0.0"

[[deps.ClimateModels]]
deps = ["CFTime", "CSV", "DataDeps", "DataFrames", "Dataverse", "Dates", "Downloads", "Git", "Glob", "JLD2", "OffsetArrays", "OrderedCollections", "Pkg", "Printf", "Random", "Statistics", "Suppressor", "TOML", "Test", "UUIDs"]
git-tree-sha1 = "0d1f3ac160e3b22ba4bf6e28159fd0508905f465"
uuid = "f6adb021-9183-4f40-84dc-8cea6f651bb0"
version = "0.3.11"

    [deps.ClimateModels.extensions]
    ClimateModelsCondaExt = ["Conda"]
    ClimateModelsIniFileExt = ["IniFile"]
    ClimateModelsMakieExt = ["Makie"]
    ClimateModelsNetCDFExt = ["NetCDF"]
    ClimateModelsOceananigansExt = ["Oceananigans"]
    ClimateModelsPyCallExt = ["PyCall"]
    ClimateModelsZarrExt = ["Zarr"]

    [deps.ClimateModels.weakdeps]
    Conda = "8f4d0f93-b110-5947-807f-2305c1781a2d"
    IniFile = "83e8ac13-25f8-5344-8a64-a9f2b223428f"
    Makie = "ee78f7c6-11fb-53f2-987a-cfe4a2b5a57a"
    NetCDF = "30363a11-5582-574a-97bb-aa9a979735b9"
    Oceananigans = "9e8cae18-63c1-5223-a75c-80ca9d6e9a09"
    PyCall = "438e738f-606a-5dbb-bf0a-cddfbfd45ab0"
    Zarr = "0a941bbe-ad1d-11e8-39d9-ab76183a1d99"

[[deps.CodecZlib]]
deps = ["TranscodingStreams", "Zlib_jll"]
git-tree-sha1 = "962834c22b66e32aa10f7611c08c8ca4e20749a9"
uuid = "944b1d66-785c-5afd-91f1-9de20f533193"
version = "0.7.8"

[[deps.ColorTypes]]
deps = ["FixedPointNumbers", "Random"]
git-tree-sha1 = "67e11ee83a43eb71ddc950302c53bf33f0690dfe"
uuid = "3da002f7-5984-5a60-b8a6-cbb66c0b333f"
version = "0.12.1"
weakdeps = ["StyledStrings"]

    [deps.ColorTypes.extensions]
    StyledStringsExt = "StyledStrings"

[[deps.Compat]]
deps = ["TOML", "UUIDs"]
git-tree-sha1 = "9d8a54ce4b17aa5bdce0ea5c34bc5e7c340d16ad"
uuid = "34da2185-b29b-5c13-b0c7-acf172513d20"
version = "4.18.1"
weakdeps = ["Dates", "LinearAlgebra"]

    [deps.Compat.extensions]
    CompatLinearAlgebraExt = "LinearAlgebra"

[[deps.CompilerSupportLibraries_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "e66e0078-7015-5450-92f7-15fbd957f2ae"
version = "1.3.0+1"

[[deps.ConcurrentUtilities]]
deps = ["Serialization", "Sockets"]
git-tree-sha1 = "d9d26935a0bcffc87d2613ce14c527c99fc543fd"
uuid = "f0e56b4a-5159-44fe-b623-3e5288b988bb"
version = "2.5.0"

[[deps.CoordinateTransformations]]
deps = ["LinearAlgebra", "StaticArrays"]
git-tree-sha1 = "a692f5e257d332de1e554e4566a4e5a8a72de2b2"
uuid = "150eb455-5306-5404-9cee-2592286d6298"
version = "0.6.4"

[[deps.Crayons]]
git-tree-sha1 = "249fe38abf76d48563e2f4556bebd215aa317e15"
uuid = "a8cc5b0e-0ffa-5ad4-8c14-923d3ee1735f"
version = "4.1.1"

[[deps.DataAPI]]
git-tree-sha1 = "abe83f3a2f1b857aac70ef8b269080af17764bbe"
uuid = "9a962f9c-6df0-11e9-0e5d-c546b8b5ee8a"
version = "1.16.0"

[[deps.DataDeps]]
deps = ["HTTP", "Libdl", "Reexport", "SHA", "Scratch", "p7zip_jll"]
git-tree-sha1 = "8ae085b71c462c2cb1cfedcb10c3c877ec6cf03f"
uuid = "124859b0-ceae-595e-8997-d05f6a7a8dfe"
version = "0.7.13"

[[deps.DataFrames]]
deps = ["Compat", "DataAPI", "DataStructures", "Future", "InlineStrings", "InvertedIndices", "IteratorInterfaceExtensions", "LinearAlgebra", "Markdown", "Missings", "PooledArrays", "PrecompileTools", "PrettyTables", "Printf", "Random", "Reexport", "SentinelArrays", "SortingAlgorithms", "Statistics", "TableTraits", "Tables", "Unicode"]
git-tree-sha1 = "d8928e9169ff76c6281f39a659f9bca3a573f24c"
uuid = "a93c6f00-e57d-5684-b7b6-d8193f3e46c0"
version = "1.8.1"

[[deps.DataStructures]]
deps = ["OrderedCollections"]
git-tree-sha1 = "e357641bb3e0638d353c4b29ea0e40ea644066a6"
uuid = "864edb3b-99cc-5e75-8d2d-829cb0a9cfe8"
version = "0.19.3"

[[deps.DataValueInterfaces]]
git-tree-sha1 = "bfc1187b79289637fa0ef6d4436ebdfe6905cbd6"
uuid = "e2d170a0-9d28-54be-80f0-106bbe20a464"
version = "1.0.0"

[[deps.Dataverse]]
deps = ["CSV", "CodecZlib", "DataFrames", "Downloads", "HTTP", "JSON", "NetworkOptions", "Tar", "ZipFile"]
git-tree-sha1 = "84ad16bc3c21a8d13e07c55287f7049cda236349"
uuid = "9c0b9be8-e31e-490f-90fe-77697562404d"
version = "0.2.7"

    [deps.Dataverse.extensions]
    DataverseCondaPkgExt = ["CondaPkg"]
    DataversePythonCallExt = ["PythonCall"]

    [deps.Dataverse.weakdeps]
    CondaPkg = "992eb4ea-22a4-4c89-a5bb-47a3300528ab"
    PythonCall = "6099a3de-0909-46bc-b1f4-468b9a2dfc0d"

[[deps.Dates]]
deps = ["Printf"]
uuid = "ade2ca70-3891-5945-98fb-dc099432e06a"
version = "1.11.0"

[[deps.DelaunayTriangulation]]
deps = ["AdaptivePredicates", "EnumX", "ExactPredicates", "Random"]
git-tree-sha1 = "783b21581a051ac91a3921ee37e26a23ed7f57a6"
uuid = "927a84f5-c5f4-47a5-9785-b46e178433df"
version = "1.6.5"

[[deps.Distances]]
deps = ["LinearAlgebra", "Statistics", "StatsAPI"]
git-tree-sha1 = "c7e3a542b999843086e2f29dac96a618c105be1d"
uuid = "b4f34e82-e78d-54a5-968a-f98e89d6e8f7"
version = "0.10.12"

    [deps.Distances.extensions]
    DistancesChainRulesCoreExt = "ChainRulesCore"
    DistancesSparseArraysExt = "SparseArrays"

    [deps.Distances.weakdeps]
    ChainRulesCore = "d360d2e6-b24c-11e9-a2a3-2a2ae2dbcce4"
    SparseArrays = "2f01184e-e22b-5df5-ae63-d93ebab69eaf"

[[deps.Distributed]]
deps = ["Random", "Serialization", "Sockets"]
uuid = "8ba89e20-285c-5b6f-9357-94700520ee1b"
version = "1.11.0"

[[deps.Downloads]]
deps = ["ArgTools", "FileWatching", "LibCURL", "NetworkOptions"]
uuid = "f43a241f-c20a-4ad4-852c-f6b1247861c6"
version = "1.6.0"

[[deps.EnumX]]
git-tree-sha1 = "bddad79635af6aec424f53ed8aad5d7555dc6f00"
uuid = "4e289a0a-7415-4d19-859d-a7e5c4648b56"
version = "1.0.5"

[[deps.ExactPredicates]]
deps = ["IntervalArithmetic", "Random", "StaticArrays"]
git-tree-sha1 = "83231673ea4d3d6008ac74dc5079e77ab2209d8f"
uuid = "429591f6-91af-11e9-00e2-59fbe8cec110"
version = "2.2.9"

[[deps.ExceptionUnwrapping]]
deps = ["Test"]
git-tree-sha1 = "d36f682e590a83d63d1c7dbd287573764682d12a"
uuid = "460bff9d-24e4-43bc-9d9f-a8973cb893f4"
version = "0.1.11"

[[deps.Expat_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "27af30de8b5445644e8ffe3bcb0d72049c089cf1"
uuid = "2e619515-83b5-522b-bb60-26c02a35a201"
version = "2.7.3+0"

[[deps.Extents]]
git-tree-sha1 = "b309b36a9e02fe7be71270dd8c0fd873625332b4"
uuid = "411431e0-e8b7-467b-b5e0-f676ba4f2910"
version = "0.1.6"

[[deps.FileIO]]
deps = ["Pkg", "Requires", "UUIDs"]
git-tree-sha1 = "d60eb76f37d7e5a40cc2e7c36974d864b82dc802"
uuid = "5789e2e9-d7fb-5bc7-8068-2c6fae9b9549"
version = "1.17.1"
weakdeps = ["HTTP"]

    [deps.FileIO.extensions]
    HTTPExt = "HTTP"

[[deps.FilePathsBase]]
deps = ["Compat", "Dates"]
git-tree-sha1 = "3bab2c5aa25e7840a4b065805c0cdfc01f3068d2"
uuid = "48062228-2e41-5def-b9a4-89aafe57970f"
version = "0.9.24"
weakdeps = ["Mmap", "Test"]

    [deps.FilePathsBase.extensions]
    FilePathsBaseMmapExt = "Mmap"
    FilePathsBaseTestExt = "Test"

[[deps.FileWatching]]
uuid = "7b1f6079-737a-58dc-b8bc-7a2ca5c1b5ee"
version = "1.11.0"

[[deps.FixedPointNumbers]]
deps = ["Statistics"]
git-tree-sha1 = "05882d6995ae5c12bb5f36dd2ed3f61c98cbb172"
uuid = "53c48c17-4a7d-5ca2-90c5-79b7896eea93"
version = "0.8.5"

[[deps.FortranFiles]]
git-tree-sha1 = "97069e9106dffe888562d974acf1d225cfbf8d4e"
uuid = "c58ffaec-ab22-586d-bfc5-781a99fd0b10"
version = "0.6.2"

[[deps.Future]]
deps = ["Random"]
uuid = "9fa8497b-333b-5362-9e8d-4d0656e87820"
version = "1.11.0"

[[deps.GeoFormatTypes]]
git-tree-sha1 = "7528a7956248c723d01a0a9b0447bf254bf4da52"
uuid = "68eda718-8dee-11e9-39e7-89f7f65f511f"
version = "0.4.5"

[[deps.GeoInterface]]
deps = ["DataAPI", "Extents", "GeoFormatTypes"]
git-tree-sha1 = "b7c5cdf45298877bb683bdda3f871ff7070985c4"
uuid = "cf35fbd7-0cd7-5166-be24-54bfbe79505f"
version = "1.6.0"

    [deps.GeoInterface.extensions]
    GeoInterfaceMakieExt = ["Makie", "GeometryBasics"]
    GeoInterfaceRecipesBaseExt = "RecipesBase"

    [deps.GeoInterface.weakdeps]
    GeometryBasics = "5c1252a2-5f33-56bf-86c9-59e7332b4326"
    Makie = "ee78f7c6-11fb-53f2-987a-cfe4a2b5a57a"
    RecipesBase = "3cdcf5f2-1ef4-517c-9805-6587b60abb01"

[[deps.GeometryOps]]
deps = ["AbstractTrees", "AdaptivePredicates", "CoordinateTransformations", "DataAPI", "DelaunayTriangulation", "ExactPredicates", "Extents", "GeoFormatTypes", "GeoInterface", "GeometryOpsCore", "LinearAlgebra", "Random", "SortTileRecursiveTree", "StaticArrays", "Statistics", "Tables"]
git-tree-sha1 = "9fa16be9c28d9c01bf2b5d73f7768d482c12b118"
uuid = "3251bfac-6a57-4b6d-aa61-ac1fef2975ab"
version = "0.1.31"

    [deps.GeometryOps.extensions]
    GeometryOpsDataFramesExt = "DataFrames"
    GeometryOpsFlexiJoinsExt = "FlexiJoins"
    GeometryOpsLibGEOSExt = "LibGEOS"
    GeometryOpsProjExt = "Proj"
    GeometryOpsTGGeometryExt = "TGGeometry"

    [deps.GeometryOps.weakdeps]
    DataFrames = "a93c6f00-e57d-5684-b7b6-d8193f3e46c0"
    FlexiJoins = "e37f2e79-19fa-4eb7-8510-b63b51fe0a37"
    LibGEOS = "a90b1aa1-3769-5649-ba7e-abc5a9d163eb"
    Proj = "c94c279d-25a6-4763-9509-64d165bea63e"
    TGGeometry = "d7e755d2-3c95-4bcf-9b3c-79ab1a78647b"

[[deps.GeometryOpsCore]]
deps = ["DataAPI", "GeoInterface", "StableTasks", "Tables"]
git-tree-sha1 = "69fc98947b06f8ac4279cf5bf8810373fe042be4"
uuid = "05efe853-fabf-41c8-927e-7063c8b9f013"
version = "0.1.7"

[[deps.Git]]
deps = ["Git_LFS_jll", "Git_jll", "JLLWrappers", "OpenSSH_jll"]
git-tree-sha1 = "824a1890086880696fc908fe12a17bcf61738bd8"
uuid = "d7ba0133-e1db-5d97-8f8c-041e4b3a1eb2"
version = "1.5.0"

[[deps.Git_LFS_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "bb8471f313ed941f299aa53d32a94ab3bee08844"
uuid = "020c3dae-16b3-5ae5-87b3-4cb189e250b2"
version = "3.7.0+0"

[[deps.Git_jll]]
deps = ["Artifacts", "Expat_jll", "JLLWrappers", "LibCURL_jll", "Libdl", "Libiconv_jll", "OpenSSL_jll", "PCRE2_jll", "Zlib_jll"]
git-tree-sha1 = "b6a684587ebe896d9f68ae777f648205940f0f70"
uuid = "f8c6e375-362e-5223-8a59-34ff63f689eb"
version = "2.51.3+0"

[[deps.Glob]]
git-tree-sha1 = "97285bbd5230dd766e9ef6749b80fc617126d496"
uuid = "c27321d9-0574-5035-807b-f59d2c89b15c"
version = "1.3.1"

[[deps.HTTP]]
deps = ["Base64", "CodecZlib", "ConcurrentUtilities", "Dates", "ExceptionUnwrapping", "Logging", "LoggingExtras", "MbedTLS", "NetworkOptions", "OpenSSL", "PrecompileTools", "Random", "SimpleBufferStream", "Sockets", "URIs", "UUIDs"]
git-tree-sha1 = "5e6fe50ae7f23d171f44e311c2960294aaa0beb5"
uuid = "cd3eb016-35fb-5094-929b-558a96fad6f3"
version = "1.10.19"

[[deps.HashArrayMappedTries]]
git-tree-sha1 = "2eaa69a7cab70a52b9687c8bf950a5a93ec895ae"
uuid = "076d061b-32b6-4027-95e0-9a2c6f6d7e74"
version = "0.2.0"

[[deps.Hyperscript]]
deps = ["Test"]
git-tree-sha1 = "179267cfa5e712760cd43dcae385d7ea90cc25a4"
uuid = "47d2ed2b-36de-50cf-bf87-49c2cf4b8b91"
version = "0.0.5"

[[deps.HypertextLiteral]]
deps = ["Tricks"]
git-tree-sha1 = "7134810b1afce04bbc1045ca1985fbe81ce17653"
uuid = "ac1192a8-f4b3-4bfe-ba22-af5b92cd3ab2"
version = "0.9.5"

[[deps.IOCapture]]
deps = ["Logging", "Random"]
git-tree-sha1 = "0ee181ec08df7d7c911901ea38baf16f755114dc"
uuid = "b5f81e59-6552-4d32-b1f0-c071b021bf89"
version = "1.0.0"

[[deps.InlineStrings]]
git-tree-sha1 = "8f3d257792a522b4601c24a577954b0a8cd7334d"
uuid = "842dd82b-1e85-43dc-bf29-5d0ee9dffc48"
version = "1.4.5"

    [deps.InlineStrings.extensions]
    ArrowTypesExt = "ArrowTypes"
    ParsersExt = "Parsers"

    [deps.InlineStrings.weakdeps]
    ArrowTypes = "31f734f8-188a-4ce0-8406-c8a06bd891cd"
    Parsers = "69de0a69-1ddd-5017-9359-2bf0b02dc9f0"

[[deps.InteractiveUtils]]
deps = ["Markdown"]
uuid = "b77e0a4c-d291-57a0-90e8-8db25a27a240"
version = "1.11.0"

[[deps.IntervalArithmetic]]
deps = ["CRlibm", "MacroTools", "OpenBLASConsistentFPCSR_jll", "Printf", "Random", "RoundingEmulator"]
git-tree-sha1 = "bf0210c01fb7d67c31fed97d7c1d1716b98ea689"
uuid = "d1acc4aa-44c8-5952-acd4-ba5d80a2a253"
version = "1.0.1"

    [deps.IntervalArithmetic.extensions]
    IntervalArithmeticArblibExt = "Arblib"
    IntervalArithmeticDiffRulesExt = "DiffRules"
    IntervalArithmeticForwardDiffExt = "ForwardDiff"
    IntervalArithmeticIntervalSetsExt = "IntervalSets"
    IntervalArithmeticLinearAlgebraExt = "LinearAlgebra"
    IntervalArithmeticRecipesBaseExt = "RecipesBase"
    IntervalArithmeticSparseArraysExt = "SparseArrays"

    [deps.IntervalArithmetic.weakdeps]
    Arblib = "fb37089c-8514-4489-9461-98f9c8763369"
    DiffRules = "b552c78f-8df3-52c6-915a-8e097449b14b"
    ForwardDiff = "f6369f11-7733-5829-9624-2563aa707210"
    IntervalSets = "8197267c-284f-5f27-9208-e0e47529a953"
    LinearAlgebra = "37e2e46d-f89d-539d-b4ee-838fcccc9c8e"
    RecipesBase = "3cdcf5f2-1ef4-517c-9805-6587b60abb01"
    SparseArrays = "2f01184e-e22b-5df5-ae63-d93ebab69eaf"

[[deps.InvertedIndices]]
git-tree-sha1 = "6da3c4316095de0f5ee2ebd875df8721e7e0bdbe"
uuid = "41ab1584-1d38-5bbf-9106-f11c6c58b48f"
version = "1.3.1"

[[deps.IteratorInterfaceExtensions]]
git-tree-sha1 = "a3f24677c21f5bbe9d2a714f95dcd58337fb2856"
uuid = "82899510-4779-5014-852e-03e436cf321d"
version = "1.0.0"

[[deps.JLD2]]
deps = ["ChunkCodecLibZlib", "ChunkCodecLibZstd", "FileIO", "MacroTools", "Mmap", "OrderedCollections", "PrecompileTools", "ScopedValues"]
git-tree-sha1 = "da2e9b4d1abbebdcca0aa68afa0aa272102baad7"
uuid = "033835bb-8acc-5ee8-8aae-3f567f8a3819"
version = "0.6.2"

    [deps.JLD2.extensions]
    UnPackExt = "UnPack"

    [deps.JLD2.weakdeps]
    UnPack = "3a884ed6-31ef-47d7-9d2a-63182c4928ed"

[[deps.JLLWrappers]]
deps = ["Artifacts", "Preferences"]
git-tree-sha1 = "0533e564aae234aff59ab625543145446d8b6ec2"
uuid = "692b3bcd-3c85-4b1f-b108-f13ce0eb3210"
version = "1.7.1"

[[deps.JSON]]
deps = ["Dates", "Mmap", "Parsers", "Unicode"]
git-tree-sha1 = "31e996f0a15c7b280ba9f76636b3ff9e2ae58c9a"
uuid = "682c06a0-de6a-54ab-a142-c8b1cf79cde6"
version = "0.21.4"

[[deps.JuliaSyntaxHighlighting]]
deps = ["StyledStrings"]
uuid = "ac6e5ff7-fb65-4e79-a425-ec3bc9c03011"
version = "1.12.0"

[[deps.LaTeXStrings]]
git-tree-sha1 = "dda21b8cbd6a6c40d9d02a73230f9d70fed6918c"
uuid = "b964fa9f-0449-5b57-a5c2-d3ea65f4040f"
version = "1.4.0"

[[deps.LazyArtifacts]]
deps = ["Artifacts", "Pkg"]
uuid = "4af54fe1-eca0-43a8-85a7-787d91b784e3"
version = "1.11.0"

[[deps.LibCURL]]
deps = ["LibCURL_jll", "MozillaCACerts_jll"]
uuid = "b27032c2-a3e7-50c8-80cd-2d36dbcbfd21"
version = "0.6.4"

[[deps.LibCURL_jll]]
deps = ["Artifacts", "LibSSH2_jll", "Libdl", "OpenSSL_jll", "Zlib_jll", "nghttp2_jll"]
uuid = "deac9b47-8bc7-5906-a0fe-35ac56dc84c0"
version = "8.11.1+1"

[[deps.LibGit2]]
deps = ["LibGit2_jll", "NetworkOptions", "Printf", "SHA"]
uuid = "76f85450-5226-5b5a-8eaa-529ad045b433"
version = "1.11.0"

[[deps.LibGit2_jll]]
deps = ["Artifacts", "LibSSH2_jll", "Libdl", "OpenSSL_jll"]
uuid = "e37daf67-58a4-590a-8e99-b0245dd2ffc5"
version = "1.9.0+0"

[[deps.LibSSH2_jll]]
deps = ["Artifacts", "Libdl", "OpenSSL_jll"]
uuid = "29816b5a-b9ab-546f-933c-edad1886dfa8"
version = "1.11.3+1"

[[deps.Libdl]]
uuid = "8f399da3-3557-5675-b5ff-fb832c97cbdb"
version = "1.11.0"

[[deps.Libiconv_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "be484f5c92fad0bd8acfef35fe017900b0b73809"
uuid = "94ce4f54-9a6c-5748-9c1c-f9c7231a4531"
version = "1.18.0+0"

[[deps.LinearAlgebra]]
deps = ["Libdl", "OpenBLAS_jll", "libblastrampoline_jll"]
uuid = "37e2e46d-f89d-539d-b4ee-838fcccc9c8e"
version = "1.12.0"

[[deps.Logging]]
uuid = "56ddb016-857b-54e1-b83d-db4d58db5568"
version = "1.11.0"

[[deps.LoggingExtras]]
deps = ["Dates", "Logging"]
git-tree-sha1 = "f00544d95982ea270145636c181ceda21c4e2575"
uuid = "e6f89c97-d47a-5376-807f-9c37f3926c36"
version = "1.2.0"

[[deps.MIMEs]]
git-tree-sha1 = "c64d943587f7187e751162b3b84445bbbd79f691"
uuid = "6c6e2e6c-3030-632d-7369-2d6c69616d65"
version = "1.1.0"

[[deps.MITgcm]]
deps = ["ClimateModels", "DataDeps", "Dataverse", "Dates", "Distributed", "FortranFiles", "Glob", "MeshArrays", "Printf", "Scratch", "SharedArrays", "SparseArrays", "Statistics", "StyledStrings", "UUIDs"]
git-tree-sha1 = "38eaef144504c69593362406b4425d3d2d1a0c9a"
uuid = "dce5fa8e-68ce-4431-a242-9469c69627a0"
version = "0.5.11"

    [deps.MITgcm.extensions]
    MITgcmNetCDFExt = ["NetCDF"]

    [deps.MITgcm.weakdeps]
    NetCDF = "30363a11-5582-574a-97bb-aa9a979735b9"

[[deps.MacroTools]]
git-tree-sha1 = "1e0228a030642014fe5cfe68c2c0a818f9e3f522"
uuid = "1914dd2f-81c6-5fcd-8719-6d5c9610ff09"
version = "0.5.16"

[[deps.Markdown]]
deps = ["Base64", "JuliaSyntaxHighlighting", "StyledStrings"]
uuid = "d6f4376e-aef5-505a-96c1-9c027394607a"
version = "1.11.0"

[[deps.MbedTLS]]
deps = ["Dates", "MbedTLS_jll", "MozillaCACerts_jll", "NetworkOptions", "Random", "Sockets"]
git-tree-sha1 = "c067a280ddc25f196b5e7df3877c6b226d390aaf"
uuid = "739be429-bea8-5141-9913-cc70e7f3736d"
version = "1.1.9"

[[deps.MbedTLS_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "3cce3511ca2c6f87b19c34ffc623417ed2798cbd"
uuid = "c8ffd9c3-330d-5841-b78e-0817d7145fa1"
version = "2.28.10+0"

[[deps.MeshArrays]]
deps = ["CatViews", "Dates", "Distributed", "Glob", "LazyArtifacts", "NearestNeighbors", "Pkg", "Printf", "SharedArrays", "SparseArrays", "Statistics", "Unitful"]
git-tree-sha1 = "fea0859c1406389b3e6657b5376b577ffffdd03a"
uuid = "cb8c808f-1acf-59a3-9d2b-6e38d009f683"
version = "0.3.23"

    [deps.MeshArrays.extensions]
    MeshArraysDataDepsExt = ["DataDeps"]
    MeshArraysGeoJSONExt = ["GeoJSON"]
    MeshArraysJLD2Ext = ["JLD2"]
    MeshArraysMakieExt = ["Makie"]
    MeshArraysProjExt = ["Proj"]
    MeshArraysShapefileExt = ["Shapefile"]

    [deps.MeshArrays.weakdeps]
    DataDeps = "124859b0-ceae-595e-8997-d05f6a7a8dfe"
    GeoJSON = "61d90e0f-e114-555e-ac52-39dfb47a3ef9"
    JLD2 = "033835bb-8acc-5ee8-8aae-3f567f8a3819"
    Makie = "ee78f7c6-11fb-53f2-987a-cfe4a2b5a57a"
    Proj = "c94c279d-25a6-4763-9509-64d165bea63e"
    Shapefile = "8e980c4a-a4fe-5da2-b3a7-4b4b0353a2f4"

[[deps.Missings]]
deps = ["DataAPI"]
git-tree-sha1 = "ec4f7fbeab05d7747bdf98eb74d130a2a2ed298d"
uuid = "e1d29d7a-bbdc-5cf2-9ac0-f12de2c33e28"
version = "1.2.0"

[[deps.Mmap]]
uuid = "a63ad114-7e13-5084-954f-fe012c677804"
version = "1.11.0"

[[deps.MozillaCACerts_jll]]
uuid = "14a3606d-f60d-562e-9121-12d972cd8159"
version = "2025.5.20"

[[deps.NearestNeighbors]]
deps = ["Distances", "StaticArrays"]
git-tree-sha1 = "ca7e18198a166a1f3eb92a3650d53d94ed8ca8a1"
uuid = "b8a86587-4115-5ab1-83bc-aa920d37bbce"
version = "0.4.22"

[[deps.NetworkOptions]]
uuid = "ca575930-c2e3-43a9-ace4-1e988b2c1908"
version = "1.3.0"

[[deps.OffsetArrays]]
git-tree-sha1 = "117432e406b5c023f665fa73dc26e79ec3630151"
uuid = "6fe1bfb0-de20-5000-8ca7-80f57d26f881"
version = "1.17.0"

    [deps.OffsetArrays.extensions]
    OffsetArraysAdaptExt = "Adapt"

    [deps.OffsetArrays.weakdeps]
    Adapt = "79e6a3ab-5dfb-504d-930d-738a2a938a0e"

[[deps.OpenBLASConsistentFPCSR_jll]]
deps = ["Artifacts", "CompilerSupportLibraries_jll", "JLLWrappers", "Libdl"]
git-tree-sha1 = "567515ca155d0020a45b05175449b499c63e7015"
uuid = "6cdc7f73-28fd-5e50-80fb-958a8875b1af"
version = "0.3.29+0"

[[deps.OpenBLAS_jll]]
deps = ["Artifacts", "CompilerSupportLibraries_jll", "Libdl"]
uuid = "4536629a-c528-5b80-bd46-f80d51c5b363"
version = "0.3.29+0"

[[deps.OpenSSH_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "OpenSSL_jll", "Zlib_jll"]
git-tree-sha1 = "301412a644646fdc0ad67d0a87487466b491e53d"
uuid = "9bd350c2-7e96-507f-8002-3f2e150b4e1b"
version = "10.2.1+0"

[[deps.OpenSSL]]
deps = ["BitFlags", "Dates", "MozillaCACerts_jll", "NetworkOptions", "OpenSSL_jll", "Sockets"]
git-tree-sha1 = "386b47442468acfb1add94bf2d85365dea10cbab"
uuid = "4d8831e6-92b7-49fb-bdf8-b643e874388c"
version = "1.6.0"

[[deps.OpenSSL_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "458c3c95-2e84-50aa-8efc-19380b2a3a95"
version = "3.5.1+0"

[[deps.OrderedCollections]]
git-tree-sha1 = "05868e21324cede2207c6f0f466b4bfef6d5e7ee"
uuid = "bac558e1-5e72-5ebc-8fee-abe8a469f55d"
version = "1.8.1"

[[deps.PCRE2_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "efcefdf7-47ab-520b-bdef-62a2eaa19f15"
version = "10.44.0+1"

[[deps.Parsers]]
deps = ["Dates", "PrecompileTools", "UUIDs"]
git-tree-sha1 = "7d2f8f21da5db6a806faf7b9b292296da42b2810"
uuid = "69de0a69-1ddd-5017-9359-2bf0b02dc9f0"
version = "2.8.3"

[[deps.Pkg]]
deps = ["Artifacts", "Dates", "Downloads", "FileWatching", "LibGit2", "Libdl", "Logging", "Markdown", "Printf", "Random", "SHA", "TOML", "Tar", "UUIDs", "p7zip_jll"]
uuid = "44cfe95a-1eb2-52ea-b672-e2afdf69b78f"
version = "1.12.0"
weakdeps = ["REPL"]

    [deps.Pkg.extensions]
    REPLExt = "REPL"

[[deps.PlutoUI]]
deps = ["AbstractPlutoDingetjes", "Base64", "ColorTypes", "Dates", "Downloads", "FixedPointNumbers", "Hyperscript", "HypertextLiteral", "IOCapture", "InteractiveUtils", "JSON", "Logging", "MIMEs", "Markdown", "Random", "Reexport", "URIs", "UUIDs"]
git-tree-sha1 = "3faff84e6f97a7f18e0dd24373daa229fd358db5"
uuid = "7f904dfe-b85e-4ff6-b463-dae2292396a8"
version = "0.7.73"

[[deps.PooledArrays]]
deps = ["DataAPI", "Future"]
git-tree-sha1 = "36d8b4b899628fb92c2749eb488d884a926614d3"
uuid = "2dfb63ee-cc39-5dd5-95bd-886bf059d720"
version = "1.4.3"

[[deps.PrecompileTools]]
deps = ["Preferences"]
git-tree-sha1 = "07a921781cab75691315adc645096ed5e370cb77"
uuid = "aea7be01-6a6a-4083-8856-8a6e6704d82a"
version = "1.3.3"

[[deps.Preferences]]
deps = ["TOML"]
git-tree-sha1 = "0f27480397253da18fe2c12a4ba4eb9eb208bf3d"
uuid = "21216c6a-2e73-6563-6e65-726566657250"
version = "1.5.0"

[[deps.PrettyTables]]
deps = ["Crayons", "LaTeXStrings", "Markdown", "PrecompileTools", "Printf", "REPL", "Reexport", "StringManipulation", "Tables"]
git-tree-sha1 = "6b8e2f0bae3f678811678065c09571c1619da219"
uuid = "08abe8d2-0d0c-5749-adfa-8a2ac140af0d"
version = "3.1.0"

[[deps.Printf]]
deps = ["Unicode"]
uuid = "de0858da-6303-5e67-8744-51eddeeeb8d7"
version = "1.11.0"

[[deps.REPL]]
deps = ["InteractiveUtils", "JuliaSyntaxHighlighting", "Markdown", "Sockets", "StyledStrings", "Unicode"]
uuid = "3fa0cd96-eef1-5676-8a61-b3b8758bbffb"
version = "1.11.0"

[[deps.Random]]
deps = ["SHA"]
uuid = "9a3f8284-a2c9-5f02-9a11-845980a1fd5c"
version = "1.11.0"

[[deps.Reexport]]
git-tree-sha1 = "45e428421666073eab6f2da5c9d310d99bb12f9b"
uuid = "189a3867-3050-52da-a836-e630ba90ab69"
version = "1.2.2"

[[deps.Requires]]
deps = ["UUIDs"]
git-tree-sha1 = "62389eeff14780bfe55195b7204c0d8738436d64"
uuid = "ae029012-a4dd-5104-9daa-d747884805df"
version = "1.3.1"

[[deps.RoundingEmulator]]
git-tree-sha1 = "40b9edad2e5287e05bd413a38f61a8ff55b9557b"
uuid = "5eaf0fd0-dfba-4ccb-bf02-d820a40db705"
version = "0.2.1"

[[deps.SHA]]
uuid = "ea8e919c-243c-51af-8825-aaa63cd721ce"
version = "0.7.0"

[[deps.ScopedValues]]
deps = ["HashArrayMappedTries", "Logging"]
git-tree-sha1 = "c3b2323466378a2ba15bea4b2f73b081e022f473"
uuid = "7e506255-f358-4e82-b7e4-beb19740aa63"
version = "1.5.0"

[[deps.Scratch]]
deps = ["Dates"]
git-tree-sha1 = "9b81b8393e50b7d4e6d0a9f14e192294d3b7c109"
uuid = "6c6a2e73-6563-6170-7368-637461726353"
version = "1.3.0"

[[deps.SentinelArrays]]
deps = ["Dates", "Random"]
git-tree-sha1 = "712fb0231ee6f9120e005ccd56297abbc053e7e0"
uuid = "91c51154-3ec4-41a3-a24f-3f23e20d615c"
version = "1.4.8"

[[deps.Serialization]]
uuid = "9e88b42a-f829-5b0c-bbe9-9e923198166b"
version = "1.11.0"

[[deps.SharedArrays]]
deps = ["Distributed", "Mmap", "Random", "Serialization"]
uuid = "1a1011a3-84de-559e-8e89-a11a2f7dc383"
version = "1.11.0"

[[deps.SimpleBufferStream]]
git-tree-sha1 = "f305871d2f381d21527c770d4788c06c097c9bc1"
uuid = "777ac1f9-54b0-4bf8-805c-2214025038e7"
version = "1.2.0"

[[deps.Sockets]]
uuid = "6462fe0b-24de-5631-8697-dd941f90decc"
version = "1.11.0"

[[deps.SortTileRecursiveTree]]
deps = ["AbstractTrees", "Extents", "GeoInterface"]
git-tree-sha1 = "f9aa6616a9b3bd01f93f27c010f1d25fc5a094a9"
uuid = "746ee33f-1797-42c2-866d-db2fce69d14d"
version = "0.1.4"

[[deps.SortingAlgorithms]]
deps = ["DataStructures"]
git-tree-sha1 = "64d974c2e6fdf07f8155b5b2ca2ffa9069b608d9"
uuid = "a2af1166-a08f-5f64-846c-94a0d3cef48c"
version = "1.2.2"

[[deps.SparseArrays]]
deps = ["Libdl", "LinearAlgebra", "Random", "Serialization", "SuiteSparse_jll"]
uuid = "2f01184e-e22b-5df5-ae63-d93ebab69eaf"
version = "1.12.0"

[[deps.StableTasks]]
git-tree-sha1 = "c4f6610f85cb965bee5bfafa64cbeeda55a4e0b2"
uuid = "91464d47-22a1-43fe-8b7f-2d57ee82463f"
version = "0.1.7"

[[deps.StaticArrays]]
deps = ["LinearAlgebra", "PrecompileTools", "Random", "StaticArraysCore"]
git-tree-sha1 = "b8693004b385c842357406e3af647701fe783f98"
uuid = "90137ffa-7385-5640-81b9-e52037218182"
version = "1.9.15"

    [deps.StaticArrays.extensions]
    StaticArraysChainRulesCoreExt = "ChainRulesCore"
    StaticArraysStatisticsExt = "Statistics"

    [deps.StaticArrays.weakdeps]
    ChainRulesCore = "d360d2e6-b24c-11e9-a2a3-2a2ae2dbcce4"
    Statistics = "10745b16-79ce-11e8-11f9-7d13ad32a3b2"

[[deps.StaticArraysCore]]
git-tree-sha1 = "6ab403037779dae8c514bad259f32a447262455a"
uuid = "1e83bf80-4336-4d27-bf5d-d5a4f845583c"
version = "1.4.4"

[[deps.Statistics]]
deps = ["LinearAlgebra"]
git-tree-sha1 = "ae3bb1eb3bba077cd276bc5cfc337cc65c3075c0"
uuid = "10745b16-79ce-11e8-11f9-7d13ad32a3b2"
version = "1.11.1"
weakdeps = ["SparseArrays"]

    [deps.Statistics.extensions]
    SparseArraysExt = ["SparseArrays"]

[[deps.StatsAPI]]
deps = ["LinearAlgebra"]
git-tree-sha1 = "9d72a13a3f4dd3795a195ac5a44d7d6ff5f552ff"
uuid = "82ae8749-77ed-4fe6-ae5f-f523153014b0"
version = "1.7.1"

[[deps.StringManipulation]]
deps = ["PrecompileTools"]
git-tree-sha1 = "725421ae8e530ec29bcbdddbe91ff8053421d023"
uuid = "892a3eda-7b42-436c-8928-eab12a02cf0e"
version = "0.4.1"

[[deps.StyledStrings]]
uuid = "f489334b-da3d-4c2e-b8f0-e476e12c162b"
version = "1.11.0"

[[deps.SuiteSparse_jll]]
deps = ["Artifacts", "Libdl", "libblastrampoline_jll"]
uuid = "bea87d4a-7f5b-5778-9afe-8cc45184846c"
version = "7.8.3+2"

[[deps.Suppressor]]
deps = ["Logging"]
git-tree-sha1 = "6dbb5b635c5437c68c28c2ac9e39b87138f37c0a"
uuid = "fd094767-a336-5f1f-9728-57cf17d0bbfb"
version = "0.2.8"

[[deps.TOML]]
deps = ["Dates"]
uuid = "fa267f1f-6049-4f14-aa54-33bafae1ed76"
version = "1.0.3"

[[deps.TableTraits]]
deps = ["IteratorInterfaceExtensions"]
git-tree-sha1 = "c06b2f539df1c6efa794486abfb6ed2022561a39"
uuid = "3783bdb8-4a98-5b6b-af9a-565f29a5fe9c"
version = "1.0.1"

[[deps.Tables]]
deps = ["DataAPI", "DataValueInterfaces", "IteratorInterfaceExtensions", "OrderedCollections", "TableTraits"]
git-tree-sha1 = "f2c1efbc8f3a609aadf318094f8fc5204bdaf344"
uuid = "bd369af6-aec1-5ad0-b16a-f7cc5008161c"
version = "1.12.1"

[[deps.Tar]]
deps = ["ArgTools", "SHA"]
uuid = "a4e569a6-e804-4fa4-b0f3-eef7a1d5b13e"
version = "1.10.0"

[[deps.Test]]
deps = ["InteractiveUtils", "Logging", "Random", "Serialization"]
uuid = "8dfed614-e22c-5e08-85e1-65c5234f0b40"
version = "1.11.0"

[[deps.TranscodingStreams]]
git-tree-sha1 = "0c45878dcfdcfa8480052b6ab162cdd138781742"
uuid = "3bb67fe8-82b1-5028-8e26-92a6c54297fa"
version = "0.11.3"

[[deps.Tricks]]
git-tree-sha1 = "311349fd1c93a31f783f977a71e8b062a57d4101"
uuid = "410a4b4d-49e4-4fbc-ab6d-cb71b17b3775"
version = "0.1.13"

[[deps.URIs]]
git-tree-sha1 = "bef26fb046d031353ef97a82e3fdb6afe7f21b1a"
uuid = "5c2747f8-b7ea-4ff2-ba2e-563bfd36b1d4"
version = "1.6.1"

[[deps.UUIDs]]
deps = ["Random", "SHA"]
uuid = "cf7118a7-6976-5b1a-9a39-7adc72f591a4"
version = "1.11.0"

[[deps.Unicode]]
uuid = "4ec0a83e-493e-50e2-b9ac-8f72acf5a8f5"
version = "1.11.0"

[[deps.Unitful]]
deps = ["Dates", "LinearAlgebra", "Random"]
git-tree-sha1 = "83360bda12f61c250835830cc40b64f487cc2230"
uuid = "1986cc42-f94f-5a68-af5c-568840ba703d"
version = "1.25.1"

    [deps.Unitful.extensions]
    ConstructionBaseUnitfulExt = "ConstructionBase"
    ForwardDiffExt = "ForwardDiff"
    InverseFunctionsUnitfulExt = "InverseFunctions"
    LatexifyExt = ["Latexify", "LaTeXStrings"]
    PrintfExt = "Printf"

    [deps.Unitful.weakdeps]
    ConstructionBase = "187b0558-2788-49d3-abe0-74a17ed4e7c9"
    ForwardDiff = "f6369f11-7733-5829-9624-2563aa707210"
    InverseFunctions = "3587e190-3f89-42d0-90ee-14403ec27112"
    LaTeXStrings = "b964fa9f-0449-5b57-a5c2-d3ea65f4040f"
    Latexify = "23fbe1c1-3f47-55db-b15f-69d7ec21a316"
    Printf = "de0858da-6303-5e67-8744-51eddeeeb8d7"

[[deps.WeakRefStrings]]
deps = ["DataAPI", "InlineStrings", "Parsers"]
git-tree-sha1 = "b1be2855ed9ed8eac54e5caff2afcdb442d52c23"
uuid = "ea10d353-3f73-51f8-a26c-33c1cb351aa5"
version = "1.4.2"

[[deps.WorkerUtilities]]
git-tree-sha1 = "cd1659ba0d57b71a464a29e64dbc67cfe83d54e7"
uuid = "76eceee3-57b5-4d4a-8e66-0e911cebbf60"
version = "1.6.1"

[[deps.ZipFile]]
deps = ["Libdl", "Printf", "Zlib_jll"]
git-tree-sha1 = "f492b7fe1698e623024e873244f10d89c95c340a"
uuid = "a5390f91-8eb1-5f08-bee0-b1d1ffed6cea"
version = "0.10.1"

[[deps.Zlib_jll]]
deps = ["Libdl"]
uuid = "83775a58-1f1d-513f-b197-d71354ab007a"
version = "1.3.1+2"

[[deps.Zstd_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "446b23e73536f84e8037f5dce465e92275f6a308"
uuid = "3161d3a3-bdf6-5164-811a-617609db77b4"
version = "1.5.7+1"

[[deps.libblastrampoline_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "8e850b90-86db-534c-a0d3-1478176c7d93"
version = "5.15.0+0"

[[deps.nghttp2_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "8e850ede-7688-5339-a07c-302acd2aaf8d"
version = "1.64.0+1"

[[deps.p7zip_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "3f19e933-33d8-53b3-aaab-bd5110c3b7a0"
version = "17.5.0+2"
"""

# ╔═╡ Cell order:
# ╟─6da3809b-4bf4-45d6-bbd2-382312e947b7
# ╟─7aca186c-efd7-4312-afd2-37fdc1a2bbb6
# ╟─3e5aac32-c82a-45e2-8d7b-05cf96a4ce3c
# ╠═1b5b7cdc-9528-4afc-b282-2828da34ae2f
# ╟─7bace9d8-4e9e-451d-aef4-5093c62b08e7
# ╠═17dd67f4-0c98-408a-9cf9-b02e9a00b5a3
# ╠═2fd20020-8b29-4262-92f3-aabf9bc609fa
# ╟─53b62449-7006-4f2a-a337-0b0e45d82f5e
# ╟─65dcb624-7f27-4247-bc73-446bf069a67a
# ╟─7d0afd5e-fa88-4f40-a1f0-1a2a944341cd
# ╟─a248a2b1-ed63-48a0-922c-b00a3dbfc294
# ╟─87bb371f-c08e-4a31-b7bf-56b58fe9cbda
# ╠═f112ecfa-8a66-4ae1-aac3-7b354bb5863f
# ╟─915a5f42-9a31-4837-a891-d22f8825930b
# ╟─3c2f45b3-b3c3-4704-ba0e-c1ada7c48c13
# ╟─8f3ffcf6-9aa8-4882-b055-583ff2bd82f4
# ╠═a5fcb7cb-a86c-4f19-a3b9-1e280ce6159b
# ╠═ac64b4cf-a6b5-4dad-a90c-bf1f483690e5
# ╟─523311f3-8d0a-4b57-b376-da945c4f20a7
# ╟─4208e083-2971-4b1c-bc72-5df1caaf3d36
# ╟─199f9dfc-c393-47f9-95d8-c0a783ae47e6
# ╠═540641e2-14e7-47da-9f8e-ac37cd499886
# ╟─3549d2f2-ab0f-44a0-a0dd-b9aa64354f31
# ╟─798f41ed-cd77-4890-9148-be278df52ba4
# ╟─6a354d1f-eba1-40f6-bdfc-1ecd8d96a0ef
# ╟─1effaab8-8237-4236-ad75-755795f898fc
# ╟─98496427-b774-4d9b-b386-ffe5ea775d37
# ╟─fc2e3116-8a7e-405d-ad6d-bdf5d4ae45c9
# ╟─f4d0e938-e544-4add-be50-7942df72a3b8
# ╟─5f5e6be3-a7ad-4772-a52c-acb765106593
# ╟─a4418a2e-831b-416f-8348-e282b08cfa2d
# ╟─5365dc71-0ddc-4d40-8f38-7c494b3305d9
# ╟─b494d021-3727-4736-8031-9587d0d542a2
# ╟─e4d61551-c0e7-4a61-b38e-774897862fed
# ╟─50272674-6935-4239-89bf-29997341148b
# ╟─63c58499-c353-4124-9270-6c22b0e6bcd7
# ╟─9405eb97-1130-427d-a9f4-e45efea361b2
# ╟─00000000-0000-0000-0000-000000000001
# ╟─00000000-0000-0000-0000-000000000002
