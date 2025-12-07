
"""
struct NamedPolygon
    geometry::Vector{Tuple{Float32, Float32}}
    name::String
    points::Vector{}
end

Use example:

```
pol=MeshArrays.Dataset("oceans_geojson1")
g=pol[1].geometry
p=MeshArrays.GI.coordinates(g)[1]
n=pol[1].name
npol1=MeshArrays.NamedPolygon(p,n,p)
apol1=MeshArrays.polyarray("anonymous",[npol1 npol1])

using CairoMakie
MeshArraysMakieExt = Base.get_extension(MeshArrays, :MeshArraysMakieExt)
lines(MeshArraysMakieExt.pol_to_Makie(apol1))
```
"""
struct NamedPolygon
    geometry::Vector{Union{Tuple,Vector}}
    name::String
    points::Vector{}
end

"""
    polyarray

polyarray data structure.

```
struct polyarray
   name::String
   f::Array{NamedPolygon}
end
```

"""
struct polyarray
   name::String
   f::Array{NamedPolygon}
end

"""

```
using DataFrames
function write_polygons_to_json(pols,nams)
        np=length(pols)
        df = DataFrame(geometry=pols)
        df[!,:id] = nams
        df[!, :name] = nams

        fn = tempname()*".json"
        GeoJSON.write(fn, df)
end

pols,nams=MeshArrays.to_pol(apol1)
write_polygons_to_json(pols,nams)
```
"""
function to_pol(apol1)
    nk=length(apol1.f)
	pols=Array{GI.Wrappers.Polygon}(undef, nk)
    nams=Array{String}(undef, nk)
	for kk in 1:nk
        line=GI.LineString(apol1.f[kk].geometry)
		pols[kk]=GI.Polygon(line)
        name=apol1.f[kk].name
		nams[kk]=apol1.f[kk].name
	end
    pols,nams
end

