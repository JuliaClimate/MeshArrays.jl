
"""
struct NamedPolygon
    geometry::Vector{Tuple{Float32, Float32}}
    name::String
    points::Vector{}
end
"""
struct NamedPolygon
    geometry::Vector{Union{Tuple,Vector}}
    name::String
    points::Vector{}
end

"""
struct polyarray
   name::String
   f::Array{NamedPolygon}
end

```
using MeshArrays, GeoJSON, DataDeps
pol=MeshArrays.Dataset("oceans_geojson1")

write(pol,tempname()*".json")

using CairoMakie
lines(pol)
```
"""
struct polyarray
   name::String
   data::Array{NamedPolygon}
end

##

to_polyarray(pol::polyarray) = pol

"""
     to_polyarray(pol)

Convert polygon data into a `polyarray`.

```
fil="countries.geojson"
MeshArrays.to_polyarray(GeoJSON.read(fil))
```
"""
function to_polyarray(pol)
    nk=length(pol)
    npol=NamedPolygon[]
    for k in 1:nk
        g=pol[k].geometry
        for i in MeshArrays.GI.coordinates(g)
            n=pol[k].name
            push!(npol,MeshArrays.NamedPolygon(i,n,i))
        end
    end
    MeshArrays.polyarray("anonymous",npol)
end

##

"""
    to_Polygon(pa::polyarray)

Convert `polyarray` into a `GI.Polygon` array.

```
pol_P,nams=MeshArrays.to_Polygon(pol)
```
"""
function to_Polygon(pa::polyarray)
    nk=length(pa.data)
	pols=Array{GI.Wrappers.Polygon}(undef, nk)
    nams=Array{String}(undef, nk)
	for kk in 1:nk
        line=GI.LineString(pa.data[kk].geometry)
		pols[kk]=GI.Polygon(line)
        name=pa.data[kk].name
		nams[kk]=pa.data[kk].name
	end
    pols,nams
end

