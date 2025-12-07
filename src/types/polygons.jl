
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
pol_a=MeshArrays.to_polyarray(pol)

using CairoMakie
MeshArraysMakieExt = Base.get_extension(MeshArrays, :MeshArraysMakieExt)
lines(MeshArraysMakieExt.pol_to_Makie(pol_a))
```
"""
struct polyarray
   name::String
   data::Array{NamedPolygon}
end

##

function to_polyarray(pol)
    nk=length(pol)
    npol=Array{NamedPolygon}(undef, nk)
    for k in 1:nk
        g=pol[k].geometry
        p=MeshArrays.GI.coordinates(g)[1]
        n=pol[k].name
        npol[k]=MeshArrays.NamedPolygon(p,n,p)
    end
    MeshArrays.polyarray("anonymous",npol)
end

##

"""
    to_Polygon(pa::polyarray)

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

pols,nams=MeshArrays.to_Polygon(pol_a)
write_polygons_to_json(pols,nams)
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

