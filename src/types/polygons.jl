
"""
struct NamedPolygon
    geometry::Vector{Tuple{Float32, Float32}}
    name::String
    points::Vector{}
end


```
pol=MeshArrays.Dataset("oceans_geojson1")
p=pol[1].geometry.coordinates[1]
n=pol[1].name
npol1=MeshArrays.NamedPolygon(p,n,p)
apol1=MeshArrays.polyarray("anonymous",[npol1 npol1])
```
"""
struct NamedPolygon
    geometry::Vector{Tuple{Float32, Float32}}
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
