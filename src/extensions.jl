
function read_jld2 end
function write_jld2 end
function ProjAxis end
function grid_lines! end

"""
    read_shp(fil; format=:polyarray)

Call `Shapefile.Table` and `Shapefile.shapes`
and return `polyarray` (default).

- format `:polyarray` (default) : return `polyarray`
- format `:Shapefile` : return vector of `name,geometry`` named tuples.
- format `:coord` : convert to `GI.coordinates`.

```
using MeshArrays, DataDeps, Shapefile
fil=MeshArrays.Dataset("countries_shp1")

using CairoMakie
lines(pol)
```
"""
function read_shp end

"""
    read_json(fil; format=:polyarray)

Call `GeoJSON.read` and return `polyarray` (default)

- format `:polyarray`` (default): return `polyarray`
- format `:coord` : convert to `GI.coordinates`.
- format `:GeoJSON` : return FeatureCollection (vector of name,geom pairs)

```
import MeshArrays, DataDeps, GeoJSON
pol=MeshArrays.Dataset("oceans_geojson1")
```
"""
function read_json end

"""
    within_pol(pol; ID=1)

Generate a `name,rule` pair to test if location 
    `lon,lat` is within polygon `pol[ID].geometry`.

```
using MeshArrays, DataDeps, GeoJSON, GeometryOps
fil=MeshArrays.Dataset("oceans_geojson1",do_read=false)
pol=MeshArrays.read_json(fil,format=:GeoJSON)

name,rule=MeshArrays.within_pol(pol; ID=11)
rule(-30,40)
```
"""
function within_pol end

