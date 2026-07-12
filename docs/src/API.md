# API documentation

## 1. Data Structures 

By default, the `MeshArray` type is an alias to [`MeshArrays.gcmarray`](@ref).

```@docs
AbstractMeshArray
MeshArrays.gcmarray
gcmgrid
varmeta
gridpath
gridmask
```

### Internals

```@docs
MeshArrays.InnerArray
MeshArrays.OuterArray
```

### More

```@docs
MeshArrays.gcmvector
MeshArrays.gcmfaces
```

## 2. Grids And I/O

Grid definitions :

```@docs
GridSpec
MeshArrays.GridSpec_default
MeshArrays.GridSpec_MITgcm
MeshArrays.GridSpec_ones
```

Read / write methods :

```@docs
MeshArrays.read
MeshArrays.read!
MeshArrays.write
```

Loading grids :

```@docs
GridLoad
MeshArrays.GridLoad_default
MeshArrays.GridLoad_ones
GridLoadVar
```

### Functionalities

```@docs
Tiles
Tiles!
exchange
```

## 3. Interpolation

```@docs
interpolation_setup
Interpolate
InterpolationFactors
knn
StereographicProjection
```

## 4. Vector Fields

```@docs
curl
convergence
gradient
ScalarPotential
VectorPotential
ThroughFlow
UVtoTransport
UVtoUEVN
MeshArrays.EkmanTrsp
MeshArrays.calc_bolus
```

### Integration Paths

```@docs
LatitudeCircles
MeshArrays.LatitudeCircle
Transect
edge_path
MeshArrays.edge_mask
```

## 5. Integration Loop

```@docs
Integration.loops
Integration.streamlined_loop
Integration.define_regions
Integration.define_sums
```

## 6. Grid Data Sets

```@docs
MeshArrays.Dataset
```

### Internals

```@docs
MeshArrays.mydatadep
```

## 7. Polygons

```@docs
MeshArrays.within_pol
MeshArrays.read_json
MeshArrays.read_shp
MeshArrays.NamedPolygon
MeshArrays.polyarray
MeshArrays.to_Polygon
MeshArrays.to_polyarray
```

## 8. Other

```@docs
isosurface
land_mask
demo.ocean_basins
demo.extended_basin
```

### Masked Array Utilities

```@docs
nansum
nanmean
nanmax
nanmin
```
