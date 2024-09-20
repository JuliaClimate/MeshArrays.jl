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

More : 

```@docs
MeshArrays.gcmvector
MeshArrays.gcmfaces
```

## 2. Grids And I/O

```@docs
simple_periodic_domain
GridSpec
GridLoad
GridLoadVar
```

More : 

```@docs
Tiles
Tiles!
exchange
MeshArrays.read
MeshArrays.read!
MeshArrays.write
UnitGrid
```

## 3. Interpolation

```@docs
interpolation_setup
Interpolate
InterpolationFactors
StereographicProjection
knn
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
```

## 5. Integration

```@docs
Integration.loops
```

## 6. Other

```@docs
LatitudeCircles
Transect
demo.ocean_basins
demo.extended_basin
isosurface
MA_datadep
```
