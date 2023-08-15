# API documentation

## 1. Data Structures 

By default, the `MeshArray` type is an alias to [`MeshArrays.gcmarray`](@ref).

```@docs
AbstractMeshArray
MeshArrays.gcmarray
gcmgrid
varmeta
```

More : 

```@docs
MeshArrays.gcmvector
MeshArrays.gcmfaces
```

## 2. Grids And I/O

```@docs
UnitGrid
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
```

## 3. Interpolation

```@docs
knn
Interpolate
InterpolationFactors
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
```

## 5. Other

```@docs
LatitudeCircles
Transect
isosurface
```
