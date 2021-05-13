# API documentation

## 1. Data Structures 

```@docs
AbstractMeshArray
MeshArrays.gcmarray
MeshArrays.gcmvector
MeshArrays.gcmfaces
gcmgrid
varmeta
```

## 2. Grids And I/O

```@docs
GridSpec
GridOfOnes
simple_periodic_domain
GridLoad
MeshArrays.read
MeshArrays.write
exchange
Tiles
```

## 3. Interapolation, etc

```@docs
knn
Interpolate
InterpolationFactors
StereographicProjection
```

## 4. Vector Fields

```@docs
gradient
convergence
smooth
ScalarPotential
VectorPotential
LatitudeCircles
ThroughFlow
```

## 5. Various

```@docs
MeshArrays.getindexetc
MeshArrays.nFacesEtc
findall
mask
```
