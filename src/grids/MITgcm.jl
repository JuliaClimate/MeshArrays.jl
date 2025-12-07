

## GridSpec function with category argument:

"""
    GridSpec_MITgcm(category="PeriodicDomain",path=tempname(); np=nothing, ID=:unknown)

See the main `GridSpec` documentation.
"""
function GridSpec_MITgcm(category="PeriodicDomain", 
    path=tempname(); np=nothing, ID=:unknown)

if category=="LatLonCap"
    nFaces=5
    grTopo="LatLonCap"
    np === nothing ? np=90 : np
    ioSize=[np np*13]
    facesSize=[(np, np*3), (np, np*3), (np, np), (np*3, np), (np*3, np)]
    if np==270
        ioPrec=Float32
    else
        ioPrec=Float64
    end
elseif category=="CubeSphere"
    nFaces=6
    grTopo="CubeSphere"
    np === nothing ? np=32 : np
    ioSize=[np np*nFaces]
    facesSize=[(np, np), (np, np), (np, np), (np, np), (np, np), (np, np)]
    ioPrec=Float32
elseif category=="PeriodicChannel"
    nFaces=1
    grTopo="PeriodicChannel"
    ioSize=[360 160]
    facesSize=[(360, 160)]
    ioPrec=Float32
elseif category=="PeriodicDomain"
    nFaces=4
    grTopo="PeriodicDomain"
    ioSize=[80 42]
    facesSize=[(40, 21), (40, 21), (40, 21), (40, 21)]
    ioPrec=Float32
elseif ID==:unknown
    error("unknown category case")
end

if ID==:unknown
    gcmgrid(path, grTopo, nFaces, facesSize, ioSize, ioPrec, read, write)
elseif ID==:LLC90
    np = 90
    GridSpec("LatLonCap", MeshArrays.Dataset("GRID_LLC90"), np=np)
elseif ID==:LLC270
    np = 270
    GridSpec("LatLonCap", MeshArrays.Dataset("GRID_LLC270"), np=np)
elseif ID==:CS32
    np = 32
    GridSpec("CubeSphere", MeshArrays.Dataset("GRID_CS32"), np=np)
elseif ID==:onedegree
    GridSpec("PeriodicChannel", MeshArrays.Dataset("GRID_LL360"))
elseif ID==:PeriodicDomain
    GridSpec_MITgcm("PeriodicDomain")
else
    error("unknwown grid")
end

end


