
import Plots: heatmap, contour, contourf

"""
    heatmap(fld::MeshArray; args...)

Calls `heatmap` for one elementary array per subplot
"""
function heatmap(x::MeshArray; args...)
    n=x.grid.nFaces
    isa(x,MeshArrays.gcmarray) ? n=length(x.fIndex) : nothing

    m=Int(sqrt(n))
    j=collect(1:n)
    j=reshape(collect(j),(m,m))
    j=reverse(permutedims(j,(2,1)); dims=2)[:]

    p=()
    for i=1:n; p=(p...,heatmap(x.f[j[i]]; args...)); end
    plot(p...)
end

"""
    contour(fld::MeshArray; args...)

Calls `contour` for one elementary array per subplot
"""
function contour(x::MeshArray; args...)
    n=x.grid.nFaces
    isa(x,MeshArrays.gcmarray) ? n=length(x.fIndex) : nothing
    p=()
    for i=1:n; p=(p...,contour(x.f[i]; args...)); end
    plot(p...)
end
