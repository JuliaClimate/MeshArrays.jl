
import Plots: heatmap, contour, contourf

"""
    heatmap(fld::MeshArray; args...)

Calls `heatmap` for one elementary array per subplot
"""
function heatmap(x::MeshArray; args...)
    n=length(x.fIndex); p=()
    for i=1:n; p=(p...,heatmap(x.f[i]; args...)); end;
    plot(p...)
end

"""
    contour(fld::MeshArray; args...)

Calls `contour` for one elementary array per subplot
"""
function contour(x::MeshArray; args...)
    n=length(x.fIndex); p=()
    for i=1:n; p=(p...,contour(x.f[i]; args...)); end;
    plot(p...)
end
