
using Plots, Dates; gr();

export qwckplot

##  qwckplot function

"""
    qwckplot(fld::AbstractMeshArray)

Call qwckplot(fld::AbstractMeshArray,ttl::String) with current date for title. Example:

```
!isdir("GRID_LLC90") ? error("missing files") : nothing
GridVariables=GCMGridLoad(GCMGridSpec("LLC90"))
qwckplot(GridVariables["Depth"])
```
"""
function qwckplot(fld::AbstractMeshArray)
    tmp1=Dates.now()
    tmp1=Dates.format(tmp1, "yyyy-mm-dd HH:MM:SS")
    qwckplot(fld,"Plotted at time "*tmp1)
end

"""
    qwckplot(fld::AbstractMeshArray,ttl::String)

Plot input using convert2array and heatmap + add title
"""
function qwckplot(fld::AbstractMeshArray,ttl::String)
    arr=MeshArrays.convert2array(fld)
    arr=permutedims(arr,[2 1])
    #This uses Plots.jl:
    p=heatmap(arr,title=ttl)
end
