
using Plots, Dates

export qwckplot

##  qwckplot function

"""
    qwckplot(fld::MeshArray)

Call qwckplot(fld::MeshArray) with date as title. Example:

```
!isdir("GRID_LLC90") ? error("missing files") : nothing
GridVariables=GridLoad(GridSpec("LLC90"))
qwckplot(GridVariables["Depth"])
```
"""
function qwckplot(fld::MeshArray)
    tmp1=Dates.now()
    tmp1=Dates.format(tmp1, "yyyy-mm-dd HH:MM:SS")
    qwckplot(fld,"Plotted at time "*tmp1)
end

"""
    qwckplot(fld::MeshArray,ttl::String)

Plot input using convert2array and heatmap + add title
"""
function qwckplot(fld::MeshArray,ttl::String)
    arr=convert2array(fld)
    arr=permutedims(arr,[2 1])
    #This uses Plots.jl:
    p=heatmap(arr,title=ttl)
end

"""
    qwckplot(fld::MeshArray,clims::NTuple{2, Number})

Plot input using convert2array and heatmap w. chosen clims
"""
function qwckplot(fld::MeshArray,clims::NTuple{2, Number})
    arr=convert2array(fld)
    arr=permutedims(arr,[2 1])
    #This uses Plots.jl:
    p=heatmap(arr,clims=clims)
end
