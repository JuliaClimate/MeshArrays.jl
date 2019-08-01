
using Plots, Dates; gr();

export qwckplot

##  qwckplot function

"""
    qwckplot(fld::gcmfaces)

Call qwckplot(fld::gcmfaces,ttl::String) with current date for title. Example:

```
!isdir("GRID_LLC90") ? error("missing files") : nothing
mygrid=GCMGridLoad(GCMGridSpec("LLC90"))
qwckplot(mygrid["Depth"])
```
"""
function qwckplot(fld::gcmfaces)
    tmp1=Dates.now()
    tmp1=Dates.format(tmp1, "yyyy-mm-dd HH:MM:SS")
    qwckplot(fld,"Plotted at time "*tmp1)
end

"""
    qwckplot(fld::gcmfaces,ttl::String)

Plot input using convert2array and heatmap + add title
"""
function qwckplot(fld::gcmfaces,ttl::String)
    arr=convert2array(fld)
    arr=permutedims(arr,[2 1])
    #This uses Plots.jl:
    p=heatmap(arr,title=ttl)
end
