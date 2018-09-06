
using Plots, Dates; gr();

export qwckplot

##  qwckplot function
#
#examples:
#  GCMGridLoad()
#  qwckplot(MeshArrays.Depth)
#  qwckplot(MeshArrays.Depth,"Ocean Bottom Depth")

function qwckplot(fld::gcmfaces)
    tmp1=Dates.now()
    tmp1=Dates.format(tmp1, "yyyy-mm-dd HH:MM:SS")
    qwckplot(fld,"Plotted at time "*tmp1)
end

function qwckplot(fld::gcmfaces,ttl::String)
    arr=convert2array(fld)
    arr=permutedims(arr,[2 1])
    #This uses Plots.jl:
    p=heatmap(arr,title=ttl)
end
