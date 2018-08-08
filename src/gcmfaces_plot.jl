# gcmfaces_plot.jl
#
#	First Draft Implementation
#
# gaelforget (https://github.com/gaelforget/gcmfaces_jl)
# Julia 0.6.2
# Created: 02.01.18
# Last Edit: 02.01.18

##  Quick Plot In Array Format  ##
function qwckplot(fld)
    arr=convert2array(fld);
    arr=permutedims(arr,[2 1]);
    #This uses Plots.jl:
    p=heatmap(arr)
end
