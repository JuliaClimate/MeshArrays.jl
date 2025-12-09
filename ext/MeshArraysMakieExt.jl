
module MeshArraysMakieExt

using MeshArrays, Makie

import MeshArrays: GI, polyarray, gridpath
import MeshArrays: plot_examples, ProjAxis, grid_lines!

import Makie: plot, plot!, scatter, scatter!, surface!
import Makie: lines, lines!, heatmap, heatmap!, contour!, contourf!

LineString=Makie.LineString
Observable=Makie.Observable
GeometryBasics=Makie.GeometryBasics

##

include("MeshArraysMakieExt_recipes.jl")
include("MeshArraysMakieExt_PrAxis.jl")
include("MeshArraysMakieExt_split.jl")
include("MeshArraysMakieExt_examples.jl")

end
