## Module associated with the MeshArrays package

module MeshArrays

using Pkg, Pkg.Artifacts
thistoml=joinpath(dirname(pathof(MeshArrays)), "..", "Project.toml")
thisversion=Pkg.TOML.parsefile(thistoml)["version"]

p=dirname(pathof(MeshArrays))
artifact_toml = joinpath(p, "../Artifacts.toml")

GRID_LLC90_hash = artifact_hash("GRID_LLC90", artifact_toml)
GRID_LLC90 = joinpath(artifact_path(GRID_LLC90_hash)*"/","GRID_LLC90-1.1/")
GRID_LL360_hash = artifact_hash("GRID_LL360", artifact_toml)
GRID_LL360 = joinpath(artifact_path(GRID_LL360_hash)*"/","GRID_LL360-1.0/")
GRID_CS32_hash = artifact_hash("GRID_CS32", artifact_toml)
GRID_CS32 = joinpath(artifact_path(GRID_CS32_hash)*"/","GRID_CS32-1.1/")

include("Types.jl");
include("Grids.jl");
include("Operations.jl");
include("Exchanges.jl");
include("ReadWrite.jl");
include("Solvers.jl");
include("Interpolation.jl");

export AbstractMeshArray, MeshArray, InnerArray, OuterArray
export gcmgrid, exchange, gradient, convergence, smooth, mask
export simple_periodic_domain, GridSpec, GridLoad, GridOfOnes, GridAddWS!
export Tiles, Interpolate, InterpolationFactors, knn
export ScalarPotential, VectorPotential, ThroughFlow
export StereographicProjection, LatitudeCircles
export maximum, minimum, mean, std

#The following exch_UV differs from normal exchange; incl. exch_UV_N
export exch_UV

#ggebbie: think about using AbstractMeshArray as input type.
"""
    function maximum
    Author ggebbie
    Compute maximum value of gcmgrid type 
# Input
- `x::MeshArrays.gcmarray{Float32,1,Array{Float32,2}}`: input of gcmarray type
- `dryval::Float32`: land value to be eliminated in calculation
# Output
- `xmax::`: maximum value of 2D field
"""
function maximum(x::gcmarray{Float32,1,Array{Float32,2}},dryval::Float32)

    isdry(z) = (z == dryval)
    
    #  vector list of nonzero elements
    xcount = [sum(count(!isdry,x[i])) for i in eachindex(x)]
    if sum(xcount) > 0
        xmax = maximum([maximum(filter(!isdry,x[i])) for i in eachindex(x) if xcount[i] > 0])
    else
        xmax = NaN32
    end
    return xmax
end

"""
    function minimum
    Author ggebbie
    Compute minimum value of gcmgrid type 
# Input
- `x::MeshArrays.gcmarray{Float32,1,Array{Float32,2}}`: input of gcmarray type
- `dryval::Float32`: land value to be eliminated in calculation
# Output
- `xmin::Float32`: minimum value of 2D field
"""
function minimum(x::gcmarray{Float32,1,Array{Float32,2}},dryval::Float32)
    xmin = -maximum(-x,dryval)
    return xmin
end

"""
    function mean
    Author: ggebbie
    Mean of gcmgrid type 
# Input
- `x::MeshArrays.gcmarray{Float32,1,Array{Float32,2}}`: input of gcmarray type
- `dryval::Float32`: land value (doesn't work for NaN32)
# Output
- `xbar::Float32`: mean value (unweighted)
"""
function mean(x::gcmarray{Float32,1,Array{Float32,2}},dryval::Float32)

    isdry(z) = (z == dryval)
    
    #  vector list of nonzero elements
    xcount = [sum(count(!isdry,x[i])) for i in eachindex(x)]
    if sum(xcount) > 0
        # don't assume 0 on land
        xsum = sum([sum(filter(!isdry,x[i])) for i in eachindex(x) if xcount[i] > 0])
        xbar = xsum/sum(xcount)
    else
        xbar = NaN32
    end
    return xbar
end

"""
    function std
    Author: ggebbie
    Compute standard deviation of gcmgrid type 
# Input
- `x::MeshArrays.gcmarray{Float32,1,Array{Float32,2}}`: input of gcmarray type
- `xbar::Float32`: mean value
- `dryval::Float32`: land value (doesn't work for NaN32)
# Output
- `σx::Float32`: standard deviation 
"""
function std(x::gcmarray{Float32,1,Array{Float32,2}},xbar::Float32,dryval)

    isdry(z) = (z == dryval)

    # x prime = fluctuation
    x′ = x .- xbar
    
    #  vector list of nonzero elements
    xcount = [sum(count(!isdry,x′[i])) for i in eachindex(x′)]
    if sum(xcount) > 0
        # don't assume 0 on land
        x²sum = sum([sum(filter(!isdry,x′[i]).^2) for i in eachindex(x) if xcount[i] > 0])
        σx = sqrt(x²sum/(sum(xcount)-1))
    else
        xbar = NaN32
    end
    return σx
end

end # module
