## Module associated with the MeshArrays package

module MeshArrays

using LazyArtifacts
p=dirname(pathof(MeshArrays))
artifact_toml = joinpath(p, "../Artifacts.toml")

function read_polygons end
function read_shp end
function read_json end
function read_JLD2 end
function write_JLD2 end
function ProjAxis end
function grid_lines! end

include("Types.jl")
include("Grids.jl")
include("GridPaths.jl")
include("Operations.jl")
include("Exchanges.jl")
include("ReadWrite.jl")
include("Solvers.jl")
include("ReIndexing.jl")
include("Interpolation.jl")
include("VerticalDimension.jl")
include("Grids_simple.jl")
include("demo.jl")
include("Integration.jl")

export AbstractMeshArray, MeshArray, MeshArray_wh
export gcmgrid, varmeta, gridpath, gridmask
export GridSpec, GridLoad, GridLoadVar, Grids_simple
export exchange, Tiles, Tiles!, Interpolate, InterpolationFactors, knn, interpolation_setup
#The following exch_UV differs from normal exchange; incl. exch_UV_N
export exch_UV
export nansum, nanmean, nanmax, nanmin
export demo, land_mask, edge_path
export Integration

#export InnerArray, OuterArray
#export smooth, mask
#export GridOfOnes 
#export GridAddWS!

export ScalarPotential, VectorPotential
export UVtoUEVN, UVtoTransport, UVtoTransport!
export gradient, curl, convergence
export LatitudeCircles, Transect, ThroughFlow
export StereographicProjection, isosurface

#export location_is_out, NeighborTileIndices_dpdo, NeighborTileIndices_cs, RelocationFunctions_cs
#export update_location_cs!, update_location_llc!, update_location_dpdo!

"""
    Dataset(d::String)

- get folder name for dataset `d`. 
- download dataset if needed.

"""
function Dataset(d::String; do_read=true, verbose=false)
    try
        da=mydatadep(d)
        fi=mydatafile(d)
        file=joinpath(da,fi)

        if do_read==true
            if occursin(".geojson",file)&&file[end-7:end]==".geojson"
                verbose ? println("read_json") : nothing
                read_json(file)
            elseif occursin(".json",file)&&file[end-4:end]==".json"
                verbose ? println("read_json") : nothing
                read_json(file)
            elseif occursin(".shp",file)&&file[end-3:end]==".shp"
                verbose ? println("read_shp") : nothing
                read_shp(file)
            else
                verbose ? println("unknown read method") : nothing
                file
            end
        else
            file
        end

    catch
        d
    end
end

function mydatafile(nam="countries_shp1")
    if nam=="countries_shp1"
        "ne_110m_admin_0_countries.shp"
    elseif nam=="countries_geojson1"
        "countries.geojson"
    elseif nam=="basemap_jpg1"
        "Blue_Marble_Next_Generation_%2B_topography_%2B_bathymetry.jpg"
    elseif nam=="interp_halfdeg"
        "interp_coeffs_halfdeg.jld2"
    elseif nam=="oceans_geojson1"
        "ocean_basins_res1000_20251109_GF.json"
    else
        error("unknown data dependency")
    end
end

function mydatadep end

MA_datadep=mydatadep
export MA_datadep

function plot_examples end; export plot_examples
examples_plot=plot_examples

end # module
