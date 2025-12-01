

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
        if d=="GRID_LLC90"
            MeshArrays.GRID_LLC90                
        elseif d=="GRID_CS32"
            MeshArrays.GRID_CS32
        elseif d=="GRID_LLC270"
            MeshArrays.GRID_LLC270
        elseif d=="GRID_LL360"
            MeshArrays.GRID_LL360
        else
            d
        end
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



