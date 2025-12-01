

"""
    Dataset(d::String; do_read=true, verbose=false)

- get folder name for dataset `d`. 
- download dataset if needed.
- read content if `do_read`

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
            artifact"GRID_LLC90"
            hash = artifact_hash("GRID_LLC90", artifact_toml)
            joinpath(artifact_path(hash)*"/","GRID_LLC90-1.1/")
        elseif d=="GRID_CS32"
            artifact"GRID_CS32"
            hash = artifact_hash("GRID_CS32", artifact_toml)
            joinpath(artifact_path(hash)*"/","GRID_CS32-1.1/")
        elseif d=="GRID_LLC270"
            artifact"GRID_LLC270"
            hash = artifact_hash("GRID_LLC270", artifact_toml)
            joinpath(artifact_path(hash)*"/","GRID_LLC270-1.0.0/")
        elseif d=="GRID_LL360"
            artifact"GRID_LL360"
            hash = artifact_hash("GRID_LL360", artifact_toml)
            joinpath(artifact_path(hash)*"/","GRID_LL360-1.0/")
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

"""
	read_polygons(file="countries.geojson")

Call `read_json` or `read_shp` depending on file extension.

"""
function read_polygons(file::String)
	if !isfile(file)
		error("file not found ($file)")
	elseif occursin(".geojson",file)&&file[end-7:end]==".geojson"
		read_json(file)
	elseif occursin(".shp",file)&&file[end-3:end]==".shp"
		read_shp(file)
	else
		error("unknown file extension ($file)")
	end
end

