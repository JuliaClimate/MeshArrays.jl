
############################################################
#                                                          #
#   Reading & Converting Polygons from GeoJSON             #
#                                                          #
############################################################

module MeshArraysGeoJSONExt

    using GeoJSON
    import MeshArrays: read_json, to_polyarray, polyarray, to_Polygon
    import Base: write

    """
        read_json(fil; format=1)

    Call `GeoJSON.read` and return `polyarray` (default)

    - format 1 (default): return `polyarray`
    - format 0.1 : convert to `GI.coordinates`.
    - format 0.2 : return FeatureCollection (vector of name,geom pairs)

    ```
    import MeshArrays, DataDeps, GeoJSON
    pol=MeshArrays.Dataset("oceans_geojson1")
    ```
    """
    function read_json(fil; format=1)
        if format==1
            to_polyarray(GeoJSON.read(fil))
        elseif format==1.2
            GeoJSON.read(fil)
        elseif format==0.1
            jsonbytes = read(fil)
            geoms = GeoJSON.read(jsonbytes)
            #GeoJSON.GI.coordinates(geoms)
            [GeoJSON.GI.coordinates(a) for a in geoms]
        else
            error("unknown format")
        end
    end

    """
        write(pa::polyarray,file=tempname()*".json")

    ```
    import MeshArrays, DataDeps, GeoJSON
    pol=MeshArrays.Dataset("oceans_geojson1")
    write(pol,tempname()*".json")
    ```
    """
    function write(pa::polyarray,file=tempname()*".json")
        pol_P,nams=to_Polygon(pa)
        df=(geom=pol_P,id=nams,name=nams) #warning : use geom, not geometry, for name
        st=GeoJSON._lower(df,geometrycolumn=:geom)
        GeoJSON.JSON3.write(file, st)
    end

end #module MeshArraysGeoJSONExt