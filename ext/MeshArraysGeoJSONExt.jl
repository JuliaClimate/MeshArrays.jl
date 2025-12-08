
############################################################
#                                                          #
#   Reading & Converting Polygons from GeoJSON             #
#                                                          #
############################################################

module MeshArraysGeoJSONExt

    using GeoJSON
    import MeshArrays: read_json, to_polyarray, polyarray, to_Polygon
    import Base: write

    function read_json(fil; format=:polyarray)
        if format==:polyarray
            to_polyarray(GeoJSON.read(fil))
        elseif format==:GeoJSON
            GeoJSON.read(fil)
        elseif format==:coord
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