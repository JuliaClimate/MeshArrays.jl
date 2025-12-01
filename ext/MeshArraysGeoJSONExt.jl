
############################################################
#                                                          #
#   Reading & Converting Polygons from GeoJSON             #
#                                                          #
############################################################

module MeshArraysGeoJSONExt

    using GeoJSON
    import MeshArrays: read_json

    """
        read_json(fil; native=false)

    - (default) call `GeoJSON.read` then convert to `GI.coordinates`.
    - if `native` then skip conversion, just call `GeoJSON.read`.

    ```
    pol=MeshArrays.Dataset("countries_geojson1",do_read=false)
    pol=MeshArrays.read_json(pol,native=true)
    ```
    """
    function read_json(fil; native=false)
        if native
            GeoJSON.read(fil)
        else
            jsonbytes = read(fil)
            geoms = GeoJSON.read(jsonbytes)
            #GeoJSON.GI.coordinates(geoms)
            [GeoJSON.GI.coordinates(a) for a in geoms]
        end
    end

end #module MeshArraysGeoJSONExt